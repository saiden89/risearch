//! Parallel Suffix Array Search for RNA Interaction Discovery
//!
//! This module implements a parallel SA traversal algorithm that searches both
//! query and target suffix arrays simultaneously. Following the C RIsearch2
//! algorithm, the target SA is built on the COMPLEMENT of the target sequence,
//! enabling exact character matching to find complementary RNA base pairs.
//!
//! # Algorithm Overview
//!
//! The C algorithm builds target SA on complement(target). This means:
//! - Query A matching target_comp A → original target has U → A-U pair ✓
//! - Query G matching target_comp A → original target has U → G-U wobble ✓
//!
//! ```text
//! Level 0:  [entire query SA] × [entire target_complement SA]
//!               ↓
//! Level k:  Partition both SA intervals by base (A,C,G,U,N)
//!           Recurse on:
//!             - Same character pairs (canonical base-pairing via complement)
//!             - Wobble pairs: query G with target_comp A, query U with target_comp C
//!               ↓
//! Level s:  When depth == seed_length, collect matches
//! ```

use crate::config::SeedConfig;
use crate::types::{Base, Interval};

/// Partitioned SA intervals for ACGU bases, stored in iteration order.
///
/// `ranges[0]` = A, `ranges[1]` = C, `ranges[2]` = G, `ranges[3]` = U.
/// Each range is a (start, end) pair of SA indices.
#[derive(Debug, Clone, Copy)]
struct BaseIntervals {
    ranges: [(usize, usize); 4],
}

impl Default for BaseIntervals {
    fn default() -> Self {
        BaseIntervals {
            ranges: [(0, 0); 4],
        }
    }
}

impl BaseIntervals {
    /// Get interval for a specific base (for tests)
    #[cfg(test)]
    fn get(&self, base: Base) -> Interval {
        let (s, e) = match base {
            Base::A => self.ranges[0],
            Base::C => self.ranges[1],
            Base::G => self.ranges[2],
            Base::U => self.ranges[3],
            _ => (0, 0),
        };
        Interval::new(s, e)
    }

    #[inline]
    fn any_non_empty(&self) -> bool {
        let r = &self.ranges;
        r[0].0 < r[0].1 || r[1].0 < r[1].1 || r[2].0 < r[2].1 || r[3].0 < r[3].1
    }
}

/// A seed match found by parallel SA search
#[derive(Debug, Clone)]
pub struct SeedMatch {
    /// Interval in query SA containing matching suffixes
    pub(crate) query_interval: Interval,
    /// Interval in target SA containing matching suffixes
    pub(crate) target_interval: Interval,
    /// Seed length for this match.
    pub(crate) seed_len: usize,
}

/// Invariant context for the search (mirrors C's use of globals).
/// Passed by reference to avoid recomputing or chasing pointers each call.
struct SearchCtx<'a> {
    q_sa: &'a [u32],
    q_seq: &'a [Base],
    t_sa: &'a [u32],
    t_seq: &'a [Base],
    min_len: usize,
    max_len: usize,
    max_mm: usize,
    min_prefix: usize,
    min_suffix: usize,
    allow_wobble: bool,
}

/// Compact stack frame for iterative DFS.
/// 20 bytes vs 304-byte stack frame + 96 bytes of register spill in recursive version.
#[derive(Clone, Copy)]
#[repr(C)]
struct Frame {
    ql: u32,
    qr: u32,
    sl: u32,
    sr: u32,
    depth: u8,
    msm: u8,
    mc: u8,
}

/// Parallel Suffix Array searcher
/// - Query SA is built on the query sequence as-is
/// - Target SA is built on the COMPLEMENT of the target sequence
/// - Same-character matching finds complementary base pairs
pub struct SeedSearcher<'a> {
    query_sa: &'a [u32],
    query_seq: &'a [Base],
    target_comp_sa: &'a [u32],
    target_comp_seq: &'a [Base],
    seed_config: &'a SeedConfig,
}

impl<'a> SeedSearcher<'a> {
    /// Create a new parallel SA searcher
    ///
    /// IMPORTANT: `target_comp_sa` and `target_comp_seq` should be built on the
    /// COMPLEMENT (not reverse complement) of the target sequence.
    pub fn new(
        query_sa: &'a [u32],
        query_seq: &'a [Base],
        target_comp_sa: &'a [u32],
        target_comp_seq: &'a [Base],
        seed_config: &'a SeedConfig,
    ) -> Self {
        Self {
            query_sa,
            query_seq,
            target_comp_sa,
            target_comp_seq,
            seed_config,
        }
    }

    /// Find seeds for a range of lengths in a single traversal.
    pub fn search_length_range(
        &self,
        min_len: usize,
        max_len: usize,
        results: &mut Vec<SeedMatch>,
    ) {
        let ctx = SearchCtx {
            q_sa: self.query_sa,
            q_seq: self.query_seq,
            t_sa: self.target_comp_sa,
            t_seq: self.target_comp_seq,
            min_len,
            max_len,
            max_mm: self.seed_config.mismatch.max_mismatches,
            min_prefix: self.seed_config.mismatch.min_prefix_matches,
            min_suffix: self.seed_config.mismatch.min_suffix_matches,
            allow_wobble: self.seed_config.allows_wobble(),
        };
        search_iterative(
            &ctx,
            self.query_sa.len(),
            self.target_comp_sa.len(),
            results,
        );
    }
}

/// Iterative parallel SA search using explicit stack.
///
/// Replaces the recursive version to eliminate function-call overhead:
/// each "call" is a 20-byte Frame push instead of a 304-byte stack frame
/// with 12 callee-saved register spills.
fn search_iterative(ctx: &SearchCtx, q_len: usize, t_len: usize, results: &mut Vec<SeedMatch>) {
    let mut stack = Vec::with_capacity(512);
    stack.push(Frame {
        ql: 0,
        qr: q_len as u32,
        sl: 0,
        sr: t_len as u32,
        depth: 0,
        msm: 0,
        mc: 0,
    });

    while let Some(f) = stack.pop() {
        let ql = f.ql as usize;
        let qr = f.qr as usize;
        let sl = f.sl as usize;
        let sr = f.sr as usize;
        let depth = f.depth as usize;
        let msm = f.msm as usize;
        let mc = f.mc as usize;

        // Record match if within length range and valid
        if depth >= ctx.min_len
            && depth <= ctx.max_len
            && (mc == 0 || (mc <= ctx.max_mm && msm >= ctx.min_suffix && msm < depth))
        {
            results.push(SeedMatch {
                query_interval: Interval::new(ql, qr),
                target_interval: Interval::new(sl, sr),
                seed_len: depth,
            });
        }

        if depth >= ctx.max_len {
            continue;
        }

        // Prune: can't accumulate enough suffix matches
        if mc > 0 && ctx.min_suffix > 0 {
            let max_possible = msm + (ctx.max_len - depth);
            if max_possible < ctx.min_suffix {
                continue;
            }
        }

        // Prune: no suffix long enough to reach min_len
        if depth < ctx.min_len
            && (!has_suffix_len_at_least(ctx.q_sa, ctx.q_seq, ql, qr, ctx.min_len)
                || !has_suffix_len_at_least(ctx.t_sa, ctx.t_seq, sl, sr, ctx.min_len))
        {
            continue;
        }

        // Partition both SA intervals by base at current depth
        let qb = partition_interval(ctx.q_sa, ctx.q_seq, ql, qr, depth);
        let sb = partition_interval(ctx.t_sa, ctx.t_seq, sl, sr, depth);

        if !qb.any_non_empty() || !sb.any_non_empty() {
            continue;
        }

        let d1 = (depth + 1) as u8;
        let mc_u8 = mc as u8;
        let msm_p1 = (msm + 1) as u8;
        let mc_p1 = (mc + 1) as u8;

        let can_mm = ctx.max_mm > 0
            && mc < ctx.max_mm
            && (depth + 1) > ctx.min_prefix
            && msm < ctx.max_len
            && ctx.max_len - (depth + 1) >= ctx.min_suffix;

        // Push children in reverse order so forward-order is processed first (LIFO)
        for qi in (0..4usize).rev() {
            let (q_lo, q_hi) = qb.ranges[qi];
            if q_lo >= q_hi {
                continue;
            }
            let q_lo = q_lo as u32;
            let q_hi = q_hi as u32;

            for si in (0..4usize).rev() {
                let (s_lo, s_hi) = sb.ranges[si];
                if s_lo >= s_hi {
                    continue;
                }

                if qi == si || (ctx.allow_wobble && si == qi + 2) {
                    stack.push(Frame {
                        ql: q_lo,
                        qr: q_hi,
                        sl: s_lo as u32,
                        sr: s_hi as u32,
                        depth: d1,
                        msm: msm_p1,
                        mc: mc_u8,
                    });
                } else if can_mm {
                    stack.push(Frame {
                        ql: q_lo,
                        qr: q_hi,
                        sl: s_lo as u32,
                        sr: s_hi as u32,
                        depth: d1,
                        msm: 0,
                        mc: mc_p1,
                    });
                }
            }
        }
    }
}

/// Threshold below which linear scan beats 5 binary searches.
/// At n=32: linear ≈ 32×104 = 3328 cycles vs binary ≈ 5×5×200 = 5000 cycles.
const SCAN_THRESHOLD: usize = 32;

/// Partition an SA interval by base at given offset.
///
/// Uses linear scan for small intervals (sequential SA reads + fewer total
/// memory accesses) and binary search for large intervals.
#[inline]
fn partition_interval(
    sa: &[u32],
    seq: &[Base],
    start: usize,
    end: usize,
    offset: usize,
) -> BaseIntervals {
    let n = end - start;
    if n == 0 {
        return BaseIntervals::default();
    }

    if n <= SCAN_THRESHOLD {
        return partition_scan(sa, seq, start, end, offset);
    }

    partition_by_base(sa, seq, start, end, offset)
}

/// Linear scan partition for small SA intervals.
///
/// SA entries are sorted by suffix, so bases at `offset` appear in discriminant
/// order: Gap(0) < A(1) < G(2) < C(3) < U(4) < N(5). We scan once and record
/// the index where each new base first appears.
#[inline]
fn partition_scan(
    sa: &[u32],
    seq: &[Base],
    start: usize,
    end: usize,
    offset: usize,
) -> BaseIntervals {
    let seq_ptr = seq.as_ptr();

    // boundaries[b] = first index where base discriminant >= b+1
    // Index 0 = first A, 1 = first G, 2 = first C, 3 = first U, 4 = first N
    let mut boundaries = [end; 5];
    let mut last_disc: u8 = 0;

    for i in start..end {
        // SAFETY: same sentinel guarantee as partition_by_base
        let disc = unsafe { *seq_ptr.add(sa[i] as usize + offset) } as u8;
        if disc != last_disc {
            // Fill all boundary slots from last_disc+1 through disc
            for b in (last_disc as usize)..disc as usize {
                if boundaries[b] == end {
                    boundaries[b] = i;
                }
            }
            last_disc = disc;
        }
    }

    let [a_start, g_start, c_start, u_start, n_start] = boundaries;
    BaseIntervals {
        ranges: [
            (a_start, g_start), // A
            (c_start, u_start), // C
            (g_start, c_start), // G
            (u_start, n_start), // U
        ],
    }
}

/// Check if any suffix in the SA interval has at least `min_len` real bases.
///
/// Sequences include a trailing sentinel Gap byte, so the original (pre-sentinel)
/// length is `seq.len() - 1`. Scans from both ends with early exit; terminates
/// in 0-2 steps for typical genomic workloads.
#[inline]
fn has_suffix_len_at_least(
    sa: &[u32],
    seq: &[Base],
    start: usize,
    end: usize,
    min_len: usize,
) -> bool {
    if min_len == 0 {
        return true;
    }
    let original_len = seq.len() - 1; // exclude sentinel
    let offset = min_len - 1;

    let mut lo = start;
    while lo < end {
        if (sa[lo] as usize) + offset < original_len {
            return true;
        }
        lo += 1;
    }
    false
}

/// Partition a valid SA range by base, returning ranges in [A, C, G, U] order.
///
/// Uses unsafe pointer access (bounds guaranteed by sentinel padding).
/// Base discriminant ordering in SA: Gap(0) < A(1) < G(2) < C(3) < U(4) < N(5)
#[inline]
fn partition_by_base(
    sa: &[u32],
    seq: &[Base],
    valid_start: usize,
    valid_end: usize,
    offset: usize,
) -> BaseIntervals {
    let sa_slice = &sa[valid_start..valid_end];
    let seq_ptr = seq.as_ptr();

    // SAFETY: Sentinel Gap bytes appended after each sequence guarantee that
    // sa[i] + offset < seq.len() for all SA entries, since short suffixes
    // read the sentinel at seq[len] which is within bounds.
    // Base discriminant ordering: Gap(0) < A(1) < G(2) < C(3) < U(4) < N(5)
    let a_start = valid_start
        + sa_slice.partition_point(|&idx| unsafe { *seq_ptr.add(idx as usize + offset) < Base::A });
    let g_start = valid_start
        + sa_slice.partition_point(|&idx| unsafe { *seq_ptr.add(idx as usize + offset) < Base::G });
    let c_start = valid_start
        + sa_slice.partition_point(|&idx| unsafe { *seq_ptr.add(idx as usize + offset) < Base::C });
    let u_start = valid_start
        + sa_slice.partition_point(|&idx| unsafe { *seq_ptr.add(idx as usize + offset) < Base::U });
    let n_start = valid_start
        + sa_slice.partition_point(|&idx| unsafe { *seq_ptr.add(idx as usize + offset) < Base::N });

    // Ranges in [A, C, G, U] iteration order (matching BASES constant)
    BaseIntervals {
        ranges: [
            (a_start, g_start), // A: [a_start, g_start)
            (c_start, u_start), // C: [c_start, u_start)
            (g_start, c_start), // G: [g_start, c_start)
            (u_start, n_start), // U: [u_start, n_start)
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{MismatchSpec, SeedConfig, SeedSpec};
    use crate::index::sa::SuffixArray;
    use crate::seq::Sequence;

    /// Create a default SeedConfig for testing with specified wobble policy
    fn test_seed_config(allow_wobble: bool) -> SeedConfig {
        SeedConfig::with_wobble(SeedSpec::LengthOnly(6), MismatchSpec::exact(), allow_wobble)
    }

    #[test]
    fn test_base_complement() {
        use crate::types::Base;
        let seq = [Base::A, Base::C, Base::G, Base::U];
        let comp: Vec<_> = seq.iter().map(|b| b.complement()).collect();
        assert_eq!(comp, vec![Base::U, Base::G, Base::C, Base::A]);

        let seq2 = [Base::A, Base::A, Base::A, Base::A];
        let comp2: Vec<_> = seq2.iter().map(|b| b.complement()).collect();
        assert_eq!(comp2, vec![Base::U, Base::U, Base::U, Base::U]);
    }

    #[test]
    fn test_partition_basic() {
        use crate::types::Base;
        let seq_bases = vec![Base::A, Base::C, Base::G, Base::U];
        let seq = Sequence::from(seq_bases.clone());
        let sa = SuffixArray::try_from(&seq).expect("SA construction failed");
        // Append sentinel after SA construction
        let mut padded: Vec<Base> = seq.iter().copied().collect();
        padded.push(Base::Gap);
        let seq = Sequence::from(padded);
        let interval = Interval::new(0, sa.len());
        let parts = partition_interval(&sa, &seq, interval.start, interval.end, 0);

        // Each base should have exactly one entry
        assert_eq!(parts.get(Base::A).len(), 1);
        assert_eq!(parts.get(Base::C).len(), 1);
        assert_eq!(parts.get(Base::G).len(), 1);
        assert_eq!(parts.get(Base::U).len(), 1);
    }

    #[test]
    fn test_homopolymer_debug() {
        use crate::types::Base;
        // Debug test for the aaa/ttt case
        let query_bases = vec![Base::A, Base::A, Base::A];
        let target_bases = vec![Base::U, Base::U, Base::U];

        // Step 1: Check complement
        let target_comp_bases: Vec<_> = target_bases.iter().map(|b| b.complement()).collect();
        eprintln!("Query: {:?}", query_bases);
        eprintln!("Target: {:?}", target_bases);
        eprintln!("Target complement: {:?}", target_comp_bases);
        assert_eq!(
            &target_comp_bases,
            &[Base::A, Base::A, Base::A],
            "complement(UUU) should be AAA"
        );

        let query = Sequence::from(query_bases);
        let _target = Sequence::from(target_bases);
        let target_comp = Sequence::from(target_comp_bases);

        // Step 2: Build SAs on unpadded sequences, then append sentinel
        let q_sa = SuffixArray::try_from(&query).expect("query SA construction failed");
        let t_sa = SuffixArray::try_from(&target_comp).expect("target SA construction failed");
        eprintln!("Query SA: {:?}", q_sa);
        eprintln!("Target comp SA: {:?}", t_sa);

        let mut q_padded: Vec<Base> = query.iter().copied().collect();
        q_padded.push(Base::Gap);
        let q_seq = Sequence::from(q_padded);
        let mut t_padded: Vec<Base> = target_comp.iter().copied().collect();
        t_padded.push(Base::Gap);
        let t_seq = Sequence::from(t_padded);

        // Step 3: Create searcher and find seeds
        let seed_args = test_seed_config(false);
        let searcher = SeedSearcher::new(&q_sa, &q_seq, &t_sa, &t_seq, &seed_args);
        let mut matches = Vec::new();
        searcher.search_length_range(3, 3, &mut matches);
        eprintln!("Raw matches: {:?}", matches);

        // The key assertion
        assert!(
            !matches.is_empty(),
            "Should find at least one 3bp match for aaa vs aaa"
        );
    }
}
