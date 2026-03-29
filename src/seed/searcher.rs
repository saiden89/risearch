//! Parallel Suffix Array seed search — direct port of C's sa_parallel_match_neg.
//!
//! Three functions mirroring C exactly:
//! - `sa_search_left`: binary search for base boundary in SA
//! - `partition_interval`: partition SA interval into [A, C, G, N, U] ranges
//! - `recurse`: recursive parallel SA traversal with inline match/mismatch

use crate::config::SeedConfig;
use crate::types::Base;
use std::ops::Range;

mod partition;
mod singleton;

use partition::partition_interval_into;
use singleton::{recurse_q_singleton, recurse_s_singleton, recurse_singleton};

// Raw Base discriminant values (from #[repr(u8)] Base enum).
// Since the enum order is Gap < A < C < G < N < U, these also act as the SA sort ranks.
const BASE_A: u8 = Base::A as u8; // 1
const BASE_C: u8 = Base::C as u8; // 2
const BASE_G: u8 = Base::G as u8; // 3
const BASE_U: u8 = Base::U as u8; // 5

/// Map a raw base u8 to its partition slot index. Returns None for N and Gap.
///
/// Partition layout: A=0, C=1, G=2, N=3 (unused), U=4
/// Slot i spans `int[i]..int[i+1]` in the partition output array.
#[inline(always)]
fn base_to_slot(base: u8) -> Option<usize> {
    match base {
        b if b == BASE_A => Some(0),
        b if b == BASE_C => Some(1),
        b if b == BASE_G => Some(2),
        b if b == BASE_U => Some(4),
        _ => None,
    }
}

/// Returns true if (q_slot, t_slot) is a valid RNA match pair.
#[inline(always)]
fn is_match_slot<const WOBBLE: bool>(q_slot: usize, t_slot: usize) -> bool {
    matches!((q_slot, t_slot), (0, 4) | (1, 2) | (2, 1) | (4, 0))
        || (WOBBLE && matches!((q_slot, t_slot), (2, 4) | (4, 2)))
}

/// Read the base discriminant at `sa[sa_idx].pos + offset` from the sequence.
///
/// # Safety
/// SA_CHAR_PADDING sentinel entries guarantee `suffix_pos + offset` is within
/// bounds for any valid SA index and depth up to max_len.
#[inline(always)]
fn sa_char(sa: &[u64], seq: &[Base], sa_idx: usize, offset: usize) -> u8 {
    let suffix_pos = sa_suffix_pos(sa, sa_idx);
    unsafe { *seq.get_unchecked(suffix_pos + offset) as u8 }
}

/// Extract the suffix position from a position-only SA entry.
#[inline(always)]
fn sa_suffix_pos(sa: &[u64], sa_idx: usize) -> usize {
    unsafe { *sa.get_unchecked(sa_idx) as usize }
}

// =============================================================================
// SEED MATCH
// =============================================================================

/// A seed match found by parallel SA search.
#[derive(Debug, Clone)]
pub struct SeedMatch {
    pub(crate) query_interval: Range<usize>,
    pub(crate) target_interval: Range<usize>,
    pub(crate) seed_len: usize,
}

// =============================================================================
// SEED SEARCHER
// =============================================================================

/// Parallel SA searcher.
///
/// Query SA is built on the query as-is.
/// Target SA is built on the COMPLEMENT of the target.
/// To mirror C parity, seed matching uses RNA pairing classes against this
/// complemented target alphabet.
pub struct SeedSearcher<'a> {
    q_sa: &'a [u64],
    q_seq: &'a [Base],
    q_sa_start: usize,
    q_sa_len: usize,
    t_sa: &'a [u64],
    t_seq: &'a [Base],
    t_sa_len: usize,
    cfg: &'a SeedConfig,
}

impl<'a> SeedSearcher<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        q_sa: &'a [u64],
        q_seq: &'a [Base],
        q_sa_start: usize,
        q_sa_len: usize,
        t_sa: &'a [u64],
        t_seq: &'a [Base],
        t_sa_len: usize,
        cfg: &'a SeedConfig,
    ) -> Self {
        Self {
            q_sa,
            q_seq,
            q_sa_start,
            q_sa_len,
            t_sa,
            t_seq,
            t_sa_len,
            cfg,
        }
    }

    /// Stream seed matches for a length range.
    pub fn for_each_length_range<F>(&self, min_len: usize, max_len: usize, mut on_match: F)
    where
        F: FnMut(SeedMatch),
    {
        let mut ctx = SeedingContext {
            q_sa: self.q_sa,
            q_seq: self.q_seq,
            t_sa: self.t_sa,
            t_seq: self.t_seq,
            min_len,
            max_len,
            max_mm: self.cfg.mismatch.max_mismatches,
            min_prefix: self.cfg.mismatch.min_prefix_matches,
            min_suffix: self.cfg.mismatch.min_suffix_matches,
            on_match: &mut on_match,
        };

        let q = self.q_sa_start..self.q_sa_len;
        let s = 0..self.t_sa_len;
        if self.cfg.seed_wobble {
            recurse::<_, true>(&mut ctx, q, s, 0, 0, 0);
        } else {
            recurse::<_, false>(&mut ctx, q, s, 0, 0, 0);
        }
    }

    /// Collect all matches into a Vec.
    pub fn search_length_range(
        &self,
        min_len: usize,
        max_len: usize,
        results: &mut Vec<SeedMatch>,
    ) {
        self.for_each_length_range(min_len, max_len, |m| results.push(m));
    }
}

// =============================================================================
// RECURSIVE SEARCH — mirrors C's sa_parallel_match_neg
// =============================================================================

/// Recursive parallel SA search.
///
/// In complement-transformed target space, C-equivalent seed matches are:
/// - canonical: A↔U, G↔C, C↔G, U↔A
/// - wobble (optional): G↔U, U↔G
/// - everything else is mismatch
struct SeedingContext<'a, F: FnMut(SeedMatch)> {
    q_sa: &'a [u64],
    q_seq: &'a [Base],
    t_sa: &'a [u64],
    t_seq: &'a [Base],
    min_len: usize,
    max_len: usize,
    max_mm: usize,
    min_prefix: usize,
    min_suffix: usize,
    on_match: &'a mut F,
}

fn recurse<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut SeedingContext<'_, F>,
    q: Range<usize>,
    s: Range<usize>,
    depth: usize,
    match_streak: usize,
    mm_count: usize,
) {
    let (ql, qr) = (q.start, q.end);
    let (sl, sr) = (s.start, s.end);

    // Emit match if within valid length range.
    // When mm_count > 0, mirror C's `last_match_count < depth` where `depth`
    // is the minimum seed length (not current recursion depth). This keeps
    // mismatch seeds non-overlapping with pure-match seeds.
    if depth >= ctx.min_len
        && depth <= ctx.max_len
        && (mm_count == 0 || (match_streak >= ctx.min_suffix && match_streak < ctx.min_len))
    {
        (ctx.on_match)(SeedMatch {
            query_interval: ql..qr,
            target_interval: sl..sr,
            seed_len: depth,
        });
    }

    if depth >= ctx.max_len {
        return;
    }

    // Prune: can't accumulate enough suffix matches in remaining depth
    if mm_count > 0 && ctx.min_suffix > 0 {
        let max_possible = match_streak + (ctx.max_len - depth);
        if max_possible < ctx.min_suffix {
            return;
        }
    }

    // Singleton fast paths: when one or both SA intervals have a single entry,
    // skip partition overhead and compare characters directly.
    if qr - ql == 1 && sr - sl == 1 {
        // Emission at this depth was already handled above; singleton continues from here.
        recurse_singleton::<F, WOBBLE>(ctx, ql, sl, depth, match_streak, mm_count, false);
        return;
    }
    if qr - ql == 1 {
        recurse_q_singleton::<F, WOBBLE>(ctx, ql, sl..sr, depth, match_streak, mm_count);
        return;
    }
    if sr - sl == 1 {
        recurse_s_singleton::<F, WOBBLE>(ctx, ql..qr, sl, depth, match_streak, mm_count);
        return;
    }

    // Partition both SA intervals by base at current depth.
    // C-style layout from `sa_search_interval("acgnu")`:
    // [A_start, C_start, G_start, N_start, U_start, end]
    // Slot i spans int[i]..int[i+1]; N slot (3) is never matched.
    let mut qi = [0usize; 6];
    let mut si = [0usize; 6];
    partition_interval_into(ctx.q_sa, ctx.q_seq, ql, qr, depth, &mut qi);
    partition_interval_into(ctx.t_sa, ctx.t_seq, sl, sr, depth, &mut si);

    // Match C's `sa_search_interval` early return: no acgnu class in either side.
    if qi[0] == qr || si[0] == sr {
        return;
    }

    let d1 = depth + 1;

    // Can we introduce a mismatch at this position?
    let can_mm = ctx.max_mm > 0
        && mm_count < ctx.max_mm
        && d1 > ctx.min_prefix
        && match_streak < ctx.min_len
        && ctx.max_len - d1 >= ctx.min_suffix;

    let ms = match_streak + 1;

    // Match branches in C query-class order: A-U, C-G, G-C, U-A (+ wobble G-U, U-G).
    // G and U are outer-guarded to skip both canonical + wobble calls when the slot is empty.
    recurse_if_nonempty::<F, WOBBLE>(ctx, qi[0]..qi[1], si[4]..si[5], d1, ms, mm_count); // A-U
    recurse_if_nonempty::<F, WOBBLE>(ctx, qi[1]..qi[2], si[2]..si[3], d1, ms, mm_count); // C-G
    if qi[2] < qi[3] {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qi[2]..qi[3], si[1]..si[2], d1, ms, mm_count); // G-C
        if WOBBLE {
            recurse_if_nonempty::<F, WOBBLE>(ctx, qi[2]..qi[3], si[4]..si[5], d1, ms, mm_count); // G-U
        }
    }
    if qi[4] < qi[5] {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qi[4]..qi[5], si[0]..si[1], d1, ms, mm_count); // U-A
        if WOBBLE {
            recurse_if_nonempty::<F, WOBBLE>(ctx, qi[4]..qi[5], si[2]..si[3], d1, ms, mm_count); // U-G
        }
    }

    if !can_mm {
        return;
    }

    let mm1 = mm_count + 1;

    // Mismatch branches: all (q_slot, t_slot) pairs that are not canonical/wobble matches.
    if qi[0] < qi[1] {
        // q=A (matches U): mismatches A, C, G
        recurse_if_nonempty::<F, WOBBLE>(ctx, qi[0]..qi[1], si[0]..si[1], d1, 0, mm1);
        recurse_if_nonempty::<F, WOBBLE>(ctx, qi[0]..qi[1], si[1]..si[2], d1, 0, mm1);
        recurse_if_nonempty::<F, WOBBLE>(ctx, qi[0]..qi[1], si[2]..si[3], d1, 0, mm1);
    }
    if qi[1] < qi[2] {
        // q=C (matches G): mismatches A, C, U
        recurse_if_nonempty::<F, WOBBLE>(ctx, qi[1]..qi[2], si[0]..si[1], d1, 0, mm1);
        recurse_if_nonempty::<F, WOBBLE>(ctx, qi[1]..qi[2], si[1]..si[2], d1, 0, mm1);
        recurse_if_nonempty::<F, WOBBLE>(ctx, qi[1]..qi[2], si[4]..si[5], d1, 0, mm1);
    }
    if qi[2] < qi[3] {
        // q=G (matches C, wobble U): mismatches A, G, and U when no wobble
        recurse_if_nonempty::<F, WOBBLE>(ctx, qi[2]..qi[3], si[0]..si[1], d1, 0, mm1);
        recurse_if_nonempty::<F, WOBBLE>(ctx, qi[2]..qi[3], si[2]..si[3], d1, 0, mm1);
        if !WOBBLE {
            recurse_if_nonempty::<F, WOBBLE>(ctx, qi[2]..qi[3], si[4]..si[5], d1, 0, mm1);
        }
    }
    if qi[4] < qi[5] {
        // q=U (matches A, wobble G): mismatches C, U, and G when no wobble
        if !WOBBLE {
            recurse_if_nonempty::<F, WOBBLE>(ctx, qi[4]..qi[5], si[2]..si[3], d1, 0, mm1);
        }
        recurse_if_nonempty::<F, WOBBLE>(ctx, qi[4]..qi[5], si[1]..si[2], d1, 0, mm1);
        recurse_if_nonempty::<F, WOBBLE>(ctx, qi[4]..qi[5], si[4]..si[5], d1, 0, mm1);
    }
}

#[inline(always)]
fn recurse_if_nonempty<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut SeedingContext<'_, F>,
    q: Range<usize>,
    s: Range<usize>,
    depth: usize,
    match_streak: usize,
    mm_count: usize,
) {
    if q.start < q.end && s.start < s.end {
        recurse::<F, WOBBLE>(ctx, q, s, depth, match_streak, mm_count);
    }
}

#[cfg(test)]
fn partition_interval(
    sa: &[u64],
    seq: &[Base],
    start: usize,
    end: usize,
    offset: usize,
) -> [usize; 6] {
    partition::partition_interval(sa, seq, start, end, offset)
}

// =============================================================================
// TESTS
// =============================================================================

#[cfg(test)]
mod tests;
