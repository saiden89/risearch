//! Parallel Suffix Array seed search.
//!
//! Recursive traversal of two suffix arrays with base-pairing constraints.
//! Singleton fast-paths avoid partition overhead when one or both SA intervals
//! have a single entry.

use crate::config::SeedConfig;
use crate::types::Base;
use std::ops::Range;

mod partition;
mod singleton;

use partition::partition_interval_into;
use singleton::{recurse_q_singleton, recurse_s_singleton, recurse_singleton};

// Raw Base discriminant values (from #[repr(u8)] Base enum).
// Enum order: Gap(0) < A(1) < C(2) < G(3) < N(4) < U(5).
const BASE_A: u8 = Base::A as u8; // 1
const BASE_C: u8 = Base::C as u8; // 2
const BASE_G: u8 = Base::G as u8; // 3
const BASE_U: u8 = Base::U as u8; // 5

/// The four matchable RNA bases indexed alongside SLOTS.
const BASES: [Base; 4] = [Base::A, Base::C, Base::G, Base::U];

/// Partition slot for each base in BASES.
/// Partition layout: `[A=0, C=1, G=2, N=3, U=4]`; slot range = `p[slot]..p[slot+1]`.
const SLOTS: [usize; 4] = [0, 1, 2, 4];

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

/// True for the four matchable bases (A, C, G, U); false for N and Gap.
#[inline(always)]
fn is_valid_base(c: u8) -> bool {
    matches!(c, BASE_A | BASE_C | BASE_G | BASE_U)
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
/// Seed matching uses RNA pairing classes against this complemented target alphabet.
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
// RECURSIVE SEARCH
// =============================================================================

/// Shared context threaded through the recursive search.
///
/// In complement-transformed target space, canonical seed matches are:
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

    // Prune: can't accumulate enough suffix matches in remaining depth.
    if mm_count > 0 && ctx.min_suffix > 0 {
        let max_possible = match_streak + (ctx.max_len - depth);
        if max_possible < ctx.min_suffix {
            return;
        }
    }

    // Singleton fast paths: skip partition overhead when one or both intervals
    // have a single entry. Emission at this depth was already handled above.
    if qr - ql == 1 && sr - sl == 1 {
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
    let mut qi = [0usize; 6];
    let mut si = [0usize; 6];
    partition_interval_into(ctx.q_sa, ctx.q_seq, ql, qr, depth, &mut qi);
    partition_interval_into(ctx.t_sa, ctx.t_seq, sl, sr, depth, &mut si);

    // No ACGU bases on either side — nothing to pair.
    if qi[0] == qr || si[0] == sr {
        return;
    }

    let d1 = depth + 1;
    let ms = match_streak + 1;
    let can_mm = ctx.max_mm > 0
        && mm_count < ctx.max_mm
        && d1 > ctx.min_prefix
        && match_streak < ctx.min_len
        && ctx.max_len - d1 >= ctx.min_suffix;

    // Pair dispatch: iterate all (query_base, target_base) combinations.
    // LLVM unrolls and constant-folds pair_type().is_match() for each (i,j).
    for i in 0..4 {
        let qs = SLOTS[i];
        if qi[qs] >= qi[qs + 1] {
            continue;
        }

        for j in 0..4 {
            let ts = SLOTS[j];
            if si[ts] >= si[ts + 1] {
                continue;
            }

            if BASES[i].pair_type(BASES[j]).is_match(WOBBLE) {
                recurse::<F, WOBBLE>(ctx, qi[qs]..qi[qs + 1], si[ts]..si[ts + 1], d1, ms, mm_count);
            } else if can_mm {
                recurse::<F, WOBBLE>(
                    ctx,
                    qi[qs]..qi[qs + 1],
                    si[ts]..si[ts + 1],
                    d1,
                    0,
                    mm_count + 1,
                );
            }
        }
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
