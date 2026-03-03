//! Parallel Suffix Array seed search — direct port of C's sa_parallel_match_neg.
//!
//! Three functions mirroring C exactly:
//! - `sa_search_left`: binary search for base boundary in SA
//! - `partition_interval`: partition SA interval into [A, C, G, N, U] ranges
//! - `recurse`: recursive parallel SA traversal with inline match/mismatch

use crate::config::SeedConfig;
use crate::types::{Base, Interval};

mod partition;
mod singleton;

use partition::partition_interval_into;
use singleton::{recurse_q_singleton, recurse_s_singleton, recurse_singleton};

// Raw Base discriminant values (from #[repr(u8)] Base enum), as returned by sa_char().
const BASE_A: u8 = Base::A as u8;
const BASE_G: u8 = Base::G as u8;
const BASE_C: u8 = Base::C as u8;
const BASE_U: u8 = Base::U as u8;

// Rank-ordered values for SA partitioning.
// Rank order: Gap(0) < A(1) < C(2) < G(3) < N(4) < U(5).
// Mirrors C's XSTRM + "acgnu" character class ordering.
const RANK_A: u8 = 1;
const RANK_C: u8 = 2;
const RANK_G: u8 = 3;
const RANK_N: u8 = 4;
const RANK_U: u8 = 5;

/// Map raw Base discriminant to sorted rank for SA partitioning.
#[inline(always)]
fn sa_char_rank(sa: &[u64], seq: &[Base], sa_idx: usize, offset: usize) -> u8 {
    match sa_char(sa, seq, sa_idx, offset) {
        BASE_A => RANK_A,
        BASE_C => RANK_C,
        BASE_G => RANK_G,
        BASE_U => RANK_U,
        5 => RANK_N, // Base::N — rarely hit, not worth a const
        _ => 0,      // Gap/sentinel
    }
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
    pub(crate) query_interval: Interval,
    pub(crate) target_interval: Interval,
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
        let mut ctx = RecurseCtx {
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

        if self.cfg.allows_wobble() {
            recurse::<_, true>(
                &mut ctx,
                self.q_sa_start,
                self.q_sa_len,
                0,
                self.t_sa_len,
                0,
                0,
                0,
            );
        } else {
            recurse::<_, false>(
                &mut ctx,
                self.q_sa_start,
                self.q_sa_len,
                0,
                self.t_sa_len,
                0,
                0,
                0,
            );
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
struct RecurseCtx<'a, F: FnMut(SeedMatch)> {
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

#[allow(clippy::too_many_arguments)]
fn recurse<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut RecurseCtx<'_, F>,
    ql: usize,
    qr: usize,
    sl: usize,
    sr: usize,
    depth: usize,
    match_streak: usize, // consecutive matches at tail (for suffix constraint)
    mm_count: usize,     // mismatches accumulated so far
) {
    // Emit match if within valid length range.
    // When mm_count > 0, also require: enough trailing matches (suffix constraint)
    // and match_streak < depth (ensures at least one mismatch actually occurred,
    // preventing pure-match seeds from being re-emitted on the mismatch path).
    if depth >= ctx.min_len
        && depth <= ctx.max_len
        && (mm_count == 0 || (match_streak >= ctx.min_suffix && match_streak < depth))
    {
        (ctx.on_match)(SeedMatch {
            query_interval: Interval::new(ql, qr),
            target_interval: Interval::new(sl, sr),
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
        recurse_singleton::<F, WOBBLE>(ctx, ql, sl, depth, match_streak, mm_count);
        return;
    }
    if qr - ql == 1 {
        recurse_q_singleton::<F, WOBBLE>(ctx, ql, sl, sr, depth, match_streak, mm_count);
        return;
    }
    if sr - sl == 1 {
        recurse_s_singleton::<F, WOBBLE>(ctx, ql, qr, sl, depth, match_streak, mm_count);
        return;
    }

    // Partition both SA intervals by base at current depth.
    // C-style layout from `sa_search_interval("acgnu")`:
    // [A_start, C_start, G_start, N_start, U_start, end]
    // N (rank 4) boundaries are computed but never matched — N can't form
    // any valid base pair, so its partition slots serve only as G's upper bound.
    let mut qint = [0usize; 6];
    let mut sint = [0usize; 6];
    partition_interval_into(ctx.q_sa, ctx.q_seq, ql, qr, depth, &mut qint);
    partition_interval_into(ctx.t_sa, ctx.t_seq, sl, sr, depth, &mut sint);

    // Match C's `sa_search_interval` early return: no acgnu class in either side.
    if qint[0] == qr || sint[0] == sr {
        return;
    }

    let d1 = depth + 1;

    // Can we introduce a mismatch at this position?
    let can_mm = ctx.max_mm > 0
        && mm_count < ctx.max_mm
        && d1 > ctx.min_prefix
        && ctx.max_len - d1 >= ctx.min_suffix;

    let (qa_lo, qa_hi) = (qint[0], qint[1]);
    let (qc_lo, qc_hi) = (qint[1], qint[2]);
    let (qg_lo, qg_hi) = (qint[2], qint[3]);
    let (qu_lo, qu_hi) = (qint[4], qint[5]);

    let (sa_lo, sa_hi) = (sint[0], sint[1]);
    let (sc_lo, sc_hi) = (sint[1], sint[2]);
    let (sg_lo, sg_hi) = (sint[2], sint[3]);
    let (su_lo, su_hi) = (sint[4], sint[5]);

    // Match branches in C query-class order: A, C, G, U.
    recurse_if_nonempty::<F, WOBBLE>(
        ctx,
        qa_lo,
        qa_hi,
        su_lo,
        su_hi,
        d1,
        match_streak + 1,
        mm_count,
    ); // A-U
    recurse_if_nonempty::<F, WOBBLE>(
        ctx,
        qc_lo,
        qc_hi,
        sg_lo,
        sg_hi,
        d1,
        match_streak + 1,
        mm_count,
    ); // C-G
    if qg_lo < qg_hi {
        recurse_if_nonempty::<F, WOBBLE>(
            ctx,
            qg_lo,
            qg_hi,
            sc_lo,
            sc_hi,
            d1,
            match_streak + 1,
            mm_count,
        ); // G-C
        if WOBBLE {
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                qg_lo,
                qg_hi,
                su_lo,
                su_hi,
                d1,
                match_streak + 1,
                mm_count,
            ); // G-U wobble
        }
    }
    if qu_lo < qu_hi {
        recurse_if_nonempty::<F, WOBBLE>(
            ctx,
            qu_lo,
            qu_hi,
            sa_lo,
            sa_hi,
            d1,
            match_streak + 1,
            mm_count,
        ); // U-A
        if WOBBLE {
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                qu_lo,
                qu_hi,
                sg_lo,
                sg_hi,
                d1,
                match_streak + 1,
                mm_count,
            ); // U-G wobble
        }
    }

    // Mismatch branches.
    if !can_mm {
        return;
    }

    // q = A (matches only U)
    if qa_lo < qa_hi {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qa_lo, qa_hi, sa_lo, sa_hi, d1, 0, mm_count + 1);
        recurse_if_nonempty::<F, WOBBLE>(ctx, qa_lo, qa_hi, sc_lo, sc_hi, d1, 0, mm_count + 1);
        recurse_if_nonempty::<F, WOBBLE>(ctx, qa_lo, qa_hi, sg_lo, sg_hi, d1, 0, mm_count + 1);
    }

    // q = C (matches G)
    if qc_lo < qc_hi {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qc_lo, qc_hi, sa_lo, sa_hi, d1, 0, mm_count + 1);
        recurse_if_nonempty::<F, WOBBLE>(ctx, qc_lo, qc_hi, sc_lo, sc_hi, d1, 0, mm_count + 1);
        recurse_if_nonempty::<F, WOBBLE>(ctx, qc_lo, qc_hi, su_lo, su_hi, d1, 0, mm_count + 1);
    }

    // q = G (matches C and wobble U)
    if qg_lo < qg_hi {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qg_lo, qg_hi, sa_lo, sa_hi, d1, 0, mm_count + 1);
        recurse_if_nonempty::<F, WOBBLE>(ctx, qg_lo, qg_hi, sg_lo, sg_hi, d1, 0, mm_count + 1);
        if !WOBBLE {
            recurse_if_nonempty::<F, WOBBLE>(ctx, qg_lo, qg_hi, su_lo, su_hi, d1, 0, mm_count + 1);
        }
    }

    // q = U (matches A and wobble G)
    if qu_lo < qu_hi {
        if !WOBBLE {
            recurse_if_nonempty::<F, WOBBLE>(ctx, qu_lo, qu_hi, sg_lo, sg_hi, d1, 0, mm_count + 1);
        }
        recurse_if_nonempty::<F, WOBBLE>(ctx, qu_lo, qu_hi, sc_lo, sc_hi, d1, 0, mm_count + 1);
        recurse_if_nonempty::<F, WOBBLE>(ctx, qu_lo, qu_hi, su_lo, su_hi, d1, 0, mm_count + 1);
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn recurse_if_nonempty<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut RecurseCtx<'_, F>,
    ql: usize,
    qr: usize,
    sl: usize,
    sr: usize,
    depth: usize,
    match_streak: usize,
    mm_count: usize,
) {
    if ql < qr && sl < sr {
        recurse::<F, WOBBLE>(ctx, ql, qr, sl, sr, depth, match_streak, mm_count);
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
