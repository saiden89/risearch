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

const USE_SINGLETON_FASTPATH: bool = true;
const LINEAR_PARTITION_CUTOFF: usize = 1024;

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
                0, // depth
                0, // match streak for suffix
                0, // mismatch count
            );
        } else {
            recurse::<_, false>(
                &mut ctx,
                self.q_sa_start,
                self.q_sa_len,
                0,
                self.t_sa_len,
                0, // depth
                0, // match streak for suffix
                0, // mismatch count
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

fn recurse<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut RecurseCtx<'_, F>,
    ql: usize,
    qr: usize,
    sl: usize,
    sr: usize,
    depth: usize,
    msm: usize, // consecutive match streak (suffix matches)
    mc: usize,  // mismatch count
) {
    // Emit match if within valid length range
    if depth >= ctx.min_len
        && depth <= ctx.max_len
        && (mc == 0 || (msm >= ctx.min_suffix && msm < depth))
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

    // Prune: can't accumulate enough suffix matches
    if mc > 0 && ctx.min_suffix > 0 {
        let max_possible = msm + (ctx.max_len - depth);
        if max_possible < ctx.min_suffix {
            return;
        }
    }

    if USE_SINGLETON_FASTPATH {
        if qr - ql == 1 && sr - sl == 1 {
            recurse_singleton::<F, WOBBLE>(ctx, ql, sl, depth, msm, mc);
            return;
        }
        if qr - ql == 1 {
            recurse_q_singleton::<F, WOBBLE>(ctx, ql, sl, sr, depth, msm, mc);
            return;
        }
        if sr - sl == 1 {
            recurse_s_singleton::<F, WOBBLE>(ctx, ql, qr, sl, depth, msm, mc);
            return;
        }
    }

    // Partition both SA intervals by base at current depth.
    // C-style layout from `sa_search_interval("acgnu")`:
    // [A_start, C_start, G_start, N_start, U_start, end]
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
        && mc < ctx.max_mm
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

    macro_rules! rec {
        ($q_lo:expr, $q_hi:expr, $s_lo:expr, $s_hi:expr, $next_msm:expr, $next_mc:expr) => {
            if $q_lo < $q_hi && $s_lo < $s_hi {
                recurse::<F, WOBBLE>(ctx, $q_lo, $q_hi, $s_lo, $s_hi, d1, $next_msm, $next_mc);
            }
        };
    }

    // Match branches in C query-class order: A, C, G, U.
    rec!(qa_lo, qa_hi, su_lo, su_hi, msm + 1, mc); // A-U
    rec!(qc_lo, qc_hi, sg_lo, sg_hi, msm + 1, mc); // C-G
    if qg_lo < qg_hi {
        rec!(qg_lo, qg_hi, sc_lo, sc_hi, msm + 1, mc); // G-C
        if WOBBLE {
            rec!(qg_lo, qg_hi, su_lo, su_hi, msm + 1, mc); // G-U wobble
        }
    }
    if qu_lo < qu_hi {
        rec!(qu_lo, qu_hi, sa_lo, sa_hi, msm + 1, mc); // U-A
        if WOBBLE {
            rec!(qu_lo, qu_hi, sg_lo, sg_hi, msm + 1, mc); // U-G wobble
        }
    }

    // Mismatch branches.
    if !can_mm {
        return;
    }

    // q = A (matches only U)
    if qa_lo < qa_hi {
        rec!(qa_lo, qa_hi, sa_lo, sa_hi, 0, mc + 1);
        rec!(qa_lo, qa_hi, sc_lo, sc_hi, 0, mc + 1);
        rec!(qa_lo, qa_hi, sg_lo, sg_hi, 0, mc + 1);
    }

    // q = C (matches G)
    if qc_lo < qc_hi {
        rec!(qc_lo, qc_hi, sa_lo, sa_hi, 0, mc + 1);
        rec!(qc_lo, qc_hi, sc_lo, sc_hi, 0, mc + 1);
        rec!(qc_lo, qc_hi, su_lo, su_hi, 0, mc + 1);
    }

    // q = G (matches C and wobble U)
    if qg_lo < qg_hi {
        rec!(qg_lo, qg_hi, sa_lo, sa_hi, 0, mc + 1);
        rec!(qg_lo, qg_hi, sg_lo, sg_hi, 0, mc + 1);
        if !WOBBLE {
            rec!(qg_lo, qg_hi, su_lo, su_hi, 0, mc + 1);
        }
    }

    // q = U (matches A and wobble G)
    if qu_lo < qu_hi {
        if !WOBBLE {
            rec!(qu_lo, qu_hi, sg_lo, sg_hi, 0, mc + 1);
        }
        rec!(qu_lo, qu_hi, sc_lo, sc_hi, 0, mc + 1);
        rec!(qu_lo, qu_hi, su_lo, su_hi, 0, mc + 1);
    }
}

#[cfg(test)]
#[inline(always)]
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
