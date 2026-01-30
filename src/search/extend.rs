//! DP extension algorithm for seed-and-extend search.
//!
//! This module contains the core algorithmic logic for extending seed matches
//! using dynamic programming, including maximality checks and alignment building.

use log::trace;
use smallvec::SmallVec;

use crate::alignment::{Alignment, Pairing};
use crate::dp::{self, DpExtension};
use crate::dsm::{DsmModel, pair_mat};
use crate::seed::{SeedHit, build_seed_alignment};
use crate::seq::Sequence;
use crate::types::Base;

use super::{FilterReason, MAX_DP_EXT, SearchContext, SearchStage};

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

#[inline]
fn base_or_gap(seq: &Sequence, pos: usize) -> Base {
    if pos >= seq.len() {
        Base::Gap
    } else {
        seq[pos]
    }
}

#[inline]
fn left_base(seq: &Sequence, anchor: usize, offset: usize) -> Base {
    if offset > anchor {
        Base::Gap
    } else {
        seq[anchor - offset]
    }
}

// =============================================================================
// EXTENSION RESULT
// =============================================================================

/// Result of DP extension phase. Stores scores and optional traces.
/// Alignment is built lazily only when output format requires it.
#[derive(Debug, Clone)]
pub(super) struct ExtensionResult {
    pub score: f64,
    pub l_q: usize,
    pub l_t: usize,
    pub r_q: usize,
    pub r_t: usize,
    pub left_ext: Option<DpExtension>,
    pub right_ext: Option<DpExtension>,
    pub q_pos: usize,
    pub t_match_end: usize,
    pub seed_len: usize,
}

// =============================================================================
// CORE EXTENSION ALGORITHM
// =============================================================================

/// Extend a seed match using DP and return extension result.
///
/// Performs maximality checks to skip non-maximal seeds, then extends
/// in both directions using the DSM energy model.
///
/// # Complexity
/// - Time: O(max_extension²) for DP
/// - Space: O(max_extension) for trace storage
pub(super) fn extend_seed<M: DsmModel>(
    ctx: &mut SearchContext<'_, '_, M>,
    q_seq: &Sequence,
    t_seq: &Sequence,
    candidate: &SeedHit,
    interval: crate::registry::SeedInterval,
) -> Option<ExtensionResult> {
    let q_pos = candidate.query_pos;
    let t_pos = candidate.target_start;
    let len = candidate.len;
    let opts = &ctx.args.extend;

    // Use pre-computed interval bounds for maximality constraint
    let (interval_start, interval_end) = (interval.start, interval.end);

    let query = q_seq;
    let target = t_seq;
    let pair_mat = pair_mat(ctx.args.seed.allows_wobble());

    trace!(
        "{} ENTERING extend_seed q_pos={} t_pos={} len={} delta_g={}",
        SearchStage::Extend,
        q_pos,
        t_pos,
        len,
        opts.delta_g
    );

    // MAXIMALITY CHECK: Left
    if q_pos > interval_start && t_pos + len < t_seq.len() {
        let q_prev = query[q_pos - 1].idx();
        let t_next = target[t_pos + len].idx();
        let p_class = pair_mat[q_prev][t_next];

        trace!(
            "{} Left Check q_pos={} t_pos={} len={} q_prev={} t_next={} pair={}",
            SearchStage::Extend,
            q_pos,
            t_pos,
            len,
            q_prev,
            t_next,
            p_class
        );

        if p_class != 0 {
            if opts.no_max_prune {
                trace!(
                    "SEED: Non-maximal Left-Ext: q_pos={} t_pos={} len={} pair={} (not pruning due to --no-max-prune)",
                    q_pos, t_pos, len, p_class
                );
            } else {
                ctx.stats.record_filter(FilterReason::MaximalityLeft);
                trace!(
                    "{} FILTERED reason={:?} q_pos={} t_pos={} len={} pair={}",
                    SearchStage::Extend,
                    FilterReason::MaximalityLeft,
                    q_pos,
                    t_pos,
                    len,
                    p_class
                );
                return None;
            }
        }
    }

    // MAXIMALITY CHECK: Right
    if q_pos + len < interval_end && t_pos > 0 {
        let q_next = query[q_pos + len].idx();
        let t_prev = target[t_pos - 1].idx();
        let p_class = pair_mat[q_next][t_prev];

        trace!(
            "{} Right Check q_pos={} t_pos={} len={} q_next={} t_prev={} pair={}",
            SearchStage::Extend,
            q_pos,
            t_pos,
            len,
            q_next,
            t_prev,
            p_class
        );

        if p_class != 0 {
            if opts.no_max_prune {
                trace!(
                    "SEED: Non-maximal Right-Ext: q_pos={} t_pos={} len={} pair={} (not pruning due to --no-max-prune)",
                    q_pos, t_pos, len, p_class
                );
            } else {
                ctx.stats.record_filter(FilterReason::MaximalityRight);
                trace!(
                    "{} FILTERED reason={:?} q_pos={} t_pos={} len={} pair={}",
                    SearchStage::Extend,
                    FilterReason::MaximalityRight,
                    q_pos,
                    t_pos,
                    len,
                    p_class
                );
                return None;
            }
        }
    }

    trace!(
        "{} Accepted Seed: q_pos={} t_pos={} len={} delta_g={}",
        SearchStage::Extend,
        q_pos,
        t_pos,
        len,
        opts.delta_g
    );

    let t_match_end = t_pos + len - 1;
    let max_ext = opts.max_extension as usize;
    let safe_ext = max_ext.min(MAX_DP_EXT);

    // Check extendability
    let left_extendable = q_pos > 0 && t_pos + len < t_seq.len();
    let right_extendable = q_pos + len < q_seq.len() && t_pos > 0;
    let seed_extendable = left_extendable || right_extendable;
    let do_extension = max_ext > 0 && seed_extendable;

    if !do_extension {
        trace!(
            "{} SKIPPING EXTENSION: max_ext={} extendable={} (L={} R={})",
            SearchStage::Extend,
            max_ext,
            seed_extendable,
            left_extendable,
            right_extendable
        );
    }

    // Seed energy calculation
    let seed_energy_raw = M::seed_energy(query, target, q_pos, t_match_end, len);

    if log::log_enabled!(log::Level::Trace) {
        let mut seed_int_str = String::with_capacity(len);
        for k in 0..len {
            let qc = query[q_pos + k];
            let tc = target[t_match_end - k];
            seed_int_str.push(qc.pairing_class(tc));
        }
        trace!(
            "[SEED] energy_raw={} q_pos={} t_end={} len={} interaction={}",
            seed_energy_raw, q_pos, t_match_end, len, seed_int_str
        );
    }

    // Seed-only result if not extending
    if !do_extension {
        let term_5p = M::terminal_5p(query[q_pos], target[t_match_end]);
        let term_3p = M::terminal_3p(query[q_pos + len - 1], target[t_pos]);
        let total_raw = seed_energy_raw + term_5p + term_3p;
        let final_score = M::to_kcal(total_raw);

        return Some(ExtensionResult {
            score: final_score,
            l_q: 0,
            l_t: 0,
            r_q: 0,
            r_t: 0,
            left_ext: None,
            right_ext: None,
            q_pos,
            t_match_end,
            seed_len: len,
        });
    }

    // DP extension
    let left_res = ctx
        .extender
        .extend_left(query, target, q_pos, t_match_end, safe_ext);
    let right_res = ctx
        .extender
        .extend_right(query, target, q_pos + len - 1, t_pos, safe_ext);

    let final_score = M::to_kcal(seed_energy_raw + left_res.score + right_res.score);

    trace!(
        "{} score={:.2} (seed={:.2} L={} R={}) L_len={}/{} R_len={}/{}",
        SearchStage::Extend,
        final_score,
        seed_energy_raw as f64 / -100.0,
        left_res.score,
        right_res.score,
        left_res.q_len,
        left_res.t_len,
        right_res.q_len,
        right_res.t_len
    );

    Some(ExtensionResult {
        score: final_score,
        l_q: left_res.q_len,
        l_t: left_res.t_len,
        r_q: right_res.q_len,
        r_t: right_res.t_len,
        left_ext: Some(left_res),
        right_ext: Some(right_res),
        q_pos,
        t_match_end,
        seed_len: len,
    })
}

// =============================================================================
// ALIGNMENT BUILDING
// =============================================================================

/// Build alignment from extension result.
/// Called only when output format requires alignment strings.
pub(super) fn build_alignment_from_extension(
    ext: &ExtensionResult,
    q_seq: &Sequence,
    t_seq: &Sequence,
    t_pos: usize,
) -> Alignment {
    // Seed-only case
    if ext.left_ext.is_none() && ext.right_ext.is_none() {
        let seed_alignment =
            build_seed_alignment(q_seq, t_seq, ext.q_pos, ext.t_match_end, ext.seed_len);
        return Alignment::new(&[], &seed_alignment, &[]);
    }

    let left_res = ext.left_ext.as_ref().expect("left_ext should be Some");
    let right_res = ext.right_ext.as_ref().expect("right_ext should be Some");

    // Left trace
    let mut left_alignment: SmallVec<[Pairing; 64]> = SmallVec::new();
    let mut li = left_res.q_len;
    let mut lj = left_res.t_len;
    for &step in &left_res.trace {
        let q_base = left_base(q_seq, ext.q_pos, li);
        let t_base = base_or_gap(t_seq, ext.t_match_end + lj);
        let pairing = Pairing::from_dp_op(step, q_base, t_base);

        // Update positions based on operation type
        match step {
            dp::DpOp::Match | dp::DpOp::Stop => {
                li = li.saturating_sub(1);
                lj = lj.saturating_sub(1);
            }
            dp::DpOp::GapQ => li = li.saturating_sub(1),
            dp::DpOp::GapT => lj = lj.saturating_sub(1),
        }
        left_alignment.push(pairing);
    }

    // Seed
    let seed_alignment =
        build_seed_alignment(q_seq, t_seq, ext.q_pos, ext.t_match_end, ext.seed_len);

    // Right trace
    let mut right_alignment: SmallVec<[Pairing; 64]> = SmallVec::new();
    let (mut ri, mut rj) = (0, 0);
    let q_anchor = ext.q_pos + ext.seed_len - 1;
    for &step in right_res.trace.iter().rev() {
        // Update positions BEFORE getting bases (right extension)
        match step {
            dp::DpOp::Match | dp::DpOp::Stop => {
                ri += 1;
                rj += 1;
            }
            dp::DpOp::GapQ => ri += 1,
            dp::DpOp::GapT => rj += 1,
        }

        let q_base = base_or_gap(q_seq, q_anchor + ri);
        let t_base = left_base(t_seq, t_pos, rj);
        let pairing = Pairing::from_dp_op(step, q_base, t_base);
        right_alignment.push(pairing);
    }

    Alignment::new(&left_alignment, &seed_alignment, &right_alignment)
}
