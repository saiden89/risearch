use anyhow::Result;
use log::trace;
use rayon::prelude::*;

use crate::dp;
use crate::dsm::{PAIR_MAT, PAIR_MAT_NO_GU};
use crate::seed::{SeedCandidate, build_seed_alignment};
use crate::seq::Sequence;
use crate::types::{Alignment, Base, Pairing, SeedPairingMode, Strand};

use super::{
    FilterReason, MAX_DP_EXT, SaIndex, SearchContext, SearchHit, SearchStage, THREAD_EXTENDER,
};

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

/// Run search and return hits (used by parity tests and library callers).
/// Assumes the provided `SearchArgs` already contains a valid `SeedSpec` for each query.
pub fn run_search(
    queries: &[(String, Sequence)],
    index: &SaIndex<'_>,
    opts: &crate::config::SearchArgs,
) -> Result<Vec<SearchHit>> {
    let per_query_hits: Vec<Vec<SearchHit>> = queries
        .par_iter()
        .map(|(q_id, q_seq)| {
            THREAD_EXTENDER.with(|ext| {
                let mut extender = ext.borrow_mut();
                let mut ctx = SearchContext::with_extender(index, opts, &mut extender);
                let mut seeds = Vec::new();
                let mut matches = Vec::new();
                crate::seed::find_seeds(q_seq, ctx.index.index, &ctx.args.seed, &mut seeds, &mut matches);

                let query = super::Query::new(q_id, q_seq);
                seeds
                    .iter()
                    .filter_map(|candidate| process_candidate(query, candidate, &mut ctx))
                    .collect()
            })
        })
        .collect();

    Ok(per_query_hits.into_iter().flatten().collect())
}

/// Extend a single candidate and build a `SearchHit` if the energy passes filters.
fn process_candidate(
    query: super::Query<'_>,
    candidate: &SeedCandidate,
    ctx: &mut SearchContext<'_, '_>,
) -> Option<SearchHit> {
    let t_idx = candidate.target_idx;
    let t_seq = match candidate.strand {
        Strand::Reverse => ctx.index.get_sequence_rc(t_idx),
        Strand::Forward => ctx.index.get_sequence(t_idx),
    };

    let t_start_idx = candidate.target_start;
    let seed_len = candidate.len;
    if t_start_idx + seed_len > t_seq.len() {
        ctx.stats.record_filter(FilterReason::SeedOutOfBounds);
        return None;
    }

    let ext = extend_seed(ctx, query.seq, t_seq, candidate)?;

    let score = ext.score;
    if score > ctx.args.extend.delta_g {
        ctx.stats.record_filter(FilterReason::EnergyAboveThreshold);
        return None;
    }

    Some(build_hit(query, ctx, candidate, ext))
}

/// Finalize coordinates, bookkeeping, and formatting for an accepted hit.
fn build_hit(
    query: super::Query<'_>,
    ctx: &SearchContext<'_, '_>,
    candidate: &SeedCandidate,
    ext: ExtensionResult,
) -> SearchHit {
    let seed_len = candidate.len;
    let q_pos = candidate.query_pos;
    let t_start_idx = candidate.target_start;

    let final_q_start = q_pos.saturating_sub(ext.l_q);
    let final_q_end = (q_pos + seed_len - 1) + ext.r_q;
    let final_t_start = t_start_idx.saturating_sub(ext.r_t);
    let final_t_end = (t_start_idx + seed_len - 1) + ext.l_t;

    let original_len = ctx.index.get_sequence_len(candidate.target_idx);
    let (out_t_start, out_t_end, strand_char) = match candidate.strand {
        Strand::Reverse => {
            let fwd_start = original_len - 1 - final_t_end;
            let fwd_end = original_len - 1 - final_t_start;
            (fwd_start + 1, fwd_end + 1, '-')
        }
        Strand::Forward => (final_t_start + 1, final_t_end + 1, '+'),
    };

    SearchHit {
        query_id: query.id.into(),
        target_id: ctx.index.get_id(candidate.target_idx).into(),
        q_start: final_q_start,
        q_end: final_q_end,
        t_start: final_t_start,
        t_end: final_t_end,
        output_q_start: final_q_start + 1,
        output_q_end: final_q_end + 1,
        output_t_start: out_t_start,
        output_t_end: out_t_end,
        strand: strand_char.into(),
        energy: ext.score.into(),
        alignment: ext.alignment,
        flank_5: String::new(),
        flank_3: String::new(),
    }
}

#[derive(Debug, Clone)]
pub(super) struct ExtensionResult {
    pub score: f64,
    pub alignment: Alignment,

    pub l_q: usize,
    pub l_t: usize,
    pub r_q: usize,
    pub r_t: usize,
}

pub(super) fn extend_seed(
    ctx: &mut SearchContext<'_, '_>,
    q_seq: &Sequence,
    t_seq: &Sequence,
    candidate: &SeedCandidate,
) -> Option<ExtensionResult> {
    let q_pos = candidate.query_pos;
    let t_pos = candidate.target_start;
    let len = candidate.len;
    let opts = &ctx.args.extend;

    // Get interval bounds from SeedSpec for maximality constraint
    // For -s 8: interval is [0, q_len), maximality check uses full query
    // For -s 1:8: interval is [0, 8), seeds at edges are maximal within interval
    let (interval_start, interval_end) = ctx
        .args
        .seed
        .seed
        .normalize(q_seq.len())
        .map(|(s1, e1, _)| (s1 - 1, e1)) // Convert to 0-based start, exclusive end
        .unwrap_or((0, q_seq.len())); // Fallback: full query

    // Wrap sequences for clean base access (used throughout function)
    let query = q_seq;
    let target = t_seq;

    // Select pair matrix based on wobble mode:
    // - AllowWobble: use PAIR_MAT (G-U wobble pairs are valid)
    // - Strict: use PAIR_MAT_NO_GU (G-U wobble pairs are NOT valid)
    let pair_mat = match ctx.args.seed.pairing {
        SeedPairingMode::AllowWobble => &PAIR_MAT,
        SeedPairingMode::Strict => &PAIR_MAT_NO_GU,
    };

    // MAXIMALITY CHECK
    // Skip non-maximal seeds: if the seed can be extended by a valid base pair
    // on either end, it's a sub-seed of a longer match and will

    trace!(
        "{} ENTERING extend_seed q_pos={} t_pos={} len={} delta_g={}",
        SearchStage::Extend,
        q_pos,
        t_pos,
        len,
        opts.delta_g
    );

    // 1. Left extendable? Check if q[q_pos-1] pairs with t[t_pos+len]
    // Use interval_start to constrain: seed touching interval left edge is maximal on left
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

    // 2. Right extendable? Check if q[q_pos+len] pairs with t[t_pos-1]
    // Use interval_end to constrain: seed touching interval right edge is maximal on right
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

    // Check extendability: seed is extendable if there's room on at least one side
    let left_extendable = q_pos > 0 && t_pos + len < t_seq.len();
    let right_extendable = q_pos + len < q_seq.len() && t_pos > 0;
    let seed_extendable = left_extendable || right_extendable;

    // Only enter extension if max_extension is nonzero AND seed is extendable
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

    // Seed energy calculation using centralized EnergyModel
    let seed_energy_raw = ctx.energy.seed_energy(query, target, q_pos, t_match_end, len);

    // Only build interaction string when trace logging is enabled (avoids allocation in hot path)
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

    // If not entering extension, return seed-only result
    // Need to add terminal penalties at both ends of the seed
    if !do_extension {
        // 5' terminal: Gap->q[0] paired with Gap->t[t_end]
        let term_5p = ctx.energy.terminal_5p(query[q_pos], target[t_match_end]);
        // 3' terminal: q[q_end]->Gap paired with t[t_pos]->Gap
        let term_3p = ctx.energy.terminal_3p(query[q_pos + len - 1], target[t_pos]);

        let total_raw = seed_energy_raw + term_5p + term_3p;
        let final_score = ctx.energy.to_kcal(total_raw);
        let seed_alignment = build_seed_alignment(query, target, q_pos, t_match_end, len);

        return Some(ExtensionResult {
            score: final_score,
            alignment: Alignment::new(&[], &seed_alignment, &[]),
            l_q: 0,
            l_t: 0,
            r_q: 0,
            r_t: 0,
        });
    }

    // DP Left: Extend Query Left (5'), Target Right (3')
    let left_res = ctx
        .extender
        .extend_left(query, target, q_pos, t_match_end, safe_ext);

    // DP Right: Extend Query Right (3'), Target Left (5')
    let right_res = ctx
        .extender
        .extend_right(query, target, q_pos + len - 1, t_pos, safe_ext);

    // Use OLD embedded DP scores for C parity (new dp module has differences)
    let final_score = ctx
        .energy
        .to_kcal(seed_energy_raw + left_res.score + right_res.score);

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

    // Use SmallVec to avoid heap allocation for typical extension sizes
    use smallvec::SmallVec;
    let mut left_alignment: SmallVec<[Pairing; 64]> = SmallVec::new();

    // Left Trace: replay from extension end back toward seed
    let mut li = left_res.q_len;
    let mut lj = left_res.t_len;
    for step in &left_res.trace {
        let pairing = match step {
            dp::DpOp::Match | dp::DpOp::Stop => {
                let p = Pairing::from_bases(
                    left_base(query, candidate.query_pos, li),
                    base_or_gap(target, t_match_end + lj),
                );
                li = li.saturating_sub(1);
                lj = lj.saturating_sub(1);
                p
            }
            dp::DpOp::GapQ => {
                let p = Pairing::GapTarget(left_base(query, candidate.query_pos, li));
                li = li.saturating_sub(1);
                p
            }
            dp::DpOp::GapT => {
                let p = Pairing::GapQuery(base_or_gap(target, t_match_end + lj));
                lj = lj.saturating_sub(1);
                p
            }
        };
        left_alignment.push(pairing);
    }

    // 2. Seed itself
    let seed_alignment = build_seed_alignment(query, target, q_pos, t_match_end, len);

    // Right Trace: replay reversed trace (from seed toward extension end)
    let mut right_alignment: SmallVec<[Pairing; 64]> = SmallVec::new();
    let (mut ri, mut rj) = (0, 0);
    let q_base = candidate.query_pos + len - 1;
    for step in right_res.trace.iter().rev() {
        let pairing = match step {
            dp::DpOp::Match | dp::DpOp::Stop => {
                ri += 1;
                rj += 1;
                Pairing::from_bases(base_or_gap(query, q_base + ri), left_base(target, t_pos, rj))
            }
            dp::DpOp::GapQ => {
                ri += 1;
                Pairing::GapTarget(base_or_gap(query, q_base + ri))
            }
            dp::DpOp::GapT => {
                rj += 1;
                Pairing::GapQuery(left_base(target, t_pos, rj))
            }
        };
        right_alignment.push(pairing);
    }

    let alignment = Alignment::new(&left_alignment, &seed_alignment, &right_alignment);

    Some(ExtensionResult {
        score: final_score,
        alignment,

        l_q: left_res.q_len,
        l_t: left_res.t_len,
        r_q: right_res.q_len,
        r_t: right_res.t_len,
    })
}
