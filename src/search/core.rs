use anyhow::{Result, bail};
use log::{info, trace, warn};
use rayon::prelude::*;

use crate::config::SearchArgs;
use crate::dp;
use crate::dsm::{PAIR_MAT, PAIR_MAT_NO_GU};
use crate::seed::{SeedCandidate, build_seed_alignment};
use crate::seq::Seq;
use crate::types::{Alignment, Pairing, SeedPairingMode, Strand};

use super::output::bytes_to_rna_string;
use super::{
    FilterReason, MAX_DP_EXT, SaIndex, SearchContext, SearchHit, SearchStage,
    THREAD_EXTENDER,
};

/// Core search logic - finds and deduplicates hits for all queries.
/// Uses parallel iteration over queries for multi-core utilization.
fn search_core(
    queries: &[(String, Vec<u8>)],
    index: &SaIndex<'_>,
    opts: &SearchArgs,
) -> Result<Vec<SearchHit>> {
    for (q_id, q_seq) in queries {
        if let Err(err) = opts.seed.seed.normalize(q_seq.len()) {
            bail!("Invalid seed spec for query '{}': {}", q_id, err);
        }
    }
    info!(
        "Starting search: {} queries, seed={:?}, max_ext={}, delta_g={}",
        queries.len(),
        opts.seed.seed,
        opts.extend.max_extension,
        opts.extend.delta_g
    );

    // Parallel query processing (per-query is independent).
    // Configure worker count via the global rayon thread pool.
    let per_query_hits: Vec<Vec<SearchHit>> = queries
        .par_iter()
        .map(|(q_id, q_seq)| {
            // Borrow thread-local extender for this query
            THREAD_EXTENDER.with(|ext| {
                let mut extender = ext.borrow_mut();
                let mut ctx = SearchContext::with_extender(index, opts, &mut *extender);

                trace!("{} id={} len={}", SearchStage::Input, q_id, q_seq.len());
                trace!("{} {}", SearchStage::Input, String::from_utf8_lossy(q_seq));

                let seeds = match find_seeds_for_query(q_seq, &mut ctx) {
                    Ok(s) => s,
                    Err(_) => return Vec::new(),
                };
                trace!("{} {} candidates found", SearchStage::Input, seeds.len());

                let mut hits = Vec::new();
                for candidate in &seeds {
                    trace!(
                        "{} q_pos={} t_idx={} t_start={} len={} strand={:?}",
                        SearchStage::Seed,
                        candidate.query_pos,
                        candidate.target_idx,
                        candidate.target_start,
                        candidate.len,
                        candidate.strand
                    );
                    if let Some(hit) = process_candidate(
                        crate::types::Query::new(q_id, q_seq),
                        candidate,
                        &mut ctx,
                    ) {
                        trace!(
                            "{} q={}-{} t={}-{} E={:.2}",
                            SearchStage::Output,
                            hit.q_start,
                            hit.q_end,
                            hit.t_start,
                            hit.t_end,
                            hit.energy
                        );
                        hits.push(hit);
                    }
                }
                hits
            })
        })
        .collect();

    let all_hits: Vec<SearchHit> = per_query_hits.into_iter().flatten().collect();

    info!("Search complete: {} hits", all_hits.len());

    Ok(all_hits)
}

/// Run search and return hits.
pub fn run_search(
    queries: &[(String, Vec<u8>)],
    index: &SaIndex<'_>,
    opts: &SearchArgs,
) -> Result<Vec<SearchHit>> {
    search_core(queries, index, opts)
}

fn find_seeds_for_query(
    q_seq: &[u8],
    ctx: &mut SearchContext<'_, '_>,
) -> Result<Vec<SeedCandidate>> {
    use crate::seed;

    trace!(
        "{} q_len={} spec={:?} pairing={:?}",
        SearchStage::Seed,
        q_seq.len(),
        ctx.args.seed.seed,
        ctx.args.seed.pairing
    );

    Ok(seed::find_seeds(q_seq, ctx.index.index, &ctx.args.seed))
}

fn process_candidate(
    query: crate::types::Query<'_>,
    candidate: &SeedCandidate,
    ctx: &mut SearchContext<'_, '_>,
) -> Option<SearchHit> {
    let t_idx = candidate.target_idx;
    let t_seq: &[u8] = match candidate.strand {
        Strand::Reverse => ctx.index.get_sequence_rc(t_idx),
        Strand::Forward => ctx.index.get_sequence(t_idx),
    };

    let t_start_idx = candidate.target_start;
    let seed_len = candidate.len;
    let q_pos = candidate.query_pos;

    trace!(
        "{} q_id={} t_idx={} q_pos={} t_start={} seed_len={} strand={:?}",
        SearchStage::Extend,
        query.id,
        t_idx,
        q_pos,
        t_start_idx,
        seed_len,
        candidate.strand
    );

    if t_start_idx + seed_len > t_seq.len() {
        ctx.stats.record_filter(FilterReason::SeedOutOfBounds);
        warn!(
            "{} FILTERED reason={:?} q_pos={} t_start={} seed_len={} t_len={}",
            SearchStage::Extend,
            FilterReason::SeedOutOfBounds,
            q_pos,
            t_start_idx,
            seed_len,
            t_seq.len()
        );
        return None;
    }

    // Call extend_seed
    let Some(ext) = extend_seed(ctx, query.seq, t_seq, candidate) else {
        return None;
    };

    let score = ext.score;
    if score > ctx.args.extend.delta_g {
        ctx.stats.record_filter(FilterReason::EnergyAboveThreshold);
        warn!(
            "{} FILTERED reason={:?} q_pos={} t_start={} seed_len={} score={:.2}",
            SearchStage::Extend,
            FilterReason::EnergyAboveThreshold,
            q_pos,
            t_start_idx,
            seed_len,
            score
        );
        return None;
    }

    // Final coordinates (0-based)
    let final_q_start = q_pos - ext.l_q;
    let final_q_end = (q_pos + seed_len - 1) + ext.r_q;
    let final_t_start = t_start_idx - ext.r_t;
    let final_t_end = (t_start_idx + seed_len - 1) + ext.l_t;

    // Output coordinates (1-based, strand-aware for target)
    let original_len = ctx.index.get_sequence_len(candidate.target_idx);
    let (out_t_start, out_t_end, strand_char) = match candidate.strand {
        Strand::Reverse => {
            let fwd_start = original_len - 1 - final_t_end;
            let fwd_end = original_len - 1 - final_t_start;
            (fwd_start + 1, fwd_end + 1, '-')
        }
        Strand::Forward => (final_t_start + 1, final_t_end + 1, '+'),
    };

    // Create flanking sequences
    let ctx_len = 20;
    let flank_5 = bytes_to_rna_string(
        &t_seq[final_t_start.saturating_sub(ctx_len)..final_t_start],
        true,
    );
    let flank_3 = if final_t_end + 1 < t_seq.len() {
        let t_3_end = (final_t_end + 1 + ctx_len).min(t_seq.len());
        bytes_to_rna_string(&t_seq[final_t_end + 1..t_3_end], false)
    } else {
        String::new()
    };

    Some(SearchHit {
        query_id: query.id.into(),
        target_id: ctx.index.get_id(t_idx).into(),

        q_start: final_q_start,
        q_end: final_q_end,
        t_start: final_t_start,
        t_end: final_t_end,
        output_q_start: final_q_start + 1, // 1-based for output
        output_q_end: final_q_end + 1,     // 1-based for output
        output_t_start: out_t_start,
        output_t_end: out_t_end,
        strand: strand_char.into(),
        energy: score.into(),
        alignment: ext.alignment,
        flank_5,
        flank_3,
    })
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
    q_seq: &[u8],
    t_seq: &[u8],
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
    let query = Seq::forward(q_seq);
    let target = Seq::new(t_seq, candidate.strand);

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
        let q_prev = query.base(q_pos - 1).idx();
        let t_next = target.base(t_pos + len).idx();
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
        let q_next = query.base(q_pos + len).idx();
        let t_prev = target.base(t_pos - 1).idx();
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
    let seed_energy_raw = ctx
        .energy
        .seed_energy(&query, &target, q_pos, t_match_end, len);

    // Only build interaction string when trace logging is enabled (avoids allocation in hot path)
    if log::log_enabled!(log::Level::Trace) {
        let mut seed_int_str = String::with_capacity(len);
        for k in 0..len {
            let qc = query.base(q_pos + k);
            let tc = target.base(t_match_end - k);
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
        let term_5p = ctx
            .energy
            .terminal_5p(query.base(q_pos), target.base(t_match_end));
        // 3' terminal: q[q_end]->Gap paired with t[t_pos]->Gap
        let term_3p = ctx
            .energy
            .terminal_3p(query.base(q_pos + len - 1), target.base(t_pos));

        let total_raw = seed_energy_raw + term_5p + term_3p;
        let final_score = ctx.energy.to_kcal(total_raw);
        let seed_alignment = build_seed_alignment(&query, &target, q_pos, t_match_end, len);

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
        .extend_left(&query, &target, q_pos, t_match_end, safe_ext);

    // DP Right: Extend Query Right (3'), Target Left (5')
    let right_res = ctx
        .extender
        .extend_right(&query, &target, q_pos + len - 1, t_pos, safe_ext);

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
                    query.left(candidate.query_pos, li),
                    target.base_or_gap(t_match_end + lj),
                );
                li = li.saturating_sub(1);
                lj = lj.saturating_sub(1);
                p
            }
            dp::DpOp::GapQ => {
                let p = Pairing::GapTarget(query.left(candidate.query_pos, li));
                li = li.saturating_sub(1);
                p
            }
            dp::DpOp::GapT => {
                let p = Pairing::GapQuery(target.base_or_gap(t_match_end + lj));
                lj = lj.saturating_sub(1);
                p
            }
        };
        left_alignment.push(pairing);
    }

    // 2. Seed itself
    let seed_alignment = build_seed_alignment(&query, &target, q_pos, t_match_end, len);

    // Right Trace: replay reversed trace (from seed toward extension end)
    let mut right_alignment: SmallVec<[Pairing; 64]> = SmallVec::new();
    let (mut ri, mut rj) = (0, 0);
    let q_base = candidate.query_pos + len - 1;
    for step in right_res.trace.iter().rev() {
        let pairing = match step {
            dp::DpOp::Match | dp::DpOp::Stop => {
                ri += 1;
                rj += 1;
                Pairing::from_bases(query.base_or_gap(q_base + ri), target.left(t_pos, rj))
            }
            dp::DpOp::GapQ => {
                ri += 1;
                Pairing::GapTarget(query.base_or_gap(q_base + ri))
            }
            dp::DpOp::GapT => {
                rj += 1;
                Pairing::GapQuery(target.left(t_pos, rj))
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
