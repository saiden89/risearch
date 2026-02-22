//! Search module - finds miRNA-target interactions.
//!
//! Pipeline:
//! 1) Entry over targets and queries
//! 2) Seed enumeration for one query/target pair
//! 3) Optional extension of each seed
//! 4) Materialize and emit final hits

use anyhow::{Context, Result};
use log::{info, trace};
use rayon::prelude::*;
use smallvec::SmallVec;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;

use crate::alignment::{Alignment, PairClass};
use crate::config::{
    ExtendConfig, FilterConfig, Matrix, OutputFormat, ScoreConfig, SearchArgs, SearchAxis,
};
use crate::dp::{DpConfig, DpExtender, DpView};
use crate::dsm::{pair_mat, seed_energy, terminal_3p, terminal_5p, DsmModel, T04, T99};
use crate::index::store::{TargetStore, TargetView};
use crate::registry::{QueryData, QueryRegistry};
use crate::seed::{for_each_seed_one_target, SeedHit, TargetSeedView};
use crate::seq::Sequence;
use crate::types::{Base, Energy, Interval, Strand};

// =============================================================================
// PUBLIC API
// =============================================================================

/// A search hit representing a miRNA-target interaction.
#[derive(Debug, Clone)]
pub struct SearchHit {
    pub query_idx: u32,
    pub target_idx: u32,
    pub q_start: usize,
    pub q_end: usize,
    pub t_start: usize,
    pub t_end: usize,
    pub strand: Strand,
    pub energy: Energy,
    pub seed_start: Option<usize>,
    pub seed_end: Option<usize>,
    pub alignment: Option<Alignment>,
    pub flank_5: Sequence,
    pub flank_3: Sequence,
}

/// Run search against mmap-backed target store and emit hits through callback.
pub fn run_search<F>(
    queries: &QueryRegistry,
    store: &TargetStore,
    opts: &SearchArgs,
    mut on_hit: F,
) -> Result<usize>
where
    F: FnMut(SearchHit) -> Result<()>,
{
    if store.is_empty() || queries.is_empty() {
        return Ok(0);
    }

    let axis = match opts.axis {
        SearchAxis::Auto => {
            if store.len() == 1 {
                Axis::ByQuery
            } else {
                let workers = rayon::current_num_threads().max(1);
                let target_threshold = workers.div_ceil(2);
                if store.len() >= target_threshold || queries.len() == 1 {
                    Axis::ByTarget
                } else {
                    Axis::ByPair
                }
            }
        }
        SearchAxis::Target => Axis::ByTarget,
        SearchAxis::Query => Axis::ByQuery,
        SearchAxis::Pair => Axis::ByPair,
    };

    info!(
        "Starting search: {} queries x {} targets, axis={:?}, seed={:?}, max_ext={}, delta_g={}",
        queries.len(),
        store.len(),
        axis,
        opts.seed.seed,
        opts.extend.max_extension,
        opts.filter.delta_g
    );

    let (tx, rx) = mpsc::sync_channel::<SearchHit>(HIT_CHANNEL_CAPACITY);
    let cancelled = AtomicBool::new(false);

    let total = std::thread::scope(|scope| -> Result<usize> {
        let cancel = &cancelled;
        let producer = scope.spawn(move || -> Result<()> {
            match (opts.score.matrix, axis) {
                (Matrix::T04, Axis::ByTarget) => {
                    produce_by_target::<T04>(queries, store, opts, tx, cancel)
                }
                (Matrix::T04, Axis::ByQuery) => {
                    produce_by_query::<T04>(queries, store, opts, tx, cancel)
                }
                (Matrix::T04, Axis::ByPair) => {
                    produce_by_pair::<T04>(queries, store, opts, tx, cancel)
                }
                (Matrix::T99, Axis::ByTarget) => {
                    produce_by_target::<T99>(queries, store, opts, tx, cancel)
                }
                (Matrix::T99, Axis::ByQuery) => {
                    produce_by_query::<T99>(queries, store, opts, tx, cancel)
                }
                (Matrix::T99, Axis::ByPair) => {
                    produce_by_pair::<T99>(queries, store, opts, tx, cancel)
                }
            }
        });

        let mut emitted = 0usize;
        let mut callback_err: Option<anyhow::Error> = None;

        while let Ok(hit) = rx.recv() {
            if callback_err.is_some() {
                continue;
            }
            match on_hit(hit) {
                Ok(()) => emitted += 1,
                Err(err) => {
                    callback_err = Some(err);
                    cancelled.store(true, Ordering::Relaxed);
                }
            }
        }

        match producer.join() {
            Ok(result) => result?,
            Err(_) => return Err(anyhow::anyhow!("Parallel search producer thread panicked")),
        }

        if let Some(err) = callback_err {
            return Err(err);
        }

        Ok(emitted)
    })?;

    info!("Search complete: {} hits", total);
    Ok(total)
}

// =============================================================================
// ORCHESTRATION
// =============================================================================

const HIT_CHANNEL_CAPACITY: usize = 2048;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Axis {
    ByTarget,
    ByQuery,
    ByPair,
}

/// Reusable per-run state.
struct SearchState<M: DsmModel> {
    extender: DpExtender<M>,
    dp_cfg: DpConfig,
}

impl<M: DsmModel> SearchState<M> {
    /// Build per-worker DP state from the run-level score/extension configuration.
    fn new(score_cfg: &ScoreConfig, extend_cfg: &ExtendConfig) -> Self {
        let dp_cfg = DpConfig::from((score_cfg, extend_cfg));
        Self {
            extender: DpExtender::<M>::from_config(dp_cfg),
            dp_cfg,
        }
    }
}

/// Parallel target-major producer:
/// each worker processes one target across all queries.
fn produce_by_target<M: DsmModel>(
    queries: &QueryRegistry,
    store: &TargetStore,
    opts: &SearchArgs,
    tx: mpsc::SyncSender<SearchHit>,
    cancelled: &AtomicBool,
) -> Result<()> {
    let target_ids: Vec<u32> = store.iter_meta().map(|(idx, _)| idx).collect();
    target_ids.into_par_iter().try_for_each_init(
        || (SearchState::<M>::new(&opts.score, &opts.extend), tx.clone()),
        |(state, tx), target_idx| -> Result<()> {
            if cancelled.load(Ordering::Relaxed) {
                return Ok(());
            }

            let target = store
                .target_view(target_idx as usize)
                .with_context(|| format!("Failed to load target #{}", target_idx))?;

            for (query_idx, query) in queries.entries().iter().enumerate() {
                if cancelled.load(Ordering::Relaxed) {
                    break;
                }

                let mut emit_hit = |hit: SearchHit| -> Result<()> {
                    if cancelled.load(Ordering::Relaxed) {
                        return Ok(());
                    }
                    tx.send(hit)
                        .map_err(|_| anyhow::anyhow!("Search hit channel disconnected"))?;
                    Ok(())
                };

                search_query_against_target::<M, _>(
                    query_idx as u32,
                    query,
                    &target,
                    target_idx,
                    opts,
                    state,
                    &mut emit_hit,
                )?;
            }

            Ok(())
        },
    )
}

/// Parallel query-major producer for single-target workloads.
///
/// This path avoids underutilization when only one target is present.
fn produce_by_query<M: DsmModel>(
    queries: &QueryRegistry,
    store: &TargetStore,
    opts: &SearchArgs,
    tx: mpsc::SyncSender<SearchHit>,
    cancelled: &AtomicBool,
) -> Result<()> {
    let (target_idx, _) = store
        .iter_meta()
        .next()
        .ok_or_else(|| anyhow::anyhow!("Target index is empty"))?;
    let target = store
        .target_view(target_idx as usize)
        .with_context(|| format!("Failed to load target #{}", target_idx))?;

    (0..queries.len()).into_par_iter().try_for_each_init(
        || (SearchState::<M>::new(&opts.score, &opts.extend), tx.clone()),
        |(state, tx), query_idx| {
            if cancelled.load(Ordering::Relaxed) {
                return Ok(());
            }

            let mut emit_hit = |hit: SearchHit| -> Result<()> {
                if cancelled.load(Ordering::Relaxed) {
                    return Ok(());
                }
                tx.send(hit)
                    .map_err(|_| anyhow::anyhow!("Search hit channel disconnected"))?;
                Ok(())
            };

            search_query_against_target::<M, _>(
                query_idx as u32,
                &queries.entries()[query_idx],
                &target,
                target_idx,
                opts,
                state,
                &mut emit_hit,
            )?;

            Ok(())
        },
    )
}

/// Parallel pair-major producer:
/// each worker processes one (target, query chunk) task.
fn produce_by_pair<M: DsmModel>(
    queries: &QueryRegistry,
    store: &TargetStore,
    opts: &SearchArgs,
    tx: mpsc::SyncSender<SearchHit>,
    cancelled: &AtomicBool,
) -> Result<()> {
    let query_entries = queries.entries();
    let query_count = query_entries.len();
    let target_ids: Vec<u32> = store.iter_meta().map(|(idx, _)| idx).collect();
    if query_count == 0 || target_ids.is_empty() {
        return Ok(());
    }

    let workers = rayon::current_num_threads().max(1);
    let chunks_per_target = (workers * 8).div_ceil(target_ids.len()).max(1);
    let q_chunk = query_count.div_ceil(chunks_per_target).max(1);
    let chunks_per_target_actual = query_count.div_ceil(q_chunk);

    let mut tasks = Vec::with_capacity(target_ids.len() * chunks_per_target_actual);
    for target_idx in target_ids {
        let mut q_start = 0usize;
        while q_start < query_count {
            let q_end = (q_start + q_chunk).min(query_count);
            tasks.push((target_idx, q_start, q_end));
            q_start = q_end;
        }
    }

    tasks.into_par_iter().try_for_each_init(
        || (SearchState::<M>::new(&opts.score, &opts.extend), tx.clone()),
        |(state, tx), (target_idx, q_start, q_end)| -> Result<()> {
            if cancelled.load(Ordering::Relaxed) {
                return Ok(());
            }

            let target = store
                .target_view(target_idx as usize)
                .with_context(|| format!("Failed to load target #{}", target_idx))?;

            for query_idx in q_start..q_end {
                if cancelled.load(Ordering::Relaxed) {
                    break;
                }

                let mut emit_hit = |hit: SearchHit| -> Result<()> {
                    if cancelled.load(Ordering::Relaxed) {
                        return Ok(());
                    }
                    tx.send(hit)
                        .map_err(|_| anyhow::anyhow!("Search hit channel disconnected"))?;
                    Ok(())
                };

                search_query_against_target::<M, _>(
                    query_idx as u32,
                    &query_entries[query_idx],
                    &target,
                    target_idx,
                    opts,
                    state,
                    &mut emit_hit,
                )?;
            }

            Ok(())
        },
    )
}

/// Query/target-scoped context that does not change across seed candidates.
struct QueryTargetCtx<'a> {
    query_idx: u32,
    query_bases: &'a [Base],
    seed_interval: Interval,
    include_alignment: bool,
    pair_matrix: &'static [[u8; 6]; 6],
    filter_cfg: &'a FilterConfig,
    target_len: usize,
}

/// Enumerate seeds for a single query-target pair, optionally extend them, and emit hits.
fn search_query_against_target<M, F>(
    query_idx: u32,
    query: &QueryData,
    target: &TargetView<'_>,
    target_idx: u32,
    opts: &SearchArgs,
    state: &mut SearchState<M>,
    on_hit: &mut F,
) -> Result<usize>
where
    M: DsmModel,
    F: FnMut(SearchHit) -> Result<()>,
{
    let target_seed_view = TargetSeedView {
        combined_seq: target.combined_seq,
        combined_sa: target.combined_sa,
        seq_len: target.seq_len,
    };

    let ctx = QueryTargetCtx {
        query_idx,
        query_bases: query.sequence(),
        seed_interval: query.seed_interval(),
        include_alignment: opts.output.format != OutputFormat::Minimal,
        pair_matrix: pair_mat(opts.seed.allows_wobble()),
        filter_cfg: &opts.filter,
        target_len: target.seq_len,
    };

    let t_forward = &target.combined_seq[..target.seq_len];
    let t_reverse = &target.combined_seq[target.seq_len + 1..2 * target.seq_len + 1];

    let mut emitted = 0usize;
    let mut callback_err: Option<anyhow::Error> = None;

    for_each_seed_one_target(query, target_idx, &target_seed_view, &opts.seed, |seed| {
        if callback_err.is_some() {
            return;
        }

        let target_bases = match seed.strand {
            Strand::Forward => t_forward,
            Strand::Reverse => t_reverse,
        };

        let Some(hit) =
            build_hit_from_seed::<M>(&mut state.extender, state.dp_cfg, &ctx, &seed, target_bases)
        else {
            return;
        };

        match on_hit(hit) {
            Ok(()) => emitted += 1,
            Err(err) => callback_err = Some(err),
        }
    });

    if let Some(err) = callback_err {
        return Err(err);
    }

    Ok(emitted)
}

/// Build a finalized `SearchHit` from a seed if extension and energy filters pass.
fn build_hit_from_seed<M: DsmModel>(
    extender: &mut DpExtender<M>,
    dp_cfg: DpConfig,
    ctx: &QueryTargetCtx<'_>,
    seed: &SeedHit,
    target_bases: &[Base],
) -> Option<SearchHit> {
    if seed.target_start + seed.seed_len.get() > target_bases.len() {
        return None;
    }

    let extension = compute_seed_extension::<M>(
        extender,
        dp_cfg,
        ctx.query_bases,
        target_bases,
        seed,
        ctx.seed_interval,
        ctx.filter_cfg,
        ctx.pair_matrix,
        ctx.include_alignment,
    )?;

    if extension.score > ctx.filter_cfg.delta_g {
        return None;
    }

    Some(SearchHit::new(
        ctx.query_idx,
        ctx.query_bases,
        target_bases,
        seed,
        &extension,
        ctx.include_alignment,
        ctx.target_len,
    ))
}

// =============================================================================
// EXTENSION
// =============================================================================

struct SeedExtension {
    score: f64,
    l_q: usize,
    l_t: usize,
    r_q: usize,
    r_t: usize,
    left_pairs: SmallVec<[PairClass; 64]>,
    right_pairs: SmallVec<[PairClass; 64]>,
}

/// Compute optional left/right DP extension around a seed and return extension metadata.
fn compute_seed_extension<M: DsmModel>(
    extender: &mut DpExtender<M>,
    dp_cfg: DpConfig,
    query_bases: &[Base],
    target_bases: &[Base],
    seed: &SeedHit,
    seed_interval: Interval,
    filter_cfg: &FilterConfig,
    pair_matrix: &'static [[u8; 6]; 6],
    with_traceback: bool,
) -> Option<SeedExtension> {
    let penalty = dp_cfg.penalty_raw();
    let q_pos = seed.query_pos;
    let t_pos = seed.target_start;
    let len = seed.seed_len.get();

    // Maximality check: left side.
    if q_pos > seed_interval.start && t_pos + len < target_bases.len() {
        let p_class = pair_matrix[query_bases[q_pos - 1].idx()][target_bases[t_pos + len].idx()];
        if p_class != 0 && !filter_cfg.no_max_prune {
            trace!("Filtered: non-maximal left");
            return None;
        }
    }

    // Maximality check: right side.
    if q_pos + len < seed_interval.end && t_pos > 0 {
        let p_class = pair_matrix[query_bases[q_pos + len].idx()][target_bases[t_pos - 1].idx()];
        if p_class != 0 && !filter_cfg.no_max_prune {
            trace!("Filtered: non-maximal right");
            return None;
        }
    }

    let t_match_end = t_pos + len - 1;
    let max_ext = dp_cfg.max_extension();
    let seed_e = seed_energy::<M>(query_bases, target_bases, q_pos, t_match_end, len, penalty);

    let can_extend_left = q_pos > 0 && t_pos + len < target_bases.len();
    let can_extend_right = q_pos + len < query_bases.len() && t_pos > 0;

    if max_ext == 0 || (!can_extend_left && !can_extend_right) {
        let term_5p = terminal_5p::<M>(query_bases[q_pos], target_bases[t_match_end], penalty);
        let term_3p = terminal_3p::<M>(query_bases[q_pos + len - 1], target_bases[t_pos], penalty);
        let nt_count = (2 * len) as i32;
        return Some(SeedExtension {
            score: M::to_kcal(seed_e + term_5p + term_3p + nt_count * penalty),
            l_q: 0,
            l_t: 0,
            r_q: 0,
            r_t: 0,
            left_pairs: SmallVec::new(),
            right_pairs: SmallVec::new(),
        });
    }

    let (l_score, l_q, l_t, left_pairs) = {
        let view = DpView::<M>::left(query_bases, target_bases, q_pos, t_match_end, max_ext);
        let result = extender.extend(&view);
        let pairs = if with_traceback {
            result.traceback(&view)
        } else {
            SmallVec::new()
        };
        (result.score, result.q_len, result.t_len, pairs)
    };

    let (r_score, r_q, r_t, right_pairs) = {
        let view = DpView::<M>::right(query_bases, target_bases, q_pos + len - 1, t_pos, max_ext);
        let result = extender.extend(&view);
        let pairs = if with_traceback {
            result.traceback(&view)
        } else {
            SmallVec::new()
        };
        (result.score, result.q_len, result.t_len, pairs)
    };

    let nt_count = (l_q + l_t + r_q + r_t + 2 * len) as i32;
    Some(SeedExtension {
        score: M::to_kcal(seed_e + l_score + r_score + nt_count * penalty),
        l_q,
        l_t,
        r_q,
        r_t,
        left_pairs,
        right_pairs,
    })
}

// =============================================================================
// HIT MATERIALIZATION
// =============================================================================

impl SearchHit {
    /// Materialize a public `SearchHit` from canonical internal inputs.
    fn new(
        query_idx: u32,
        query_bases: &[Base],
        target_bases: &[Base],
        seed: &SeedHit,
        ext: &SeedExtension,
        include_alignment: bool,
        original_target_len: usize,
    ) -> Self {
        let q_pos = seed.query_pos;
        let t_start = seed.target_start;
        let len = seed.seed_len.get();

        let final_q_start = q_pos.saturating_sub(ext.l_q);
        let final_q_end = (q_pos + len - 1) + ext.r_q;
        let final_t_start = t_start.saturating_sub(ext.r_t);
        let final_t_end = (t_start + len - 1) + ext.l_t;

        let (final_t_start, final_t_end, strand) = match seed.strand {
            Strand::Forward => (final_t_start, final_t_end, Strand::Forward),
            Strand::Reverse => {
                let fwd_start = original_target_len - 1 - final_t_end;
                let fwd_end = original_target_len - 1 - final_t_start;
                (fwd_start, fwd_end, Strand::Reverse)
            }
        };

        let alignment = if include_alignment {
            let t_match_end = seed.target_start + len - 1;
            let mut seed_pairs: SmallVec<[PairClass; 64]> = SmallVec::with_capacity(len);
            for i in 0..len {
                seed_pairs.push(PairClass::from_bases(
                    query_bases[q_pos + i],
                    target_bases[t_match_end - i],
                ));
            }
            Some(Alignment::new(
                &ext.left_pairs,
                &seed_pairs,
                &ext.right_pairs,
            ))
        } else {
            None
        };

        let (seed_start, seed_end) = if include_alignment {
            let start = ext.left_pairs.len();
            (Some(start), Some(start + len))
        } else {
            (None, None)
        };

        Self {
            query_idx,
            target_idx: seed.target_id.0,
            q_start: final_q_start,
            q_end: final_q_end,
            t_start: final_t_start,
            t_end: final_t_end,
            strand,
            energy: ext.score.into(),
            seed_start,
            seed_end,
            alignment,
            flank_5: Sequence::from(Vec::new()),
            flank_3: Sequence::from(Vec::new()),
        }
    }
}
