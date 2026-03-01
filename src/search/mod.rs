//! Search module - finds miRNA-target interactions.
//!
//! Pipeline:
//! 1) Parallel iteration over queries (global SA traversal per query)
//! 2) Seed enumeration across all targets via single SA traversal
//! 3) Optional extension of each seed
//! 4) Materialize and emit final hits

use anyhow::Result;
use log::info;
use rayon::prelude::*;
use smallvec::SmallVec;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;

use crate::alignment::{Alignment, PairClass};
use crate::config::{ExtendConfig, FilterConfig, Matrix, OutputFormat, ScoreConfig, SearchArgs};
use crate::dp::{DpConfig, DpExtender, DpView};
use crate::dsm::{pair_mat, stack_with_penalty, terminal_3p, terminal_5p, DsmModel, T04, T99};
use crate::index::store::{GlobalView, TargetStore};
use crate::registry::{QueryData, QueryRegistry};
use crate::seed::{for_each_seed, SeedHit};
use crate::seq::{SeqView, Sequence};
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
    mut on_chunk: F,
) -> Result<usize>
where
    F: FnMut(OutputChunk) -> Result<()>,
{
    if store.is_empty() || queries.is_empty() {
        return Ok(0);
    }

    let global = store.global_view();

    info!(
        "Starting search: {} queries x {} targets, seed={:?}, max_ext={}, delta_g={}",
        queries.len(),
        store.len(),
        opts.seed.seed,
        opts.extend.max_extension,
        opts.filter.delta_g
    );

    // Fast path for single-worker runs: avoid producer thread + channel fan-in.
    if rayon::current_num_threads() == 1 {
        let total = match opts.score.matrix {
            Matrix::T04 => {
                run_search_single_thread::<T04, _>(queries, store, &global, opts, &mut on_chunk)?
            }
            Matrix::T99 => {
                run_search_single_thread::<T99, _>(queries, store, &global, opts, &mut on_chunk)?
            }
        };
        info!("Search complete: {} hits", total);
        return Ok(total);
    }

    let (tx, rx) = mpsc::sync_channel::<OutputChunk>(CHUNK_CHANNEL_CAPACITY);
    let cancelled = AtomicBool::new(false);

    let total = std::thread::scope(|scope| -> Result<usize> {
        let cancel = &cancelled;
        let producer = scope.spawn(move || -> Result<()> {
            match opts.score.matrix {
                Matrix::T04 => produce_by_query::<T04>(queries, store, &global, opts, tx, cancel),
                Matrix::T99 => produce_by_query::<T99>(queries, store, &global, opts, tx, cancel),
            }
        });

        let mut emitted = 0usize;
        let mut callback_err: Option<anyhow::Error> = None;

        while let Ok(chunk) = rx.recv() {
            if callback_err.is_some() {
                continue;
            }
            let count = chunk.hits;
            match on_chunk(chunk) {
                Ok(()) => emitted += count,
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

const CHUNK_CHANNEL_CAPACITY: usize = 128;
const CHUNK_SIZE_THRESHOLD: usize = 64 * 1024;

/// A batch of pre-formatted hit lines ready for writing.
pub struct OutputChunk {
    pub data: Vec<u8>,
    pub hits: usize,
    /// When multi-file output is active, identifies which query produced this chunk.
    pub query_idx: Option<u32>,
}

/// Reusable per-worker state.
struct SearchState<M: DsmModel> {
    extender: DpExtender<M>,
    dp_cfg: DpConfig,
}

impl<M: DsmModel> SearchState<M> {
    fn new(score_cfg: &ScoreConfig, extend_cfg: &ExtendConfig) -> Self {
        let dp_cfg = DpConfig::from((score_cfg, extend_cfg));
        Self {
            extender: DpExtender::<M>::from_config(dp_cfg),
            dp_cfg,
        }
    }
}

struct FormatState {
    fmt_bufs: crate::output::format::OutputBuffers,
    chunk_data: Vec<u8>,
    chunk_hits: usize,
}

impl FormatState {
    fn new() -> Self {
        Self {
            fmt_bufs: crate::output::format::OutputBuffers::new(),
            chunk_data: Vec::with_capacity(CHUNK_SIZE_THRESHOLD + 1024),
            chunk_hits: 0,
        }
    }

    fn flush<F>(&mut self, query_idx: Option<u32>, on_chunk: &mut F) -> Result<()>
    where
        F: FnMut(OutputChunk) -> Result<()>,
    {
        if self.chunk_hits > 0 {
            let data = std::mem::replace(
                &mut self.chunk_data,
                Vec::with_capacity(CHUNK_SIZE_THRESHOLD + 1024),
            );
            let hits = std::mem::replace(&mut self.chunk_hits, 0);
            on_chunk(OutputChunk {
                data,
                hits,
                query_idx,
            })?;
        }
        Ok(())
    }
}

#[inline]
fn target_transformed_slices<'a>(
    global: &'a GlobalView<'_>,
    target_idx: usize,
) -> (&'a [Base], &'a [Base], usize) {
    let target_len = global.seq_lens[target_idx] as usize;
    let target_offset = global.offsets[target_idx] as usize;
    let t_fwd = &global.combined_seq[target_offset..target_offset + target_len];
    let t_rc =
        &global.combined_seq[target_offset + target_len + 1..target_offset + 2 * target_len + 1];
    (t_fwd, t_rc, target_len)
}

#[inline]
fn append_formatted_hit(
    format: &mut FormatState,
    hit: &SearchHit,
    output_format: OutputFormat,
    store: &TargetStore,
    global: &GlobalView<'_>,
    query_name: &str,
    query_seq: &[Base],
) {
    let target_name = store.get_name(hit.target_idx);
    let target_idx = hit.target_idx as usize;
    let (t_fwd, t_rc, _) = target_transformed_slices(global, target_idx);
    crate::output::format::append_hit_names_vec(
        &mut format.fmt_bufs,
        hit,
        output_format,
        &mut format.chunk_data,
        query_name,
        target_name,
        SeqView::from(query_seq),
        SeqView::from(t_fwd),
        SeqView::from(t_rc),
    );
    format.chunk_hits += 1;
}

/// Process one query end-to-end and stream formatted chunks via `on_chunk`.
fn process_query<M, FC, FO>(
    query_idx: u32,
    query: &QueryData,
    queries: &QueryRegistry,
    store: &TargetStore,
    global: &GlobalView<'_>,
    opts: &SearchArgs,
    state: &mut SearchState<M>,
    format: &mut FormatState,
    flush_after_query: bool,
    is_cancelled: FC,
    on_chunk: &mut FO,
) -> Result<usize>
where
    M: DsmModel,
    FC: Fn() -> bool,
    FO: FnMut(OutputChunk) -> Result<()>,
{
    if is_cancelled() {
        return Ok(0);
    }

    let query_name = queries.get_name(query_idx);
    let query_seq = query.sequence().as_slice();
    let chunk_qi = if opts.output.multifile {
        Some(query_idx)
    } else {
        None
    };
    let mut local_hits = 0usize;

    search_query::<M, _>(
        query_idx,
        query,
        global,
        opts,
        state,
        &mut |hit: SearchHit| -> Result<()> {
            if is_cancelled() {
                return Ok(());
            }

            append_formatted_hit(
                format,
                &hit,
                opts.output.format,
                store,
                global,
                query_name,
                query_seq,
            );
            local_hits += 1;

            if format.chunk_data.len() >= CHUNK_SIZE_THRESHOLD {
                format.flush(chunk_qi, on_chunk)?;
            }
            Ok(())
        },
    )?;

    if flush_after_query {
        format.flush(chunk_qi, on_chunk)?;
    }

    Ok(local_hits)
}

/// Single-thread execution path (no channel fan-in).
fn run_search_single_thread<M, F>(
    queries: &QueryRegistry,
    store: &TargetStore,
    global: &GlobalView<'_>,
    opts: &SearchArgs,
    on_chunk: &mut F,
) -> Result<usize>
where
    M: DsmModel,
    F: FnMut(OutputChunk) -> Result<()>,
{
    let mut state = SearchState::<M>::new(&opts.score, &opts.extend);
    let mut format = FormatState::new();
    let mut total_hits = 0usize;
    let multifile = opts.output.multifile;

    for (query_idx, query) in queries.entries().iter().enumerate() {
        total_hits += process_query::<M, _, _>(
            query_idx as u32,
            query,
            queries,
            store,
            global,
            opts,
            &mut state,
            &mut format,
            multifile,
            || false,
            on_chunk,
        )?;
    }

    format.flush(None, on_chunk)?;
    Ok(total_hits)
}

/// Parallel query-major producer: each worker processes one query against the global SA.
fn produce_by_query<M: DsmModel>(
    queries: &QueryRegistry,
    store: &TargetStore,
    global: &GlobalView<'_>,
    opts: &SearchArgs,
    tx: mpsc::SyncSender<OutputChunk>,
    cancelled: &AtomicBool,
) -> Result<()> {
    (0..queries.len()).into_par_iter().try_for_each_init(
        || {
            (
                SearchState::<M>::new(&opts.score, &opts.extend),
                FormatState::new(),
                tx.clone(),
            )
        },
        |(state, format, tx), query_idx| {
            if cancelled.load(Ordering::Relaxed) {
                return Ok(());
            }

            let mut flush_to_tx = |chunk: OutputChunk| -> Result<()> {
                tx.send(chunk)
                    .map_err(|_| anyhow::anyhow!("Search hit channel disconnected"))
            };

            process_query::<M, _, _>(
                query_idx as u32,
                &queries.entries()[query_idx as usize],
                queries,
                store,
                global,
                opts,
                state,
                format,
                true,
                || cancelled.load(Ordering::Relaxed),
                &mut flush_to_tx,
            )?;

            Ok(())
        },
    )
}

/// Enumerate seeds for a single query across all targets via the global SA,
/// optionally extend them, and emit hits.
fn search_query<M, F>(
    query_idx: u32,
    query: &QueryData,
    global: &GlobalView<'_>,
    opts: &SearchArgs,
    state: &mut SearchState<M>,
    on_hit: &mut F,
) -> Result<()>
where
    M: DsmModel,
    F: FnMut(SearchHit) -> Result<()>,
{
    let query_bases = query.sequence().as_slice();
    let seed_interval = query.seed_interval();
    let include_alignment = opts.output.format != OutputFormat::Minimal;
    let pair_matrix = pair_mat(opts.seed.allows_wobble());
    let filter_cfg = &opts.filter;

    let mut callback_err: Option<anyhow::Error> = None;

    // Unified path: `compute_seed_extension` handles both extension and
    // max_extension==0 no-extension cases.
    for_each_seed(query, global, &opts.seed, |seed| {
        if callback_err.is_some() {
            return;
        }

        let target_idx_u32 = seed.target_id.0;
        let target_idx = target_idx_u32 as usize;
        let (t_forward_trans, t_reverse_trans, target_len) =
            target_transformed_slices(global, target_idx);

        let mut seed = seed;
        seed.target_start = target_len.saturating_sub(seed.target_start + seed.seed_len.get());

        let target_trans = match seed.strand {
            Strand::Forward => t_forward_trans,
            Strand::Reverse => t_reverse_trans,
        };

        let Some(hit) = build_hit_from_seed::<M>(
            &mut state.extender,
            state.dp_cfg,
            query_idx,
            query_bases,
            seed_interval,
            include_alignment,
            pair_matrix,
            filter_cfg,
            target_len,
            &seed,
            target_trans,
        ) else {
            return;
        };

        if let Err(err) = on_hit(hit) {
            callback_err = Some(err);
        }
    });

    if let Some(err) = callback_err {
        return Err(err);
    }

    Ok(())
}

/// Build a finalized `SearchHit` from a seed if extension and energy filters pass.
fn build_hit_from_seed<M: DsmModel>(
    extender: &mut DpExtender<M>,
    dp_cfg: DpConfig,
    query_idx: u32,
    query_bases: &[Base],
    seed_interval: Interval,
    include_alignment: bool,
    pair_matrix: &'static [[u8; 6]; 6],
    filter_cfg: &FilterConfig,
    target_len: usize,
    seed: &SeedHit,
    target_trans: &[Base],
) -> Option<SearchHit> {
    if seed.target_start + seed.seed_len.get() > target_trans.len() {
        return None;
    }

    let extension = compute_seed_extension::<M>(
        extender,
        dp_cfg,
        query_bases,
        target_trans,
        seed,
        seed_interval,
        filter_cfg,
        pair_matrix,
        include_alignment,
    )?;

    if extension.score > filter_cfg.delta_g {
        return None;
    }

    Some(SearchHit::new(
        query_idx,
        query_bases,
        target_trans,
        seed,
        &extension,
        include_alignment,
        target_len,
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

#[inline(always)]
fn seed_energy_transformed<M: DsmModel>(
    query_bases: &[Base],
    target_transformed: &[Base],
    q_pos: usize,
    t_pos: usize,
    len: usize,
    penalty: i32,
) -> i32 {
    if len <= 1 {
        return 0;
    }
    let mut score = 0;
    let t_match_end = t_pos + len - 1;
    for i in 0..(len - 1) {
        score += stack_with_penalty::<M>(
            query_bases[q_pos + i].idx(),
            query_bases[q_pos + i + 1].idx(),
            target_transformed[t_match_end - i].complement().idx(),
            target_transformed[t_match_end - i - 1].complement().idx(),
            penalty,
        );
    }
    score
}

/// Compute optional left/right DP extension around a seed and return extension metadata.
fn compute_seed_extension<M: DsmModel>(
    extender: &mut DpExtender<M>,
    dp_cfg: DpConfig,
    query_bases: &[Base],
    target_trans: &[Base],
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

    // Maximality check in scoring/raw coordinate space.
    if q_pos > seed_interval.start && t_pos + len < target_trans.len() {
        let p_class =
            pair_matrix[query_bases[q_pos - 1].idx()][target_trans[t_pos + len].complement().idx()];
        if p_class != 0 && !filter_cfg.no_max_prune {
            return None;
        }
    }

    if q_pos + len < seed_interval.end && t_pos > 0 {
        let p_class =
            pair_matrix[query_bases[q_pos + len].idx()][target_trans[t_pos - 1].complement().idx()];
        if p_class != 0 && !filter_cfg.no_max_prune {
            return None;
        }
    }

    let t_match_end = t_pos + len - 1;
    let max_ext = dp_cfg.max_extension();
    let seed_e =
        seed_energy_transformed::<M>(query_bases, target_trans, q_pos, t_pos, len, penalty);

    let can_extend_left = q_pos > 0 && t_pos + len < target_trans.len();
    let can_extend_right = q_pos + len < query_bases.len() && t_pos > 0;

    if max_ext == 0 || (!can_extend_left && !can_extend_right) {
        let term_5p = terminal_5p::<M>(
            query_bases[q_pos],
            target_trans[t_match_end].complement(),
            penalty,
        );
        let term_3p = terminal_3p::<M>(
            query_bases[q_pos + len - 1],
            target_trans[t_pos].complement(),
            penalty,
        );
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
        let view = DpView::<M>::left(query_bases, target_trans, q_pos, t_match_end, max_ext);
        let result = extender.extend(&view);
        let pairs = if with_traceback {
            result.traceback(&view)
        } else {
            SmallVec::new()
        };
        (result.score, result.q_len, result.t_len, pairs)
    };

    let (r_score, r_q, r_t, right_pairs) = {
        let view = DpView::<M>::right(query_bases, target_trans, q_pos + len - 1, t_pos, max_ext);
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
        target_trans: &[Base],
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
                    target_trans[t_match_end - i].complement(),
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
