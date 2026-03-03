//! Search module - finds miRNA-target interactions.
//!
//! Pipeline:
//! 1) Parallel iteration over queries (global SA traversal per query)
//! 2) Seed enumeration across all targets via single SA traversal
//! 3) Optional extension of each seed
//! 4) Materialize and emit final hits

use anyhow::{Context, Result};
use log::info;
use rayon::prelude::*;
use smallvec::SmallVec;
use std::io::Write;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use crate::alignment::{Alignment, PairClass};
use crate::config::{ExtendConfig, FilterConfig, OutputFormat, ScoreConfig, SearchArgs};
use crate::dp::{DpConfig, DpExtender, DpView};
use crate::dsm::{pair_mat, DsmModel};
use crate::index::store::{GlobalView, TargetStore};
use crate::output::writer::{HitFormatter, OutputChunk, OutputWriter};
use crate::registry::QueryRegistry;
use crate::seed::{for_each_seed, SeedHit};
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

// =============================================================================
// ORCHESTRATION
// =============================================================================

/// Run search against mmap-backed target store and write hits to `output_path`.
///
/// In multifile mode, `output_path` is the directory where per-query files are
/// created. Otherwise it is the single output file path.
pub fn run_search(
    queries: &QueryRegistry,
    store: &TargetStore,
    opts: &SearchArgs,
    output_path: &Path,
) -> Result<usize> {
    let Some(ctx) = SearchContext::try_new(queries, store, opts) else {
        return Ok(0);
    };

    let total = if opts.output.multifile {
        std::fs::create_dir_all(output_path).with_context(|| {
            format!(
                "Failed to create output directory {:?} for --multifile",
                output_path
            )
        })?;
        run_multifile(&ctx, output_path)
    } else {
        run_single_file(&ctx, output_path)
    }?;

    info!("Search complete: {} hits", total);
    Ok(total)
}

/// Shared context for all search backends.
struct SearchContext<'a> {
    queries: &'a QueryRegistry,
    store: &'a TargetStore,
    global: GlobalView<'a>,
    opts: &'a SearchArgs,
    model: DsmModel,
}

impl<'a> SearchContext<'a> {
    fn try_new(
        queries: &'a QueryRegistry,
        store: &'a TargetStore,
        opts: &'a SearchArgs,
    ) -> Option<Self> {
        if store.is_empty() || queries.is_empty() {
            return None;
        }

        let global = store.global_view();
        let dp_cfg = DpConfig::from((&opts.score, &opts.extend));
        let model = DsmModel::new(opts.score.matrix, dp_cfg.penalty_raw());

        info!(
            "Starting search: {} queries x {} targets, seed={:?}, max_ext={}, delta_g={}",
            queries.len(),
            store.len(),
            opts.seed.seed,
            opts.extend.max_extension,
            opts.filter.delta_g
        );

        Some(Self {
            queries,
            store,
            global,
            opts,
            model,
        })
    }

    #[inline]
    fn target_slices(&self, target_idx: usize) -> (&[Base], &[Base], usize) {
        let target_len = self.global.seq_lens[target_idx] as usize;
        let target_offset = self.global.offsets[target_idx] as usize;
        let t_fwd = &self.global.combined_seq[target_offset..target_offset + target_len];
        let t_rc = &self.global.combined_seq
            [target_offset + target_len + 1..target_offset + 2 * target_len + 1];
        (t_fwd, t_rc, target_len)
    }
}

/// Reusable per-worker state.
struct SearchState {
    extender: DpExtender,
    dp_cfg: DpConfig,
}

impl SearchState {
    fn new(model: DsmModel, score_cfg: &ScoreConfig, extend_cfg: &ExtendConfig) -> Self {
        let dp_cfg = DpConfig::from((score_cfg, extend_cfg));
        Self {
            extender: DpExtender::from_config(model, dp_cfg),
            dp_cfg,
        }
    }
}

/// Single-file backend: `Mutex<OutputWriter>` + rayon `par_iter`.
fn run_single_file(ctx: &SearchContext<'_>, output_path: &Path) -> Result<usize> {
    let writer = Mutex::new(OutputWriter::new(&ctx.opts.output, output_path)?);
    let total = AtomicUsize::new(0);

    (0..ctx.queries.len()).into_par_iter().try_for_each_init(
        || {
            (
                SearchState::new(ctx.model.clone(), &ctx.opts.score, &ctx.opts.extend),
                HitFormatter::new(ctx.opts.output.format),
            )
        },
        |(state, fmt), qi| -> Result<()> {
            let emitted =
                process_query(ctx, qi as u32, state, fmt, true, || false, &mut |chunk| {
                    writer.lock().unwrap().write_chunk(&chunk)
                })?;
            total.fetch_add(emitted, Ordering::Relaxed);
            Ok(())
        },
    )?;

    writer.into_inner().unwrap().flush_all()?;
    Ok(total.load(Ordering::Relaxed))
}

/// Multifile backend: per-worker lazy file writers + rayon `par_iter`.
fn run_multifile(ctx: &SearchContext<'_>, output_dir: &Path) -> Result<usize> {
    let total = AtomicUsize::new(0);
    let ext = crate::output::output_extension(&ctx.opts.output);
    let output_paths = crate::output::writer::build_multifile_paths(ctx.queries, output_dir, ext);

    (0..ctx.queries.len()).into_par_iter().try_for_each_init(
        || {
            (
                SearchState::new(ctx.model.clone(), &ctx.opts.score, &ctx.opts.extend),
                HitFormatter::new(ctx.opts.output.format),
            )
        },
        |(state, format), query_idx| -> Result<()> {
            let file_path = &output_paths[query_idx];
            let mut writer: Option<Box<dyn Write>> = None;

            let mut flush_to_writer = |chunk: OutputChunk| -> Result<()> {
                if writer.is_none() {
                    writer = Some(
                        crate::output::open_output(Some(file_path), &ctx.opts.output)
                            .with_context(|| {
                                format!("Failed to open output file {:?}", file_path)
                            })?,
                    );
                }
                writer
                    .as_mut()
                    .expect("writer initialized")
                    .write_all(&chunk.data)
                    .with_context(|| format!("Failed to write output file {:?}", file_path))
            };

            let emitted = process_query(
                ctx,
                query_idx as u32,
                state,
                format,
                true,
                || false,
                &mut flush_to_writer,
            )?;

            if let Some(writer) = writer.as_mut() {
                writer
                    .flush()
                    .with_context(|| format!("Failed to flush output file {:?}", file_path))?;
            }

            total.fetch_add(emitted, Ordering::Relaxed);
            Ok(())
        },
    )?;

    Ok(total.load(Ordering::Relaxed))
}

/// Process one query end-to-end and stream formatted chunks via `on_chunk`.
fn process_query<FC, FO>(
    ctx: &SearchContext<'_>,
    query_idx: u32,
    state: &mut SearchState,
    format: &mut HitFormatter,
    flush_after_query: bool,
    is_cancelled: FC,
    on_chunk: &mut FO,
) -> Result<usize>
where
    FC: Fn() -> bool,
    FO: FnMut(OutputChunk) -> Result<()>,
{
    if is_cancelled() {
        return Ok(0);
    }

    let query_name = ctx.queries.get_name(query_idx);
    let query = &ctx.queries.entries()[query_idx as usize];
    let query_seq = query.sequence().as_slice();
    let mut local_hits = 0usize;
    let mut last_target_idx = None;
    let mut cached_target = None;

    search_query(ctx, query_idx, state, &mut |hit: SearchHit| -> Result<()> {
        if is_cancelled() {
            return Ok(());
        }

        let target_idx = hit.target_idx as usize;
        if last_target_idx != Some(target_idx) {
            let (t_fwd, t_rc, _) = ctx.target_slices(target_idx);
            let target_name = ctx.store.get_name(hit.target_idx);
            cached_target = Some((target_name, t_fwd, t_rc));
            last_target_idx = Some(target_idx);
        }
        let (target_name, t_fwd, t_rc) = cached_target.as_ref().unwrap();

        if let Some(chunk) = format.add_hit(&hit, query_name, target_name, query_seq, t_fwd, t_rc) {
            on_chunk(chunk)?;
        }
        local_hits += 1;
        Ok(())
    })?;

    if flush_after_query {
        if let Some(chunk) = format.flush() {
            on_chunk(chunk)?;
        }
    }

    Ok(local_hits)
}

/// Enumerate seeds for a single query across all targets via the global SA,
/// optionally extend them, and emit hits.
fn search_query<F>(
    ctx: &SearchContext<'_>,
    query_idx: u32,
    state: &mut SearchState,
    on_hit: &mut F,
) -> Result<()>
where
    F: FnMut(SearchHit) -> Result<()>,
{
    let query = &ctx.queries.entries()[query_idx as usize];
    let query_bases = query.sequence().as_slice();
    let seed_interval = query.seed_interval();
    let include_alignment = ctx.opts.output.format != OutputFormat::Minimal;
    let pair_matrix = pair_mat(ctx.opts.seed.allows_wobble());
    let filter_cfg = &ctx.opts.filter;

    let mut callback_err: Option<anyhow::Error> = None;

    // Unified path: `compute_seed_extension` handles both extension and
    // max_extension==0 no-extension cases.
    for_each_seed(query, &ctx.global, &ctx.opts.seed, |seed| {
        if callback_err.is_some() {
            return;
        }

        let target_idx_u32 = seed.target_id.0;
        let target_idx = target_idx_u32 as usize;
        let (t_forward_trans, t_reverse_trans, target_len) = ctx.target_slices(target_idx);

        let mut seed = seed;
        seed.target_start = target_len.saturating_sub(seed.target_start + seed.seed_len.get());

        let target_trans = match seed.strand {
            Strand::Forward => t_forward_trans,
            Strand::Reverse => t_reverse_trans,
        };

        let Some(hit) = build_hit_from_seed(
            state,
            &ctx.model,
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
#[allow(clippy::too_many_arguments)]
fn build_hit_from_seed(
    state: &mut SearchState,
    model: &DsmModel,
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

    let extension = compute_seed_extension(
        &mut state.extender,
        model,
        state.dp_cfg,
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

/// Compute optional left/right DP extension around a seed and return extension metadata.
#[allow(clippy::too_many_arguments)]
fn compute_seed_extension(
    extender: &mut DpExtender,
    model: &DsmModel,
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
    let seed_e = model.seed_energy(query_bases, target_trans, q_pos, t_pos, len);

    let can_extend_left = q_pos > 0 && t_pos + len < target_trans.len();
    let can_extend_right = q_pos + len < query_bases.len() && t_pos > 0;

    if max_ext == 0 || (!can_extend_left && !can_extend_right) {
        let term_5p = model.terminal_5p(query_bases[q_pos].idx(), target_trans[t_match_end].complement().idx());
        let term_3p = model.terminal_3p(query_bases[q_pos + len - 1].idx(), target_trans[t_pos].complement().idx());
        let nt_count = (2 * len) as i32;
        return Some(SeedExtension {
            score: crate::dsm::to_kcal(seed_e + term_5p + term_3p + nt_count * penalty),
            l_q: 0,
            l_t: 0,
            r_q: 0,
            r_t: 0,
            left_pairs: SmallVec::new(),
            right_pairs: SmallVec::new(),
        });
    }

    let (l_score, l_q, l_t, left_pairs) = {
        let view = DpView::left(query_bases, target_trans, q_pos, t_match_end, max_ext);
        let result = extender.extend(&view);
        let pairs = if with_traceback {
            result.traceback(&view)
        } else {
            SmallVec::new()
        };
        (result.score, result.q_len, result.t_len, pairs)
    };

    let (r_score, r_q, r_t, right_pairs) = {
        let view = DpView::right(query_bases, target_trans, q_pos + len - 1, t_pos, max_ext);
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
        score: crate::dsm::to_kcal(seed_e + l_score + r_score + nt_count * penalty),
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
