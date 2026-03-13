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
use std::ops::Range;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use crate::alignment::{Alignment, PairClass};
use crate::config::{FilterConfig, OutputFormat, SearchConfig};
use crate::dp::gotoh::Gotoh;
use crate::dp::{DpConfig, DpGrid, DpView};
use crate::dsm::ScoringTable;
use crate::index::store::{GlobalView, TargetStore};
use crate::output::format::HitCtx;
use crate::output::writer::{HitFormatter, OutputChunk, OutputWriter};
use crate::registry::QueryRegistry;
use crate::seed::{for_each_seed, SeedHit};
use crate::types::{Base, Energy, Strand};

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
    opts: &SearchConfig,
    output_path: &Path,
) -> Result<usize> {
    if store.is_empty() || queries.is_empty() {
        return Ok(0);
    }
    info!(
        "Starting search: {} queries x {} targets, seed={:?}, max_ext={}, delta_g={}",
        queries.len(),
        store.len(),
        opts.seed.seed,
        opts.extend.max_extension,
        opts.filter.delta_g
    );
    let ctx = SearchContext::new(queries, store, opts);

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
    opts: &'a SearchConfig,
}

impl<'a> SearchContext<'a> {
    fn new(
        queries: &'a QueryRegistry,
        store: &'a TargetStore,
        opts: &'a SearchConfig,
    ) -> Self {
        let global = store.global_view();
        Self { queries, store, global, opts }
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

/// Per-worker DP extension engine. Each worker owns one; `&mut` is safe because
/// only `grid` is mutated — the rest is read-only configuration.
struct ExtensionEngine {
    grid: DpGrid,
    dp_cfg: DpConfig,
    model: ScoringTable,
    gotoh_left: Gotoh,
    gotoh_right: Gotoh,
}

impl ExtensionEngine {
    fn new(opts: &SearchConfig) -> Self {
        let dp_cfg = DpConfig::from((&opts.score, &opts.extend));
        let model = ScoringTable::new(
            opts.score.matrix,
            dp_cfg.penalty_raw(),
            opts.seed.allows_wobble(),
        );
        let left_model = model.transpose();
        let gotoh_right = Gotoh::new(&model);
        let gotoh_left = Gotoh::new(&left_model);
        Self { grid: DpGrid::new(dp_cfg.max_extension()), dp_cfg, model, gotoh_left, gotoh_right }
    }
}

/// Single-file backend: `Mutex<OutputWriter>` + rayon `par_iter`.
fn run_single_file(ctx: &SearchContext<'_>, output_path: &Path) -> Result<usize> {
    let writer = Mutex::new(OutputWriter::new(&ctx.opts.output, output_path)?);
    let total = AtomicUsize::new(0);

    (0..ctx.queries.len()).into_par_iter().try_for_each_init(
        || {
            (
                ExtensionEngine::new(ctx.opts),
                HitFormatter::new(ctx.opts.output.format),
            )
        },
        |(engine, fmt), qi| -> Result<()> {
            let emitted = process_query(ctx, qi as u32, engine, fmt, &mut |chunk| {
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
                ExtensionEngine::new(ctx.opts),
                HitFormatter::new(ctx.opts.output.format),
            )
        },
        |(engine, format), query_idx| -> Result<()> {
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

            let emitted = process_query(ctx, query_idx as u32, engine, format, &mut flush_to_writer)?;

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

/// Process one query end-to-end: enumerate seeds, extend, filter, format, and
/// stream chunks via `on_chunk`.
fn process_query<FO>(
    ctx: &SearchContext<'_>,
    query_idx: u32,
    engine: &mut ExtensionEngine,
    format: &mut HitFormatter,
    on_chunk: &mut FO,
) -> Result<usize>
where
    FO: FnMut(OutputChunk) -> Result<()>,
{
    let query = &ctx.queries.entries()[query_idx as usize];
    let query_name = ctx.queries.get_name(query_idx);
    let query_seq = query.sequence().as_slice();
    let seed_interval = query.seed_interval.clone();
    let include_alignment = ctx.opts.output.format != OutputFormat::Minimal;
    let filter_cfg = &ctx.opts.filter;

    let mut local_hits = 0usize;
    let mut last_target_idx = None::<usize>;
    let mut cached_t_name = None;

    for_each_seed(query, &ctx.global, &ctx.opts.seed, |seed| {
        let target_idx = seed.target_id.0 as usize;
        let (t_fwd, t_rc, target_len) = ctx.target_slices(target_idx);
        let target_trans = match seed.strand {
            Strand::Forward => t_fwd,
            Strand::Reverse => t_rc,
        };

        let Some(hit) = build_hit_from_seed(
            engine,
            query_idx,
            query_seq,
            seed_interval.clone(),
            include_alignment,
            filter_cfg,
            target_len,
            &seed,
            target_trans,
        ) else {
            return Ok(());
        };

        if last_target_idx != Some(target_idx) {
            cached_t_name = Some(ctx.store.get_name(hit.target_idx));
            last_target_idx = Some(target_idx);
        }
        let hit_ctx = HitCtx {
            q_name: query_name,
            q_seq: query_seq,
            t_name: cached_t_name.unwrap(),
            t_fwd,
            t_rc,
        };

        if let Some(chunk) = format.add_hit(&hit, hit_ctx) {
            on_chunk(chunk)?;
        }
        local_hits += 1;
        Ok(())
    })?;

    if let Some(chunk) = format.flush() {
        on_chunk(chunk)?;
    }

    Ok(local_hits)
}

/// Build a finalized `SearchHit` from a seed if extension and energy filters pass.
fn build_hit_from_seed(
    engine: &mut ExtensionEngine,
    query_idx: u32,
    query_bases: &[Base],
    seed_interval: Range<usize>,
    include_alignment: bool,
    filter_cfg: &FilterConfig,
    target_len: usize,
    seed: &SeedHit,
    target_trans: &[Base],
) -> Option<SearchHit> {
    debug_assert!(seed.target_start + seed.len.get() <= target_trans.len());
    engine
        .extend_seed(query_bases, target_trans, seed, seed_interval, filter_cfg, include_alignment)
        .filter(|ext| ext.score <= filter_cfg.delta_g)
        .map(|ext| SearchHit::new(query_idx, query_bases, target_trans, seed, &ext, target_len))
}

// =============================================================================
// EXTENSION
// =============================================================================

type AlignmentPairs = (SmallVec<[PairClass; 64]>, SmallVec<[PairClass; 64]>);

struct SeedExtension {
    score: f64,
    l_q: usize,
    l_t: usize,
    r_q: usize,
    r_t: usize,
    pairs: Option<AlignmentPairs>,
}

impl ExtensionEngine {
    fn extend_seed(
        &mut self,
        query_bases: &[Base],
        target_trans: &[Base],
        seed: &SeedHit,
        seed_interval: Range<usize>,
        filter_cfg: &FilterConfig,
        include_alignment: bool,
    ) -> Option<SeedExtension> {
        let penalty = self.dp_cfg.penalty_raw();
        let q_start = seed.query_start;
        let t_start = seed.target_start;
        let len = seed.len.get();

        // Maximality check in scoring/raw coordinate space.
        if q_start > seed_interval.start && t_start + len < target_trans.len() {
            let q_base = query_bases[q_start - 1];
            let t_base = target_trans[t_start + len];
            if self.model.is_pair(q_base, t_base) && !filter_cfg.no_max_prune {
                return None;
            }
        }

        if q_start + len < seed_interval.end && t_start > 0 {
            let q_base = query_bases[q_start + len];
            let t_base = target_trans[t_start - 1];
            if self.model.is_pair(q_base, t_base) && !filter_cfg.no_max_prune {
                return None;
            }
        }

        let t_match_end = t_start + len - 1;
        let max_ext = self.dp_cfg.max_extension();
        let seed_e = crate::dsm::seed_energy(&self.model, query_bases, target_trans, q_start, t_start, len);

        let can_extend_left = q_start > 0 && t_start + len < target_trans.len();
        let can_extend_right = q_start + len < query_bases.len() && t_start > 0;

        if max_ext == 0 || (!can_extend_left && !can_extend_right) {
            let term_5p = self.gotoh_left.terminal(
                query_bases[q_start].idx(),
                target_trans[t_match_end].idx(),
            );
            let term_3p = self.gotoh_right.terminal(
                query_bases[q_start + len - 1].idx(),
                target_trans[t_start].idx(),
            );
            let nt_count = (2 * len) as i32;
            return Some(SeedExtension {
                score: crate::dsm::to_kcal(seed_e + term_5p + term_3p + nt_count * penalty),
                l_q: 0,
                l_t: 0,
                r_q: 0,
                r_t: 0,
                pairs: include_alignment.then(|| (SmallVec::new(), SmallVec::new())),
            });
        }

        let (l_score, l_q, l_t, l_pairs) = {
            let view = DpView::left(query_bases, target_trans, q_start, t_match_end, max_ext);
            let result = self.gotoh_left.extend(&view, &mut self.grid);
            let pairs = include_alignment.then(|| result.traceback(&view));
            (result.score, result.q_len, result.t_len, pairs)
        };

        let (r_score, r_q, r_t, r_pairs) = {
            let view = DpView::right(query_bases, target_trans, q_start + len - 1, t_start, max_ext);
            let result = self.gotoh_right.extend(&view, &mut self.grid);
            let pairs = include_alignment.then(|| result.traceback(&view));
            (result.score, result.q_len, result.t_len, pairs)
        };

        let nt_count = (l_q + l_t + r_q + r_t + 2 * len) as i32;
        Some(SeedExtension {
            score: crate::dsm::to_kcal(seed_e + l_score + r_score + nt_count * penalty),
            l_q,
            l_t,
            r_q,
            r_t,
            pairs: l_pairs.zip(r_pairs),
        })
    }
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
        original_target_len: usize,
    ) -> Self {
        let q_start = seed.query_start;
        let t_start = seed.target_start;
        let len = seed.len.get();

        let final_q_start = q_start.saturating_sub(ext.l_q);
        let final_q_end = (q_start + len - 1) + ext.r_q;
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

        let (alignment, seed_start, seed_end) = if let Some((left, right)) = &ext.pairs {
            let t_match_end = seed.target_start + len - 1;
            let mut seed_pairs: SmallVec<[PairClass; 64]> = SmallVec::with_capacity(len);
            for i in 0..len {
                seed_pairs.push(PairClass::from_bases(
                    query_bases[q_start + i],
                    target_trans[t_match_end - i].complement(),
                ));
            }
            let start = left.len();
            (Some(Alignment::new(left, &seed_pairs, right)), Some(start), Some(start + len))
        } else {
            (None, None, None)
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
        }
    }
}
