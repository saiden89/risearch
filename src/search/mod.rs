//! Search module - finds miRNA-target interactions.
//!
//! Pipeline:
//! 1) Parallel iteration over queries (global SA traversal per query)
//! 2) Seed enumeration across all targets via single SA traversal
//! 3) Optional extension of each seed
//! 4) Materialize and emit final hits

mod extension;

use anyhow::{Context, Result};
use log::info;
use rayon::prelude::*;
use smallvec::SmallVec;
use std::ops::Range;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use self::extension::ExtensionEngine;
use crate::alignment::{Alignment, PairClass};
use crate::config::{OutputFormat, SearchConfig};
use crate::dp::{DpConfig, DpView, ExtendDir};
use crate::dsm::{DsmRegistry, ScoringModel};
use crate::index::store::TargetStore;
use crate::output::writer::{build_multifile_paths, HitFormatter, OutputChunk, OutputWriter};
use crate::registry::QueryRegistry;
use crate::seed::{SeedHit, SeedingEngine};
use crate::types::{Base, Energy, Strand};


#[derive(Debug, Clone)]
pub struct SearchHit {
    pub query_idx: usize,
    pub target_idx: usize,
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

/// Run search collecting all hits into memory. Alignment data is always included.
pub fn run_search_in_memory(
    queries: &QueryRegistry,
    store: &TargetStore,
    opts: &SearchConfig,
) -> Result<Vec<SearchHit>> {
    let Some(ctx) = init_search(queries, store, opts)? else {
        return Ok(Vec::new());
    };
    let seeds = SeedingEngine::new(ctx.queries, ctx.store).run(&ctx.opts.seed);

    let hits: Vec<SearchHit> = seeds
        .into_par_iter()
        .flat_map_iter(|(qi, seeds)| {
            let mut worker = SearchWorker::new(ctx.opts, &ctx.model);
            let query = &ctx.queries.entries()[qi];
            let query_seq = query.sequence();
            let seed_interval = query.seed_interval.clone();
            seeds.into_iter().filter_map(move |seed| {
                let target_idx = seed.target_idx;
                let (t_fwd, t_rc, target_len) = ctx.store.target_slices(target_idx);
                let target_trans = match seed.strand {
                    Strand::Forward => t_fwd,
                    Strand::Reverse => t_rc,
                };

                worker.build_hit_from_seed(
                    ctx.opts,
                    qi,
                    query_seq,
                    seed_interval.clone(),
                    true,
                    target_len,
                    &seed,
                    target_trans,
                )
            })
        })
        .collect();

    Ok(hits)
}

/// Run search and write hits to `output_path` (or directory in multifile mode).
pub fn run_search(
    queries: &QueryRegistry,
    store: &TargetStore,
    opts: &SearchConfig,
    output_path: &Path,
) -> Result<()> {
    let Some(ctx) = init_search(queries, store, opts)? else {
        return Ok(());
    };
    let seed_groups = SeedingEngine::new(ctx.queries, ctx.store).run(&ctx.opts.seed);
    let total = AtomicUsize::new(0);

    if opts.output.multifile {
        std::fs::create_dir_all(output_path).with_context(|| {
            format!(
                "Failed to create output directory {:?} for --multifile",
                output_path
            )
        })?;
        let ext = opts.output.compress.extension();
        let paths = build_multifile_paths(queries, output_path, ext);

        seed_groups.into_par_iter().try_for_each_init(
            || {
                (
                    SearchWorker::new(ctx.opts, &ctx.model),
                    HitFormatter::new(ctx.opts.output.format),
                )
            },
            |(worker, fmt), (qi, seeds)| -> Result<()> {
                let (emitted, chunks) = worker.process_query_seeds(&ctx, qi, &seeds, fmt)?;
                if !chunks.is_empty() {
                    let mut w = OutputWriter::new(&ctx.opts.output, &paths[qi])?;
                    write_chunks(&mut w, &chunks)?;
                    w.flush_all()?;
                }
                total.fetch_add(emitted, Ordering::Relaxed);
                Ok(())
            },
        )?;
    } else {
        let writer = Mutex::new(OutputWriter::new(&opts.output, output_path)?);

        seed_groups.into_par_iter().try_for_each_init(
            || {
                (
                    SearchWorker::new(ctx.opts, &ctx.model),
                    HitFormatter::new(ctx.opts.output.format),
                )
            },
            |(worker, fmt), (qi, seeds)| -> Result<()> {
                let (emitted, chunks) = worker.process_query_seeds(&ctx, qi, &seeds, fmt)?;
                write_chunks(&mut *writer.lock().unwrap(), &chunks)?;
                total.fetch_add(emitted, Ordering::Relaxed);
                Ok(())
            },
        )?;

        writer.into_inner().unwrap().flush_all()?;
    }

    let total = total.load(Ordering::Relaxed);
    info!("Search complete: {} hits", total);
    Ok(())
}

struct SearchContext<'a> {
    queries: &'a QueryRegistry,
    store: &'a TargetStore,
    opts: &'a SearchConfig,
    model: ScoringModel,
}

impl<'a> SearchContext<'a> {
    fn new(
        queries: &'a QueryRegistry,
        store: &'a TargetStore,
        opts: &'a SearchConfig,
        model: ScoringModel,
    ) -> Self {
        Self {
            queries,
            store,
            opts,
            model,
        }
    }
}

fn init_search<'a>(
    queries: &'a QueryRegistry,
    store: &'a TargetStore,
    opts: &'a SearchConfig,
) -> Result<Option<SearchContext<'a>>> {
    if store.is_empty() || queries.is_empty() {
        return Ok(None);
    }

    info!(
        "Starting search: {} queries x {} targets, seed_length={:?}, max_ext={}, delta_g={}",
        queries.len(),
        store.len(),
        opts.seed.seed_length,
        opts.extend.max_extension,
        opts.filter.delta_g
    );

    let (initiation, source_table) =
        DsmRegistry::load(&opts.score.dsm_id, opts.score.temperature)?;

    let model = ScoringModel::new(&source_table, initiation, opts.score.penalty);
    Ok(Some(SearchContext::new(queries, store, opts, model)))
}

struct SearchWorker {
    extension: ExtensionEngine,
    model: ScoringModel,
}

impl SearchWorker {
    fn new(opts: &SearchConfig, model: &ScoringModel) -> Self {
        let dp_cfg = DpConfig::from(&opts.extend);
        let model = model.clone();
        Self {
            extension: ExtensionEngine::new(dp_cfg.max_extension(), &model),
            model,
        }
    }

    fn process_query_seeds(
        &mut self,
        ctx: &SearchContext<'_>,
        query_idx: usize,
        seeds: &[SeedHit],
        format: &mut HitFormatter,
    ) -> Result<(usize, Vec<OutputChunk>)> {
        let query = &ctx.queries.entries()[query_idx];
        let query_name = ctx.queries.get_name(query_idx);
        let query_seq = query.sequence();
        let seed_interval = query.seed_interval.clone();
        let include_alignment = ctx.opts.output.format != OutputFormat::Minimal;
        let mut chunks = Vec::new();
        let mut local_hits = 0usize;
        let mut last_target_idx = None::<usize>;
        let mut cached_t_name = None;

        for seed in seeds {
            let target_idx = seed.target_idx;
            let (t_fwd, t_rc, target_len) = ctx.store.target_slices(target_idx);
            let target_trans = match seed.strand {
                Strand::Forward => t_fwd,
                Strand::Reverse => t_rc,
            };

            let Some(hit) = self.build_hit_from_seed(
                ctx.opts,
                query_idx,
                query_seq,
                seed_interval.clone(),
                include_alignment,
                target_len,
                seed,
                target_trans,
            ) else {
                continue;
            };

            if last_target_idx != Some(target_idx) {
                cached_t_name = Some(ctx.store.get_name(hit.target_idx));
                last_target_idx = Some(target_idx);
            }

            if let Some(chunk) = format.add_hit(
                &hit,
                query_name,
                query_seq,
                cached_t_name.unwrap(),
                t_fwd,
                t_rc,
            ) {
                chunks.push(chunk);
            }
            local_hits += 1;
        }

        if let Some(chunk) = format.flush() {
            chunks.push(chunk);
        }

        Ok((local_hits, chunks))
    }

    fn build_hit_from_seed(
        &mut self,
        opts: &SearchConfig,
        query_idx: usize,
        query_bases: &[Base],
        seed_interval: Range<usize>,
        include_alignment: bool,
        target_len: usize,
        seed: &SeedHit,
        target_trans: &[Base],
    ) -> Option<SearchHit> {
        debug_assert!(seed.target_start + seed.len <= target_trans.len());
        let q_start = seed.query_start;
        let t_start = seed.target_start;
        let len = seed.len;
        let t_match_end = t_start + len - 1;

        if !opts.filter.no_max_prune
            && !is_maximal(
                seed,
                query_bases,
                target_trans,
                &seed_interval,
                opts.seed.seed_wobble,
            )
        {
            return None;
        }

        let duplex_score = self
            .model
            .ungapped_duplex_score(query_bases, target_trans, q_start, t_start, len);
        let max_ext = DpConfig::from(&opts.extend).max_extension();
        let view_left = DpView::new(
            query_bases,
            target_trans,
            q_start,
            t_match_end,
            ExtendDir::Left,
            max_ext,
        );
        let left = self.extension.extend(&view_left, include_alignment);

        let view_right = DpView::new(
            query_bases,
            target_trans,
            q_start + len - 1,
            t_start,
            ExtendDir::Right,
            max_ext,
        );
        let right = self.extension.extend(&view_right, include_alignment);
        let nt_count = left.q_ext + left.t_ext + right.q_ext + right.t_ext + 2 * len;
        let stacking_stability = duplex_score + left.energy + right.energy;
        let energy = self.model.binding_energy(stacking_stability, nt_count);

        (energy <= opts.filter.delta_g).then(|| {
            SearchHit::new(
                query_idx,
                query_bases,
                target_trans,
                seed,
                &left,
                &right,
                energy,
                include_alignment,
                target_len,
            )
        })
    }
}

fn write_chunks(writer: &mut OutputWriter, chunks: &[OutputChunk]) -> Result<()> {
    for chunk in chunks {
        writer.write_chunk(chunk)?;
    }
    Ok(())
}

struct ExtensionResult {
    energy: Energy,
    q_ext: usize,
    t_ext: usize,
    pairs: Option<SmallVec<[PairClass; 64]>>,
}

fn is_maximal(
    seed: &SeedHit,
    query_bases: &[Base],
    target_trans: &[Base],
    seed_interval: &Range<usize>,
    seed_wobble: bool,
) -> bool {
    let q_start = seed.query_start;
    let t_start = seed.target_start;
    let len = seed.len;

    if q_start > seed_interval.start
        && t_start + len < target_trans.len()
        && query_bases[q_start - 1].forms_pair(target_trans[t_start + len], seed_wobble)
    {
        return false;
    }

    if q_start + len < seed_interval.end
        && t_start > 0
        && query_bases[q_start + len].forms_pair(target_trans[t_start - 1], seed_wobble)
    {
        return false;
    }

    true
}

const BINDING_SITE_FLANK_LEN: usize = 20;

impl SearchHit {
    pub(crate) fn query_bases<'a>(&self, q_seq: &'a [Base]) -> &'a [Base] {
        &q_seq[self.q_start..self.q_end + 1]
    }

    pub(crate) fn target_bases<'a>(&self, t_fwd: &'a [Base], t_rc: &'a [Base]) -> &'a [Base] {
        match self.strand {
            Strand::Forward => &t_fwd[self.t_start..self.t_end + 1],
            Strand::Reverse => {
                let len = t_fwd.len();
                &t_rc[len - 1 - self.t_end..len - self.t_start]
            }
        }
    }

    pub(crate) fn target_flanks<'a>(
        &self,
        t_fwd: &'a [Base],
        t_rc: &'a [Base],
    ) -> (&'a [Base], &'a [Base]) {
        let len = t_fwd.len();
        let (oriented, start, end) = match self.strand {
            Strand::Forward => (t_fwd, self.t_start, self.t_end),
            Strand::Reverse => (t_rc, len - 1 - self.t_end, len - 1 - self.t_start),
        };

        let right_start = end + 1;
        let right_end = (right_start + BINDING_SITE_FLANK_LEN).min(oriented.len());
        let left_start = start.saturating_sub(BINDING_SITE_FLANK_LEN);

        (
            &oriented[right_start..right_end],
            &oriented[left_start..start],
        )
    }

    fn new(
        query_idx: usize,
        query_bases: &[Base],
        target_trans: &[Base],
        seed: &SeedHit,
        left: &ExtensionResult,
        right: &ExtensionResult,
        energy: Energy,
        include_alignment: bool,
        original_target_len: usize,
    ) -> Self {
        let q_start = seed.query_start;
        let t_start = seed.target_start;
        let len = seed.len;

        let final_q_start = q_start - left.q_ext;
        let final_q_end = q_start + len - 1 + right.q_ext;
        let mut final_t_start = t_start - right.t_ext;
        let mut final_t_end = t_start + len - 1 + left.t_ext;

        if seed.strand == Strand::Reverse {
            let tmp = original_target_len - 1 - final_t_end;
            final_t_end = original_target_len - 1 - final_t_start;
            final_t_start = tmp;
        }

        let (alignment, seed_start, seed_end) = if include_alignment {
            let t_match_end = t_start + len - 1;
            let mut seed_pairs: SmallVec<[PairClass; 64]> = SmallVec::with_capacity(len);
            for i in 0..len {
                seed_pairs.push(PairClass::from_view_bases(
                    query_bases[q_start + i],
                    target_trans[t_match_end - i],
                ));
            }
            let left_pairs = left.pairs.as_deref().unwrap_or(&[]);
            let right_pairs = right.pairs.as_deref().unwrap_or(&[]);
            let start = left_pairs.len();
            (
                Some(Alignment::from_parts(
                    left_pairs,
                    &seed_pairs,
                    right_pairs,
                )),
                Some(start),
                Some(start + len),
            )
        } else {
            (None, None, None)
        };

        Self {
            query_idx,
            target_idx: seed.target_idx,
            q_start: final_q_start,
            q_end: final_q_end,
            t_start: final_t_start,
            t_end: final_t_end,
            strand: seed.strand,
            energy,
            seed_start,
            seed_end,
            alignment,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;
    use std::path::PathBuf;

    use super::*;
    use crate::config::{
        ExtendConfig, FilterConfig, OutputCompression, OutputConfig, OutputFormat,
        ScoreConfig, SeedConfig,
    };
    use crate::index::store::TargetStore;
    use crate::registry::QueryRegistry;
    use crate::types::DsmId;

    fn workspace_root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
    }

    fn test_config() -> SearchConfig {
        SearchConfig {
            seed: SeedConfig {
                seed_start: None,
                seed_end: None,
                seed_length: Some(8),
                seed_wobble: true,
                max_mismatches: 0,
                min_prefix_matches: 1,
                min_suffix_matches: 0,
            },
            score: ScoreConfig {
                dsm_id: DsmId::from("t04"),
                penalty: Energy::from_kcal(3.5),
                temperature: 37,
            },
            extend: ExtendConfig {
                max_extension: 20,
            },
            filter: FilterConfig {
                delta_g: Energy::from_kcal(-10.0),
                seed_energy: Energy::from_kcal(0.0),
                no_max_prune: false,
            },
            output: OutputConfig {
                format: OutputFormat::Detailed,
                compress: OutputCompression::None,
                multifile: false,
            },
        }
    }

    fn build_store(target_fa: &std::path::Path) -> (TargetStore, tempfile::TempDir) {
        let tmpdir = tempfile::tempdir().unwrap();
        let idx = tmpdir.path().join("target.idx");
        TargetStore::build(target_fa, &idx).unwrap();
        let store = TargetStore::open(&idx).unwrap();
        (store, tmpdir)
    }

    #[test]
    fn in_memory_matches_file_hit_count() {
        let root = workspace_root();
        let query_path = root.join("tests/data/query.fa");
        let target_path = root.join("tests/data/target.fa");

        let (store, _tmp) = build_store(&target_path);
        let mut config = test_config();
        config.output.format = OutputFormat::Minimal;
        let queries = QueryRegistry::from_fasta(&query_path, &config.seed).unwrap();

        let hits = run_search_in_memory(&queries, &store, &config).unwrap();

        let out = tempfile::NamedTempFile::with_suffix(".tsv").unwrap();
        run_search(&queries, &store, &config, out.path()).unwrap();
        let file_hit_count = std::fs::read_to_string(out.path())
            .unwrap()
            .lines()
            .filter(|line| !line.trim().is_empty())
            .count();

        assert_eq!(
            hits.len(),
            file_hit_count,
            "run_search_in_memory and file output must report the same hit count"
        );
    }

    #[test]
    fn in_memory_empty_for_non_matching_target() {
        let root = workspace_root();
        let query_path = root.join("tests/data/query.fa");

        let mut target_file = tempfile::NamedTempFile::new().unwrap();
        write!(target_file, ">dummy\nAAAAAAAAAAAAAAAA\n").unwrap();
        let (store, _tmp) = build_store(target_file.path());

        let config = SearchConfig {
            filter: FilterConfig {
                delta_g: Energy::from_kcal(-100.0),
                ..test_config().filter
            },
            ..test_config()
        };
        let queries = QueryRegistry::from_fasta(&query_path, &config.seed).unwrap();
        let hits = run_search_in_memory(&queries, &store, &config).unwrap();

        assert!(
            hits.is_empty(),
            "no hits expected against a non-matching target"
        );
    }

    #[test]
    fn from_fastas_same_result_as_from_fasta() {
        let root = workspace_root();
        let query_path = root.join("tests/data/query.fa");
        let target_path = root.join("tests/data/target.fa");

        let (store, _tmp) = build_store(&target_path);
        let config = test_config();

        let all = QueryRegistry::from_fasta(&query_path, &config.seed).unwrap();
        let single = QueryRegistry::from_fastas(&[query_path.as_path()], &config.seed).unwrap();
        assert_eq!(single.len(), all.len());

        let hits_single = run_search_in_memory(&single, &store, &config).unwrap();
        let hits_all = run_search_in_memory(&all, &store, &config).unwrap();
        assert_eq!(hits_single.len(), hits_all.len());
    }
}
