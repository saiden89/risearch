//! Search module - finds miRNA-target interactions.
//!
//! Pipeline (parallel over queries; each query is processed end-to-end by one
//! worker):
//! 1) Build the query's own suffix array, traverse it against the shared target
//!    SA to enumerate seeds
//! 2) Extend each seed (DP) into a final hit
//! 3) Collapse hits that share a final bounding box to the lowest-energy one
//!    (unless `--no-dedup`), then format and emit

mod extension;

use anyhow::{bail, Context, Result};
use log::info;
use rayon::prelude::*;
use smallvec::SmallVec;
use std::collections::HashMap;
use std::ops::Range;
use std::path::Path;
use std::sync::Mutex;

use self::extension::ExtensionEngine;
use crate::alignment::{Alignment, PairClass};
use crate::config::SearchConfig;
use crate::dp::{DpConfig, DpView, ExtendDir, MAX_EXT};
use crate::dsm::{DsmRegistry, ScoringModel};
use crate::index::store::TargetRegistry;
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

/// Run search collecting all hits into memory.
pub fn run_search_in_memory(
    queries: &QueryRegistry,
    store: &TargetRegistry,
    opts: &SearchConfig,
) -> Result<Vec<SearchHit>> {
    let Some(ctx) = init_search(queries, store, opts)? else {
        return Ok(Vec::new());
    };
    let seeds = SeedingEngine::new(ctx.queries, ctx.store).run(&ctx.opts.seed)?;

    let dedup = !ctx.opts.filter.no_dedup;
    let include_alignment = ctx.opts.extend.build_alignment;
    let hits: Vec<SearchHit> = seeds
        .into_par_iter()
        .flat_map_iter(|(qi, seeds)| {
            let mut worker = SearchWorker::new(ctx.opts, &ctx.model);
            let query = &ctx.queries.entries()[qi];
            let query_seq = query.sequence();
            let seed_interval = query.seed_interval.clone();
            let produced: Vec<SearchHit> = seeds
                .into_iter()
                .filter_map(|seed| {
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
                        include_alignment,
                        target_len,
                        &seed,
                        target_trans,
                    )
                })
                .collect();
            // Dedup is exact per query: every box-mate of `qi` is in `produced`.
            let group = if dedup {
                dedup_hits(produced)
            } else {
                produced
            };
            group.into_iter()
        })
        .collect();

    Ok(hits)
}

/// Key identifying one final bounding box: a `(query, target, strand)` plus the
/// extended span. Overlapping seeds that converge on the same box collapse to
/// one row under [`hit_is_better`].
#[derive(PartialEq, Eq, Hash)]
struct BoxKey {
    query_idx: usize,
    target_idx: usize,
    strand: Strand,
    q_start: usize,
    q_end: usize,
    t_start: usize,
    t_end: usize,
}

impl BoxKey {
    fn of(hit: &SearchHit) -> Self {
        Self {
            query_idx: hit.query_idx,
            target_idx: hit.target_idx,
            strand: hit.strand,
            q_start: hit.q_start,
            q_end: hit.q_end,
            t_start: hit.t_start,
            t_end: hit.t_end,
        }
    }
}

/// Deterministic total order: a candidate beats the incumbent iff it has a more
/// negative energy, or — on an exact energy tie — a lexicographically smaller
/// pairing fingerprint. Both terms are pure functions of a hit's own fields, so
/// the surviving set is independent of insertion/scheduling order. (The tie-break
/// only distinguishes alignment-printing formats; Minimal rows that tie on energy
/// are output-identical regardless of which wins.)
fn hit_is_better(candidate: &SearchHit, current: &SearchHit) -> bool {
    use std::cmp::Ordering;
    match candidate.energy.cmp(&current.energy) {
        Ordering::Less => true,
        Ordering::Greater => false,
        Ordering::Equal => fingerprint_cmp(candidate, current) == Ordering::Less,
    }
}

/// Lexicographic compare of the two hits' pairing fingerprints without
/// allocating. Note `PairClass`'s own variant order differs from `symbol()`
/// order, so compare the mapped symbol chars. A missing alignment is an empty
/// fingerprint, which sorts before any non-empty one.
fn fingerprint_cmp(a: &SearchHit, b: &SearchHit) -> std::cmp::Ordering {
    fn syms(hit: &SearchHit) -> impl Iterator<Item = char> + '_ {
        hit.alignment
            .iter()
            .flat_map(|a| a.steps().iter().map(|p| p.symbol()))
    }
    syms(a).cmp(syms(b))
}

/// Collapse hits sharing a [`BoxKey`] to the single best one. Caller guarantees
/// all hits belong to one query, so this is globally exact for that query.
fn dedup_hits(hits: Vec<SearchHit>) -> Vec<SearchHit> {
    use std::collections::hash_map::Entry;
    let mut best: HashMap<BoxKey, SearchHit> = HashMap::with_capacity(hits.len());
    for hit in hits {
        match best.entry(BoxKey::of(&hit)) {
            Entry::Occupied(mut e) => {
                if hit_is_better(&hit, e.get()) {
                    e.insert(hit);
                }
            }
            Entry::Vacant(e) => {
                e.insert(hit);
            }
        }
    }
    best.into_values().collect()
}

/// Run search and write hits to `output_path` (or directory in multifile mode).
///
/// Parallelizes over queries: each rayon worker owns one query end-to-end —
/// seed it against the shared target SA, extend, format, and write its block —
/// then frees its per-query state. Peak RAM is bounded by the concurrent
/// workers' per-query state, not the total hit count. Output line order is not
/// stable across runs (the hit set is); callers needing order must sort.
pub fn run_search(
    queries: &QueryRegistry,
    store: &TargetRegistry,
    opts: &SearchConfig,
    output_path: &Path,
) -> Result<()> {
    let Some(ctx) = init_search(queries, store, opts)? else {
        return Ok(());
    };
    let engine = SeedingEngine::new(ctx.queries, ctx.store);

    let init = || {
        (
            SearchWorker::new(ctx.opts, &ctx.model),
            HitFormatter::new(ctx.opts.output.format),
        )
    };

    let counts: Vec<usize> = if opts.output.multifile {
        fs_err::create_dir_all(output_path).with_context(|| {
            format!(
                "Failed to create output directory {:?} for --multifile",
                output_path
            )
        })?;
        let ext = opts.output.compress.extension();
        let paths = build_multifile_paths(queries, output_path, ext);

        (0..ctx.queries.len())
            .into_par_iter()
            .map_init(init, |(worker, fmt), qi| -> Result<usize> {
                let (emitted, chunks) = process_one_query(&ctx, &engine, worker, qi, fmt)?;
                // One writer per query (matches RIsearch2); no file for an
                // empty query.
                if !chunks.is_empty() {
                    let mut writer = OutputWriter::new(&ctx.opts.output, &paths[qi])?;
                    write_chunks(&mut writer, &chunks)?;
                    writer.flush_all()?;
                }
                Ok(emitted)
            })
            .collect::<Result<Vec<usize>>>()?
    } else {
        let writer = Mutex::new(OutputWriter::new(&opts.output, output_path)?);

        let counts = (0..ctx.queries.len())
            .into_par_iter()
            .map_init(init, |(worker, fmt), qi| -> Result<usize> {
                let (emitted, chunks) = process_one_query(&ctx, &engine, worker, qi, fmt)?;
                // Format off-lock; take the writer once to append this query's
                // block (one bulk write per query).
                write_chunks(&mut writer.lock().unwrap(), &chunks)?;
                Ok(emitted)
            })
            .collect::<Result<Vec<usize>>>()?;

        writer.into_inner().unwrap().flush_all()?;
        counts
    };

    let total: usize = counts.iter().sum();
    info!("Search complete: {} hits", total);
    Ok(())
}

/// Seed one query against the shared target SA, then extend + format its hits.
/// All of a query's seeds are produced and consumed here, so the worker owns the
/// query end-to-end.
fn process_one_query(
    ctx: &SearchContext<'_>,
    engine: &SeedingEngine<'_>,
    worker: &mut SearchWorker,
    qi: usize,
    fmt: &mut HitFormatter,
) -> Result<(usize, Vec<OutputChunk>)> {
    let mut seeds: Vec<SeedHit> = Vec::new();
    engine.seed_query(qi, &ctx.opts.seed, |seed| seeds.push(seed))?;
    worker.process_query_seeds(ctx, qi, &seeds, fmt)
}

struct SearchContext<'a> {
    queries: &'a QueryRegistry,
    store: &'a TargetRegistry,
    opts: &'a SearchConfig,
    model: ScoringModel,
}

impl<'a> SearchContext<'a> {
    fn new(
        queries: &'a QueryRegistry,
        store: &'a TargetRegistry,
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
    store: &'a TargetRegistry,
    opts: &'a SearchConfig,
) -> Result<Option<SearchContext<'a>>> {
    if store.is_empty() || queries.is_empty() {
        return Ok(None);
    }

    // Unlimited extension (`-l -1`) promises to span the whole query, but the DP
    // buffers cap each side at MAX_EXT. Refuse rather than silently clamp: a
    // query longer than the cap cannot be served as requested. (Mirrors clap
    // rejecting an explicit `-l > MAX_EXT`.)
    if opts.extend.max_extension < 0 {
        for (_, q) in queries.iter() {
            let n = q.sequence().len();
            if n > MAX_EXT {
                bail!(
                    "query '{}' is {} nt; `-l -1` cannot extend across it ({} nt cap). \
                     Pass an explicit `-l <={}` to accept the cap, or shorten the query.",
                    q.name(),
                    n,
                    MAX_EXT,
                    MAX_EXT
                );
            }
        }
    }

    let max_ext = if opts.extend.max_extension < 0 {
        format!("unlimited(<={MAX_EXT})")
    } else {
        opts.extend.max_extension.to_string()
    };
    info!(
        "Starting search: {} queries x {} targets, seed_length={:?}, max_ext={}, delta_g={}",
        queries.len(),
        store.len(),
        opts.seed.seed_length,
        max_ext,
        opts.filter.delta_g
    );

    let (initiation, source_table) = DsmRegistry::load(&opts.score.dsm_id, opts.score.temperature)?;

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
        let dedup = !ctx.opts.filter.no_dedup;
        let include_alignment = ctx.opts.extend.build_alignment;

        let mut hits = Vec::new();
        for seed in seeds {
            let (t_fwd, t_rc, target_len) = ctx.store.target_slices(seed.target_idx);
            let target_trans = match seed.strand {
                Strand::Forward => t_fwd,
                Strand::Reverse => t_rc,
            };
            if let Some(hit) = self.build_hit_from_seed(
                ctx.opts,
                query_idx,
                query_seq,
                seed_interval.clone(),
                include_alignment,
                target_len,
                seed,
                target_trans,
            ) {
                hits.push(hit);
            }
        }
        if dedup {
            hits = dedup_hits(hits);
        }

        let mut chunks = Vec::new();
        for hit in &hits {
            let (t_fwd, t_rc, _) = ctx.store.target_slices(hit.target_idx);
            let t_name = ctx.store.get_name(hit.target_idx);
            if let Some(chunk) = format.add_hit(hit, query_name, query_seq, t_name, t_fwd, t_rc) {
                chunks.push(chunk);
            }
        }
        if let Some(chunk) = format.flush() {
            chunks.push(chunk);
        }
        Ok((hits.len(), chunks))
    }

    // Cohesive per-seed extension inputs (query/target buffers, seed interval,
    // and output flags); bundling them would only add indirection.
    #[allow(clippy::too_many_arguments)]
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

        let duplex_score =
            self.model
                .ungapped_duplex_score(query_bases, target_trans, q_start, t_start, len);
        let dp_cfg = DpConfig::from(&opts.extend);
        let ext_left = dp_cfg.side_cap(q_start + 1);
        let ext_right = dp_cfg.side_cap(query_bases.len() - (q_start + len - 1));
        let view_left = DpView::new(
            query_bases,
            target_trans,
            q_start,
            t_match_end,
            ExtendDir::Left,
            ext_left,
        );
        let left = self.extension.extend(&view_left, include_alignment);

        let view_right = DpView::new(
            query_bases,
            target_trans,
            q_start + len - 1,
            t_start,
            ExtendDir::Right,
            ext_right,
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

    // Cohesive hit-assembly inputs (query/target buffers, seed, both extension
    // results, and energy); bundling them would only add indirection.
    #[allow(clippy::too_many_arguments)]
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
                Some(Alignment::from_parts(left_pairs, &seed_pairs, right_pairs)),
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

    use super::*;
    use crate::config::{
        ExtendConfig, FilterConfig, OutputCompression, OutputConfig, OutputFormat, ScoreConfig,
        SeedConfig,
    };
    use crate::index::store::TargetRegistry;
    use crate::registry::QueryRegistry;
    use crate::types::DsmId;

    const QUERY_FA: &str = include_str!("../../tests/data/query.fa");
    const TARGET_FA: &str = include_str!("../../tests/data/target.fa");

    fn fixture(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
        f
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
                build_alignment: true,
            },
            filter: FilterConfig {
                delta_g: Energy::from_kcal(-10.0),
                seed_energy: Energy::from_kcal(0.0),
                no_max_prune: false,
                no_dedup: false,
            },
            output: OutputConfig {
                format: OutputFormat::Detailed,
                compress: OutputCompression::None,
                multifile: false,
            },
        }
    }

    fn build_store(target_fa: &std::path::Path) -> (TargetRegistry, tempfile::TempDir) {
        let tmpdir = tempfile::tempdir().unwrap();
        let idx = tmpdir.path().join("target.idx");
        TargetRegistry::build(target_fa, &idx, None).unwrap();
        let store = TargetRegistry::open(&idx).unwrap();
        (store, tmpdir)
    }

    #[test]
    fn in_memory_matches_file_hit_count() {
        let query_f = fixture(QUERY_FA);
        let target_f = fixture(TARGET_FA);

        let (store, _tmp) = build_store(target_f.path());
        let mut config = test_config();
        config.output.format = OutputFormat::Minimal;
        // Mirror what the CLI derives for Minimal, so this also covers dedup's
        // empty-fingerprint (first-wins) tie-break.
        config.extend.build_alignment = false;
        let queries = QueryRegistry::from_fasta(query_f.path(), &config.seed).unwrap();

        let hits = run_search_in_memory(&queries, &store, &config).unwrap();

        let out = tempfile::NamedTempFile::with_suffix(".tsv").unwrap();
        run_search(&queries, &store, &config, out.path()).unwrap();
        let file_hit_count = fs_err::read_to_string(out.path())
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
        let query_f = fixture(QUERY_FA);

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
        let queries = QueryRegistry::from_fasta(query_f.path(), &config.seed).unwrap();
        let hits = run_search_in_memory(&queries, &store, &config).unwrap();

        assert!(
            hits.is_empty(),
            "no hits expected against a non-matching target"
        );
    }

    fn run_to_lines(
        config: &SearchConfig,
        queries: &QueryRegistry,
        store: &TargetRegistry,
    ) -> Vec<String> {
        let out = tempfile::NamedTempFile::with_suffix(".tsv").unwrap();
        run_search(queries, store, config, out.path()).unwrap();
        fs_err::read_to_string(out.path())
            .unwrap()
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(str::to_string)
            .collect()
    }

    /// Split a Minimal-format row into (box key, energy). Energy is always the
    /// last tab field; everything before it identifies the bounding box.
    fn split_box_energy(line: &str) -> (&str, f64) {
        let (key, energy) = line.rsplit_once('\t').unwrap();
        (key, energy.parse().unwrap())
    }

    fn mm2_minimal_fixture() -> (
        TargetRegistry,
        tempfile::TempDir,
        QueryRegistry,
        SearchConfig,
    ) {
        let target_f = fixture(TARGET_FA);
        let (store, tmp) = build_store(target_f.path());
        let mut config = test_config();
        config.output.format = OutputFormat::Minimal;
        config.extend.build_alignment = false;
        // Parameters that produce overlapping-seed box collisions on this
        // fixture (~33% redundant rows): short seed, mismatches, wide window,
        // permissive energy.
        config.seed.seed_length = Some(6);
        config.seed.max_mismatches = 2;
        config.extend.max_extension = 30;
        config.filter.delta_g = Energy::from_kcal(-5.0);
        let query_f = fixture(QUERY_FA);
        let queries = QueryRegistry::from_fasta(query_f.path(), &config.seed).unwrap();
        (store, tmp, queries, config)
    }

    /// Default dedup keeps exactly one row per bounding box, and that row carries
    /// the minimum energy among the box's collapsed candidates. `--no-dedup`
    /// reproduces the full per-seed row set.
    #[test]
    fn dedup_keeps_one_min_energy_row_per_box() {
        use std::collections::{HashMap, HashSet};
        let (store, _tmp, queries, config) = mm2_minimal_fixture();

        let mut nodedup = config.clone();
        nodedup.filter.no_dedup = true;
        let raw = run_to_lines(&nodedup, &queries, &store);
        let dedup = run_to_lines(&config, &queries, &store);

        assert!(!dedup.is_empty(), "fixture must yield hits");
        assert!(
            dedup.len() < raw.len(),
            "dedup must remove overlapping rows"
        );

        let mut min_energy: HashMap<String, f64> = HashMap::new();
        for l in &raw {
            let (k, e) = split_box_energy(l);
            min_energy
                .entry(k.to_string())
                .and_modify(|m| {
                    if e < *m {
                        *m = e;
                    }
                })
                .or_insert(e);
        }

        let mut seen = HashSet::new();
        for l in &dedup {
            let (k, e) = split_box_energy(l);
            assert!(
                seen.insert(k.to_string()),
                "duplicate box survived dedup: {k}"
            );
            assert!(
                (e - min_energy[k]).abs() < 1e-9,
                "kept energy {e} is not the box minimum {} for {k}",
                min_energy[k]
            );
        }
        assert_eq!(
            seen.len(),
            min_energy.len(),
            "dedup must keep exactly one row per distinct box"
        );
    }

    /// The deduped set must be stable across runs. `dedup_hits` drains a HashMap,
    /// so a tie-break that depended on hash/iteration order would surface here as
    /// run-to-run drift (two HashMaps use different random seeds). Locks the
    /// order-independence parity tests can't (parity is pinned to `--no-dedup`).
    #[test]
    fn dedup_result_is_deterministic() {
        let (store, _tmp, queries, config) = mm2_minimal_fixture();
        let mut a = run_to_lines(&config, &queries, &store);
        let mut b = run_to_lines(&config, &queries, &store);
        a.sort();
        b.sort();
        assert!(!a.is_empty(), "fixture must yield hits");
        assert_eq!(a, b, "deduped set must not depend on hash/iteration order");
    }
}
