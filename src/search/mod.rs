//! Search module - finds miRNA-target interactions.
//!
//! Pipeline (parallel over queries; each query is processed end-to-end by one
//! worker):
//! 1) Build the query's own suffix array, traverse it against the shared target
//!    SA to enumerate seeds
//! 2) Extend each seed (DP) into a final hit
//! 3) Collapse hits that share a final bounding box to the lowest-energy one
//!    (unless `--no-dedup`), then hand the query's hits to a [`HitSink`]

mod extension;

use anyhow::{bail, Result};
use log::info;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::ops::Range;
use std::sync::Mutex;

use self::extension::{ExtensionEngine, SeedExtension};
use crate::alignment::Alignment;
use crate::config::{SearchConfig, UNLIMITED_EXTENSION};
use crate::dp::MAX_EXT;
use crate::dsm::ScoringModel;
use crate::index::store::TargetRegistry;
use crate::index::view::TargetView;
use crate::registry::QueryRegistry;
use crate::seed::{SeedHit, SeedingEngine};
use crate::types::{Base, Energy, Strand};

/// One accepted interaction: a span of a query paired against a span of a target.
///
/// Coordinates are 0-based and inclusive at both ends; the output formats add 1.
/// `t_start`/`t_end` always address the original FASTA sequence. Internally the
/// search uses physical duplex views ordered 3'->5' alongside the 5'->3' query;
/// that span is converted to FASTA coordinates once when the hit is assembled.
///
/// Bases and names are not carried, only indices into the query and target
/// registries. [`query`](Self::query) slices a supplied query sequence;
/// [`target`](Self::target) resolves the paired span through a [`TargetView`].
///
/// `alignment` is populated only under [`ExtendConfig::build_alignment`](crate::ExtendConfig),
/// and carries the seed's span within its own steps.
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
    // Boxed: the hit array is grown and moved per seed, and an inline `Alignment`
    // would make every element carry its column buffer.
    pub alignment: Option<Box<Alignment>>,
}

/// Where a search's hits go.
///
/// The driver runs the seed -> extend -> dedup pipeline and hands each query's
/// surviving hits to one of these; the destination (one text file, one file per
/// query, Arrow columns) is entirely the implementor's business.
///
/// `consume` is called at most once per `query_idx`, from a rayon worker, in
/// arbitrary order — hence `&self` plus interior mutability. Lock once per call,
/// never per hit, and do any expensive conversion before taking the lock.
pub trait HitSink: Sync {
    fn consume(&self, query_idx: usize, hits: Vec<SearchHit>) -> Result<()>;

    /// Finalize the destination. The driver calls this once every query has been
    /// consumed and the search succeeded — including when no hits arrived at all,
    /// which is what lets a zero-hit run truncate a stale output file. A run with
    /// nothing to search (empty target or query set) returns without consuming or
    /// flushing anything.
    fn flush(&self) -> Result<()> {
        Ok(())
    }
}

/// Collects every hit into one `Vec`. Peak RAM scales with the hit count, so
/// prefer a sink that writes through for large searches.
#[derive(Default)]
pub struct VecSink(Mutex<Vec<SearchHit>>);

impl VecSink {
    pub fn into_hits(self) -> Vec<SearchHit> {
        self.0.into_inner().unwrap()
    }
}

impl HitSink for VecSink {
    fn consume(&self, _query_idx: usize, hits: Vec<SearchHit>) -> Result<()> {
        self.0.lock().unwrap().extend(hits);
        Ok(())
    }
}

/// Key identifying one final bounding box: a `(query, target, strand)` plus the
/// extended span. Overlapping seeds that converge on the same box collapse to
/// one row under [`SearchHit::cmp`].
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

/// Collapse hits sharing a [`BoxKey`] to the single best one. Caller guarantees
/// all hits belong to one query, so this is globally exact for that query.
fn dedup_hits(hits: Vec<SearchHit>) -> Vec<SearchHit> {
    use std::collections::hash_map::Entry;
    let mut best: HashMap<BoxKey, SearchHit> = HashMap::with_capacity(hits.len());
    for hit in hits {
        match best.entry(BoxKey::of(&hit)) {
            Entry::Occupied(mut e) => {
                if hit.cmp(e.get()).is_lt() {
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

/// Run search, handing every query's hits to `sink` as they are produced.
/// Returns the number of hits handed over.
///
/// Parallelizes over queries: each rayon worker owns one query end-to-end —
/// seed it against the shared target SA, extend, dedup — then hands the hits to
/// the sink and frees its per-query state. Peak RAM is therefore bounded by the
/// concurrent workers' per-query state plus whatever the sink retains, not by the
/// total hit count.
///
/// Queries arrive at the sink in completion order, and hits within a query in
/// `HashMap` order, so nothing about the output order is stable across runs. The
/// hit *set* is. Callers needing an order must sort.
///
/// The config is validated before any query runs, so a sink that defers touching
/// its destination until first use will not have disturbed it if this returns an
/// error.
pub fn run_search(
    queries: &QueryRegistry,
    store: &TargetRegistry,
    opts: &SearchConfig,
    sink: &dyn HitSink,
) -> Result<usize> {
    opts.validate().map_err(anyhow::Error::msg)?;

    if store.is_empty() || queries.is_empty() {
        // Deliberately not flushed: an empty store or query set leaves the
        // destination untouched, where a zero-hit run truncates it.
        return Ok(0);
    }

    check_unlimited_fits(queries, opts)?;
    log_search_banner(queries, store, opts);

    let ctx = SearchContext {
        queries,
        store,
        opts,
        model: ScoringModel::load(
            &opts.score.dsm_id,
            opts.score.temperature,
            opts.score.penalty,
        )?,
    };
    let engine = SeedingEngine::new(ctx.queries, ctx.store);
    let counts: Vec<usize> = (0..ctx.queries.len())
        .into_par_iter()
        .map_init(
            || SearchWorker::new(ctx.opts, &ctx.model),
            |worker, qi| -> Result<usize> {
                let hits = worker.search_query(&ctx, &engine, qi)?;
                let emitted = hits.len();
                sink.consume(qi, hits)?;
                Ok(emitted)
            },
        )
        .collect::<Result<Vec<usize>>>()?;
    sink.flush()?;

    let total = counts.iter().sum();
    info!("Search complete: {} hits", total);
    Ok(total)
}

struct SearchContext<'a> {
    queries: &'a QueryRegistry,
    store: &'a TargetRegistry,
    opts: &'a SearchConfig,
    model: ScoringModel,
}

/// A negative `-l` is the unlimited sentinel: extend across the whole query
/// rather than a fixed length.
fn is_unlimited(max_extension: i32) -> bool {
    max_extension == UNLIMITED_EXTENSION
}

/// Unlimited extension (`-l -1`) promises to span the whole query, but the DP
/// buffers cap each side at MAX_EXT. Refuse rather than silently clamp: a query
/// longer than the cap cannot be served as requested. (Mirrors clap rejecting an
/// explicit `-l > MAX_EXT`.)
fn check_unlimited_fits(queries: &QueryRegistry, opts: &SearchConfig) -> Result<()> {
    if !is_unlimited(opts.extend.max_extension) {
        return Ok(());
    }
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
    Ok(())
}

fn log_search_banner(queries: &QueryRegistry, store: &TargetRegistry, opts: &SearchConfig) {
    let max_ext = if is_unlimited(opts.extend.max_extension) {
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
}

struct SearchWorker {
    extension: ExtensionEngine,
}

impl SearchWorker {
    fn new(opts: &SearchConfig, model: &ScoringModel) -> Self {
        let max_window = (!is_unlimited(opts.extend.max_extension))
            .then_some(opts.extend.max_extension as usize);
        Self {
            extension: ExtensionEngine::from_parts(max_window, model),
        }
    }

    /// Seed one query, extend and score every seed, and collapse the survivors.
    fn search_query(
        &mut self,
        ctx: &SearchContext<'_>,
        engine: &SeedingEngine<'_>,
        query_idx: usize,
    ) -> Result<Vec<SearchHit>> {
        let opts = ctx.opts;
        let query = ctx.queries.get(query_idx).sequence();

        // Loop-invariant: rebuilding it per seed re-enters the rkyv root twice.
        let tview = ctx.store.view();
        let mut hits = Vec::new();
        engine.seed_query(query_idx, &opts.seed, |seed| {
            debug_assert_eq!(seed.query_idx(), query_idx);
            let target = tview.target(seed.target_idx(), seed.strand());
            let target_range = seed.target_range();
            debug_assert!(target_range.end <= target.len());

            let ext = self
                .extension
                .extend_seed(query, target, &seed, opts.extend.build_alignment);

            if ext.energy <= opts.filter.delta_g {
                hits.push(SearchHit::new(&seed, ext, tview));
            }
        })?;
        // Dedup is exact per query: every box-mate of `query_idx` is in `hits`.
        Ok(if opts.filter.no_dedup {
            hits
        } else {
            dedup_hits(hits)
        })
    }
}

/// Convert a non-empty internal half-open range to public inclusive bounds.
fn inclusive_bounds(mut range: Range<usize>) -> (usize, usize) {
    let start = range.start;
    let end = range
        .next_back()
        .expect("extended duplex ranges must be non-empty");
    (start, end)
}

impl SearchHit {
    /// Compare two candidates for the same [`BoxKey`].
    ///
    /// `Less` means `self` is retained, `Greater` means `other` is retained, and
    /// `Equal` means they are indistinguishable under the deduplication policy.
    /// Energy orders first; an exact tie is broken by the alignment fingerprint.
    fn cmp(&self, other: &Self) -> Ordering {
        self.energy.cmp(&other.energy).then_with(|| {
            match (self.alignment.as_deref(), other.alignment.as_deref()) {
                (Some(a), Some(b)) => a.fingerprint_cmp(b),
                (None, None) => Ordering::Equal,
                (None, Some(_)) => Ordering::Less,
                (Some(_), None) => Ordering::Greater,
            }
        })
    }

    pub fn query<'a>(&self, q_seq: &'a [Base]) -> &'a [Base] {
        &q_seq[self.q_start..=self.q_end]
    }

    pub fn target<'a>(&self, targets: TargetView<'a>) -> &'a [Base] {
        let target = targets.target(self.target_idx, self.strand);
        &target[self.duplex_target_range(targets)]
    }

    /// Map the public inclusive FASTA coordinates into a half-open range over
    /// the strand-selected duplex target.
    pub(crate) fn duplex_target_range(&self, targets: TargetView<'_>) -> Range<usize> {
        let fasta_range = self.t_start
            ..self
                .t_end
                .checked_add(1)
                .expect("inclusive target end must fit a half-open range");
        targets.map_target_range(self.target_idx, self.strand, fasta_range)
    }

    /// Convert the extension's duplex-frame target span to FASTA coordinates and
    /// record the hit.
    fn new(seed: &SeedHit, ext: SeedExtension, targets: TargetView<'_>) -> Self {
        let SeedExtension {
            q_range,
            t_range,
            energy: binding_energy,
            alignment,
        } = ext;
        let (q_start, q_end) = inclusive_bounds(q_range);
        let fasta_t_range = targets.map_target_range(seed.target_idx(), seed.strand(), t_range);
        let (t_start, t_end) = inclusive_bounds(fasta_t_range);

        Self {
            query_idx: seed.query_idx(),
            target_idx: seed.target_idx(),
            q_start,
            q_end,
            t_start,
            t_end,
            strand: seed.strand(),
            energy: binding_energy,
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
    use crate::output::TextSink;
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
                no_max_prune: false,
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
                no_dedup: false,
            },
        }
    }

    fn test_output() -> OutputConfig {
        OutputConfig {
            format: OutputFormat::Detailed,
            compress: OutputCompression::None,
            multifile: false,
        }
    }

    fn minimal_output() -> OutputConfig {
        OutputConfig {
            format: OutputFormat::Minimal,
            ..test_output()
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
    fn text_output_emits_one_line_per_hit() {
        let query_f = fixture(QUERY_FA);
        let target_f = fixture(TARGET_FA);

        let (store, _tmp) = build_store(target_f.path());
        let mut config = test_config();
        let output = minimal_output();
        // Mirror what the CLI derives for Minimal, so this also covers dedup's
        // empty-fingerprint (first-wins) tie-break.
        config.extend.build_alignment = false;
        let queries = QueryRegistry::from_fasta(query_f.path(), &config.seed).unwrap();

        let out = tempfile::NamedTempFile::with_suffix(".tsv").unwrap();
        let hits = run_to_path(&queries, &store, &config, &output, out.path());
        let file_hit_count = fs_err::read_to_string(out.path())
            .unwrap()
            .lines()
            .filter(|line| !line.trim().is_empty())
            .count();

        assert_eq!(
            hits, file_hit_count,
            "text output must emit exactly one line per retained hit"
        );
    }

    #[test]
    fn no_hits_against_non_matching_target() {
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
        let hits = run_search(&queries, &store, &config, &VecSink::default()).unwrap();

        assert_eq!(hits, 0, "no hits expected against a non-matching target");
    }

    #[test]
    fn run_search_validates_configs_constructed_without_clap() {
        let query_f = fixture(QUERY_FA);
        let target_f = fixture(TARGET_FA);
        let (store, _tmp) = build_store(target_f.path());
        let mut config = test_config();
        let queries = QueryRegistry::from_fasta(query_f.path(), &config.seed).unwrap();

        config.score.penalty = Energy::from_kcal(-0.1);
        let err = run_search(&queries, &store, &config, &VecSink::default()).unwrap_err();
        assert!(err.to_string().contains("penalty"));

        config.score.penalty = Energy::from_kcal(0.0);
        config.extend.max_extension = -2;
        let err = run_search(&queries, &store, &config, &VecSink::default()).unwrap_err();
        assert!(err.to_string().contains("max extension"));
    }

    /// Zero-hit behaviour differs by topology and neither branch is covered by a
    /// CLI test: multifile creates no file at all, single-file truncates the
    /// output to nothing.
    #[test]
    fn zero_hit_run_writes_empty_file_but_no_multifile_entry() {
        let query_f = fixture(QUERY_FA);
        let mut target_file = tempfile::NamedTempFile::new().unwrap();
        write!(target_file, ">dummy\nAAAAAAAAAAAAAAAA\n").unwrap();
        let (store, _tmp) = build_store(target_file.path());

        let mut config = test_config();
        config.filter.delta_g = Energy::from_kcal(-100.0);
        let mut output = test_output();
        let queries = QueryRegistry::from_fasta(query_f.path(), &config.seed).unwrap();

        let out = tempfile::NamedTempFile::with_suffix(".tsv").unwrap();
        fs_err::write(out.path(), b"stale").unwrap();
        run_to_path(&queries, &store, &config, &output, out.path());
        assert_eq!(fs_err::metadata(out.path()).unwrap().len(), 0);

        output.multifile = true;
        let tmp = tempfile::tempdir().unwrap();
        // Not pre-created, so this also pins that the directory itself is made.
        let dir = tmp.path().join("multi");
        run_to_path(&queries, &store, &config, &output, &dir);
        assert_eq!(fs_err::read_dir(&dir).unwrap().count(), 0);
    }

    /// The sink is built before the search validates the config, so it must not
    /// touch its destination until the driver hands it something.
    #[test]
    fn rejected_config_leaves_no_multifile_directory() {
        let query_f = fixture(&format!(">longq\n{}\n", "A".repeat(MAX_EXT + 1)));
        let target_f = fixture(TARGET_FA);
        let (store, _tmp) = build_store(target_f.path());

        let mut config = test_config();
        config.extend.max_extension = -1;
        let output = OutputConfig {
            multifile: true,
            ..test_output()
        };
        let queries = QueryRegistry::from_fasta(query_f.path(), &config.seed).unwrap();

        let tmp = tempfile::tempdir().unwrap();
        let dir = tmp.path().join("multi");
        let sink = TextSink::new(&queries, &store, &output, &dir).unwrap();
        assert!(run_search(&queries, &store, &config, &sink).is_err());
        assert!(!dir.exists());
    }

    /// Multifile addresses its writers by index (`paths[query_idx]`) while queries
    /// that produce nothing get no file at all, so a query dropping out in the
    /// middle must not shift the others' rows into the wrong file.
    #[test]
    fn multifile_writes_each_querys_rows_to_its_own_file() {
        let query_f = fixture(">q1\nAAAA\n>q2\nCCCC\n>q3\nAAAAA\n");
        let target_f = fixture(">t\nUUUUUUUU\n");
        let (store, _tmp) = build_store(target_f.path());

        let mut config = test_config();
        let output = OutputConfig {
            multifile: true,
            ..minimal_output()
        };
        config.extend.build_alignment = false;
        config.extend.max_extension = 0;
        config.seed.seed_length = Some(4);
        config.filter.delta_g = Energy::from_kcal(100.0);
        let queries = QueryRegistry::from_fasta(query_f.path(), &config.seed).unwrap();

        let dir = tempfile::tempdir().unwrap();
        run_to_path(&queries, &store, &config, &output, dir.path());

        let mut files: Vec<(String, String)> = fs_err::read_dir(dir.path())
            .unwrap()
            .map(|entry| {
                let path = entry.unwrap().path();
                let stem = path.file_stem().unwrap().to_string_lossy().into_owned();
                (stem, fs_err::read_to_string(&path).unwrap())
            })
            .collect();
        files.sort();

        let stems: Vec<&str> = files.iter().map(|(stem, _)| stem.as_str()).collect();
        assert_eq!(stems, ["q1", "q3"], "only queries with hits get a file");
        for (stem, body) in &files {
            assert!(!body.trim().is_empty(), "{stem} file is empty");
            for line in body.lines() {
                assert_eq!(
                    line.split('\t').next().unwrap(),
                    stem,
                    "row landed in the wrong query's file"
                );
            }
        }
    }

    /// The [`HitSink`] contract the in-tree sinks depend on: every query
    /// is consumed exactly once, and every hit in a batch carries the
    /// `query_idx` handed alongside it (multifile picks its file from that index,
    /// single-file looks up the query's name and sequence).
    #[test]
    fn streaming_consumes_each_query_once_with_matching_hit_indices() {
        #[derive(Default)]
        struct Spy(Mutex<Vec<usize>>);

        impl HitSink for Spy {
            fn consume(&self, query_idx: usize, hits: Vec<SearchHit>) -> Result<()> {
                assert!(hits.iter().all(|hit| hit.query_idx == query_idx));
                self.0.lock().unwrap().push(query_idx);
                Ok(())
            }
        }

        let (store, _tmp, queries, config) = mm2_minimal_fixture();
        let spy = Spy::default();
        run_search(&queries, &store, &config, &spy).unwrap();

        let mut seen = spy.0.into_inner().unwrap();
        seen.sort_unstable();
        assert_eq!(seen, (0..queries.len()).collect::<Vec<_>>());
    }

    /// Search into a text file at `path` (a directory under `--multifile`),
    /// returning the hit count. The sink dies here, so the caller reads a file
    /// whose compression trailer is already written.
    fn run_to_path(
        queries: &QueryRegistry,
        store: &TargetRegistry,
        config: &SearchConfig,
        output: &OutputConfig,
        path: &std::path::Path,
    ) -> usize {
        let sink = TextSink::new(queries, store, output, path).unwrap();
        run_search(queries, store, config, &sink).unwrap()
    }

    fn run_to_lines(
        config: &SearchConfig,
        queries: &QueryRegistry,
        store: &TargetRegistry,
    ) -> Vec<String> {
        let out = tempfile::NamedTempFile::with_suffix(".tsv").unwrap();
        run_to_path(queries, store, config, &minimal_output(), out.path());
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
