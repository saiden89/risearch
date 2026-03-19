//! Test runners for parity testing.
//!
//! Provides `RustRunner<S>` and `ParityRunner` for orchestrating
//! parity tests between Rust and C implementations.

use log::info;
use std::fs;
use std::path::{Path, PathBuf};

use crate::support::c_runner::{CRunner, NoIndex as CNoIndex};
use crate::support::comparison::{LogTag, ParityComparator};
use crate::support::status::TEST_PARITY_MODE;
use crate::support::{init_test_logging, parse_output, workspace_root};

// =============================================================================
// TYPE-STATE MARKERS
// =============================================================================

/// Marker for an unindexed Rust runner.
struct NoIndex;

/// Marker for an indexed Rust runner.
struct Indexed {
    index_path: PathBuf,
    target_store: risearch::TargetStore,
}

// =============================================================================
// RUST RUNNER
// =============================================================================

/// Runner for the Rust risearch implementation.
///
/// Uses type-state pattern to ensure you can only search after indexing.
struct RustRunner<S> {
    target_path: PathBuf,
    state: S,
}

impl RustRunner<NoIndex> {
    /// Create a new Rust runner for the given target.
    fn new(target: &Path) -> Self {
        Self {
            target_path: target.to_path_buf(),
            state: NoIndex,
        }
    }

    /// Create an index from the target file.
    /// Consumes self and returns an indexed runner.
    /// If `index_path` is None, uses a temporary directory.
    fn create_index(self, index_path: Option<&Path>) -> RustRunner<Indexed> {
        let (index_path, _tmpdir) = match index_path {
            Some(p) => (p.to_path_buf(), None),
            None => {
                let tmpdir = tempfile::tempdir().expect("tempdir");
                let path = tmpdir.path().join("target.idx");
                // Leak tempdir to keep index alive
                (path, Some(std::mem::ManuallyDrop::new(tmpdir)))
            }
        };

        risearch::TargetStore::build_from_fasta(&self.target_path, &index_path)
            .expect("build store index");
        let target_store = risearch::TargetStore::open(&index_path).expect("open target store");

        RustRunner {
            target_path: self.target_path,
            state: Indexed {
                index_path,
                target_store,
            },
        }
    }
}

impl RustRunner<Indexed> {
    /// Search using risearch as a library.
    fn search(
        &self,
        query_path: &Path,
        args: &risearch::config::SearchConfig,
    ) -> (Vec<risearch::SearchHit>, risearch::QueryRegistry) {
        let query_registry =
            risearch::QueryRegistry::from_fasta(query_path, &args.seed).expect("read query FASTA");
        let mut search_args = args.clone();
        // Parity parser expects binding-site columns (pairing + target sequence, optional flanks).
        search_args.output.format = risearch::config::OutputFormat::BindingSite;

        let tmp = tempfile::NamedTempFile::with_suffix(".tsv").unwrap();
        risearch::search::run_search(
            &query_registry,
            &self.state.target_store,
            &search_args,
            tmp.path(),
        )
        .expect("search");
        let rust_out = fs::read_to_string(tmp.path()).expect("read output");
        let (hits, _) = parse_output(&rust_out, &query_registry, &self.state.target_store);
        (hits, query_registry)
    }

    /// Get the index path.
    #[allow(dead_code)] // Useful API for debugging
    fn index_path(&self) -> &Path {
        &self.state.index_path
    }
}

// =============================================================================
// PARITY RUNNER
// =============================================================================

/// Combined runner for parity testing between Rust and C.
///
/// Manages both runners and provides comparison methods.
pub(crate) struct ParityRunner {
    rust: RustRunner<Indexed>,
    c: CRunner<crate::support::c_runner::Indexed>,
    #[allow(dead_code)]
    tmpdir: tempfile::TempDir,
}

impl ParityRunner {
    /// Create a new parity runner for the given target.
    /// Creates indexes for both Rust and C implementations.
    pub(crate) fn new(target: &Path) -> Self {
        init_test_logging();

        let root = workspace_root();
        let tmpdir = tempfile::tempdir().expect("tempdir");

        // Create Rust index
        let rust_index_path = tmpdir.path().join("target.idx");
        let rust = RustRunner::<NoIndex>::new(target).create_index(Some(&rust_index_path));

        // Create C index
        let c_index_path = tmpdir.path().join("target.pksuf");
        let c = CRunner::<CNoIndex>::new(&root).create_index(target, &c_index_path);

        Self { rust, c, tmpdir }
    }

    /// Compare Rust and C results, returning the parsed outputs.
    #[allow(dead_code)] // Useful API for detailed analysis
    fn compare(
        &self,
        query: &Path,
        args: &[&str],
    ) -> (
        Vec<risearch::SearchHit>,
        Vec<risearch::SearchHit>,
        risearch::QueryRegistry,
    ) {
        let search_args = parse_search_args(args);
        let (rust_hits, query_registry) = self.rust.search(query, &search_args);

        // Translate args to C format
        let c_args = translate_args_for_c(args);
        let c_args_ref: Vec<&str> = c_args.iter().map(|s| s.as_str()).collect();
        let c_out = self.c.search(query, &c_args_ref);
        let (c_hits, _) = parse_output(&c_out, &query_registry, &self.rust.state.target_store);

        (rust_hits, c_hits, query_registry)
    }

    /// Run comparison and assert parity passes.
    pub(crate) fn assert_pass(&self, query: &Path, test_name: &str, args: &[&str]) {
        let (rust_recs, c_recs, query_registry) = self.compare(query, args);

        let result = ParityComparator::new(&rust_recs, &c_recs).compare();

        info!(
            "{} {} - Rust={} hits, C={} hits, exact={}, extras={}, missings={}",
            LogTag::Parity,
            test_name,
            rust_recs.len(),
            c_recs.len(),
            result.exact_matches,
            result.extras.len(),
            result.missings.len()
        );

        // Detailed debug logging (all hit types with tables)
        result.log_details_with_context(
            test_name,
            Some(&query_registry),
            Some(&self.rust.state.target_store),
        );

        if !result.is_pass(*TEST_PARITY_MODE) {
            panic!(
                "\n{} FAILED [mode={:?}]: {} ({} co-optimal, {} rust-better, {} rust-worse, {} missing, {} extra)\nRun with RUST_LOG=debug for detailed diff analysis.\n",
                LogTag::Parity,
                *TEST_PARITY_MODE,
                test_name,
                result.co_optimal.len(),
                result.rust_better.len(),
                result.rust_worse.len(),
                result.missings.len(),
                result.extras.len()
            );
        }
    }
}

// =============================================================================
// SINGLE SEQ PARITY RUNNER
// =============================================================================

/// Runner for single-sequence parity tests.
///
/// Creates temporary FASTA files from raw sequences.
#[allow(dead_code)] // Used by some parity suites, unused in others.
pub(crate) struct SingleSeqRunner {
    query_path: PathBuf,
    #[allow(dead_code)] // Kept for potential future use
    target_path: PathBuf,
    runner: ParityRunner,
    #[allow(dead_code)]
    tmpdir: tempfile::TempDir,
}

#[allow(dead_code)] // Methods are only used by specific parity suites.
impl SingleSeqRunner {
    /// Create a new single-sequence runner.
    pub(crate) fn new(query_seq: &str, target_seq: &str) -> Self {
        let tmpdir = tempfile::tempdir().expect("tempdir");

        let query_path = tmpdir.path().join("query.fa");
        let target_path = tmpdir.path().join("target.fa");

        let q_upper = query_seq.to_uppercase();
        let t_upper = target_seq.to_uppercase();

        fs::write(&query_path, format!(">query\n{}\n", q_upper)).expect("write query");
        fs::write(&target_path, format!(">target\n{}\n", t_upper)).expect("write target");

        let runner = ParityRunner::new(&target_path);

        Self {
            query_path,
            target_path,
            runner,
            tmpdir,
        }
    }

    /// Run comparison and assert parity passes.
    pub(crate) fn assert_pass(&self, test_name: &str, args: &[&str]) {
        self.runner.assert_pass(&self.query_path, test_name, args);
    }
}

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

/// Translate Rust CLI args to C-compatible args.
///
/// Handles:
/// - Rust-only flags (--no-max-prune, --experimental, --dp-band*)
/// - New seed syntax (--seed-start/end/length) → legacy -s format
/// - Flag name differences (--no-seed-wobble → --noGUseed)
fn translate_args_for_c(args: &[&str]) -> Vec<String> {
    let mut c_args: Vec<String> = Vec::new();
    let mut seed_start: Option<&str> = None;
    let mut seed_end: Option<&str> = None;
    let mut seed_length: Option<&str> = None;
    let mut has_legacy_seed = false;

    let mut iter = args.iter().copied().peekable();
    while let Some(arg) = iter.next() {
        match arg {
            // Rust-only flags - skip
            "--no-max-prune" | "--experimental" => continue,
            "--dp-band" | "--dp-band-mode" => {
                let _ = iter.next(); // skip value
                continue;
            }

            // Legacy seed format - pass through
            "-s" | "--seed" => {
                has_legacy_seed = true;
                if let Some(val) = iter.next() {
                    c_args.push("-s".to_string());
                    c_args.push(val.to_string());
                }
            }
            _ if arg.starts_with("--seed=") => {
                has_legacy_seed = true;
                c_args.push("-s".to_string());
                c_args.push(arg["--seed=".len()..].to_string());
            }

            // New seed syntax - collect for later
            "--seed-start" => seed_start = iter.next(),
            "--seed-end" => seed_end = iter.next(),
            "--seed-length" => seed_length = iter.next(),
            _ if arg.starts_with("--seed-start=") => {
                seed_start = Some(&arg["--seed-start=".len()..])
            }
            _ if arg.starts_with("--seed-end=") => seed_end = Some(&arg["--seed-end=".len()..]),
            _ if arg.starts_with("--seed-length=") => {
                seed_length = Some(&arg["--seed-length=".len()..])
            }

            // Pairing flags - translate names
            "--no-seed-wobble" | "--noGUseed" => {
                c_args.push("--noGUseed".to_string())
            }
            // Pass through everything else
            _ => c_args.push(arg.to_string()),
        }
    }

    // Convert new seed syntax to legacy format
    if !has_legacy_seed {
        let spec = match (seed_start, seed_end, seed_length) {
            (Some(s), Some(e), Some(l)) => Some(format!("{}:{}/{}", s, e, l)),
            (Some(s), Some(e), None) => Some(format!("{}:{}", s, e)),
            (None, None, Some(l)) => Some(l.to_string()),
            _ => None,
        };
        if let Some(s) = spec {
            c_args.extend(["-s".to_string(), s]);
        }
    }

    c_args
}

/// Parse CLI-style args into SearchArgs using clap.
fn parse_search_args(args: &[&str]) -> risearch::config::SearchConfig {
    use clap::Parser;

    let mut cli_args: Vec<String> = vec!["risearch".into()];
    let mut iter = args.iter().copied().peekable();
    while let Some(arg) = iter.next() {
        let mapped = match arg {
            "-p" => {
                if let Some(next) = iter.peek().copied() {
                    let fmt = match next {
                        "1" => Some("detailed"),
                        "2" => Some("cigar"),
                        "3" => Some("bindingsite"),
                        "4" => Some("minimal"),
                        _ => None,
                    };
                    if let Some(fmt) = fmt {
                        let _ = iter.next();
                        Some(format!("-f={}", fmt))
                    } else {
                        Some("-f=detailed".to_string())
                    }
                } else {
                    Some("-f=detailed".to_string())
                }
            }
            "-p1" => Some("-f=detailed".to_string()),
            "-p2" => Some("-f=cigar".to_string()),
            "-p3" => Some("-f=bindingsite".to_string()),
            "-p4" => Some("-f=minimal".to_string()),
            _ => None,
        };

        if let Some(m) = mapped {
            cli_args.push(m);
        } else {
            cli_args.push(arg.to_string());
        }
    }
    #[derive(Parser)]
    struct FakeCmd {
        #[command(flatten)]
        search: risearch::cli::args::SearchArgs,
    }

    let parsed = FakeCmd::try_parse_from(&cli_args).expect("Failed to parse search args");
    parsed.search.into()
}
