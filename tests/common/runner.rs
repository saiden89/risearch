//! Test runners for parity testing.
//!
//! Provides `RustRunner<S>` and `ParityRunner` for orchestrating
//! parity tests between Rust and C implementations.

use log::info;
use std::fs;
use std::path::{Path, PathBuf};

use crate::common::c_runner::{CRunner, NoIndex as CNoIndex};
use crate::common::comparison::{LogTag, ParityComparator};
use crate::common::status::TEST_PARITY_MODE;
use crate::common::{init_test_logging, parse_output, workspace_root};

// =============================================================================
// TYPE-STATE MARKERS
// =============================================================================

/// Marker for an unindexed Rust runner.
pub struct NoIndex;

/// Marker for an indexed Rust runner.
pub struct Indexed {
    index_path: PathBuf,
    #[allow(dead_code)]
    index_file: risearch::sa::SaIndexFile,
}

// =============================================================================
// RUST RUNNER
// =============================================================================

/// Runner for the Rust risearch implementation.
///
/// Uses type-state pattern to ensure you can only search after indexing.
pub struct RustRunner<S> {
    target_path: PathBuf,
    state: S,
}

impl RustRunner<NoIndex> {
    /// Create a new Rust runner for the given target.
    pub fn new(target: &Path) -> Self {
        Self {
            target_path: target.to_path_buf(),
            state: NoIndex,
        }
    }

    /// Create an index from the target file.
    /// Consumes self and returns an indexed runner.
    /// If `index_path` is None, uses a temporary directory.
    pub fn create_index(self, index_path: Option<&Path>) -> RustRunner<Indexed> {
        let (index_path, _tmpdir) = match index_path {
            Some(p) => (p.to_path_buf(), None),
            None => {
                let tmpdir = tempfile::tempdir().expect("tempdir");
                let path = tmpdir.path().join("target.idx");
                // Leak tempdir to keep index alive
                (path, Some(std::mem::ManuallyDrop::new(tmpdir)))
            }
        };

        risearch::sa::create_suffix_array(&self.target_path, &index_path).expect("build index");
        let index_file = risearch::sa::load_index_file(&index_path).expect("load index");

        RustRunner {
            target_path: self.target_path,
            state: Indexed {
                index_path,
                index_file,
            },
        }
    }
}

impl RustRunner<Indexed> {
    /// Search using risearch as a library.
    pub fn search(
        &self,
        query_path: &Path,
        args: &risearch::args::SearchArgs,
    ) -> Vec<risearch::SearchHit> {
        use risearch::search::SaIndex;

        let index = SaIndex {
            index: &self.state.index_file,
        };
        let queries = risearch::io::read_fasta_sequences(query_path).expect("read query FASTA");
        let hits = risearch::run_search(&queries, &index, args).expect("search");

        // Normalize to 1-based coordinates to match C output format for comparison
        hits.into_iter()
            .map(|mut h| {
                h.q_start += 1;
                h.q_end += 1;
                h
            })
            .collect()
    }

    /// Get the index path.
    #[allow(dead_code)] // Useful API for debugging
    pub fn index_path(&self) -> &Path {
        &self.state.index_path
    }
}

// =============================================================================
// PARITY RUNNER
// =============================================================================

/// Combined runner for parity testing between Rust and C.
///
/// Manages both runners and provides comparison methods.
pub struct ParityRunner {
    rust: RustRunner<Indexed>,
    c: CRunner<crate::common::c_runner::Indexed>,
    #[allow(dead_code)]
    tmpdir: tempfile::TempDir,
}

impl ParityRunner {
    /// Create a new parity runner for the given target.
    /// Creates indexes for both Rust and C implementations.
    pub fn new(target: &Path) -> Self {
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
    pub fn compare(
        &self,
        query: &Path,
        args: &[&str],
    ) -> (Vec<risearch::SearchHit>, Vec<risearch::SearchHit>) {
        let search_args = parse_search_args(args);
        let rust_hits = self.rust.search(query, &search_args);

        // Translate Rust args to C args
        let mut c_args: Vec<&str> = Vec::new();
        let mut iter = args.iter().copied().peekable();
        while let Some(arg) = iter.next() {
            match arg {
                "--no-max-prune" => continue, // Rust-only flag
                "--experimental" => continue, // Rust-only flag
                "--dp-band" => {
                    // Skip value
                    let _ = iter.next();
                    continue;
                }
                "--dp-band-mode" => {
                    let _ = iter.next();
                    continue;
                }
                "-w" | "--wobble" => c_args.push("--noGUseed"), // Map -w to C flag
                _ => c_args.push(arg),
            }
        }
        let c_out = self.c.search(query, &c_args);
        let (c_hits, _) = parse_output(&c_out);

        (rust_hits, c_hits)
    }

    /// Run comparison and assert parity passes.
    pub fn assert_pass(&self, query: &Path, test_name: &str, args: &[&str]) {
        let (rust_recs, c_recs) = self.compare(query, args);

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
        result.log_details(test_name);

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
pub struct SingleSeqRunner {
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
    pub fn new(query_seq: &str, target_seq: &str) -> Self {
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
    pub fn assert_pass(&self, test_name: &str, args: &[&str]) {
        self.runner.assert_pass(&self.query_path, test_name, args);
    }
}

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

/// Parse CLI-style args into SearchArgs using clap.
/// Always includes --no-dedup-shadow for C parity (C doesn't filter contained hits).
pub fn parse_search_args(args: &[&str]) -> risearch::args::SearchArgs {
    use clap::Parser;

    let mut cli_args = vec!["risearch"];
    cli_args.extend(args.iter().copied());
    // C doesn't do shadow dedup, so disable it for parity tests
    cli_args.push("--no-dedup-shadow");

    #[derive(Parser)]
    struct FakeCmd {
        #[command(flatten)]
        search: risearch::args::SearchArgs,
    }

    let parsed = FakeCmd::try_parse_from(&cli_args).expect("Failed to parse search args");
    parsed.search
}
