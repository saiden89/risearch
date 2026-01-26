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
        args: &risearch::config::SearchArgs,
    ) -> Vec<risearch::SearchHit> {
        use risearch::search::SaIndex;

        let index = SaIndex {
            index: &self.state.index_file,
        };
        let queries = risearch::sa::process_sequences(query_path)
            .expect("read query FASTA")
            .sequences
            .into_iter()
            .map(|s| (s.name, s.sequence))
            .collect::<Vec<_>>();
        let hits = risearch::search::run_search(&queries, &index, args).expect("search");

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
        let mut c_args: Vec<String> = Vec::new();
        let mut seed_start: Option<&str> = None;
        let mut seed_end: Option<&str> = None;
        let mut seed_length: Option<&str> = None;
        let mut has_legacy_seed = false;
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
                "-s" | "--seed" => {
                    has_legacy_seed = true;
                    if let Some(val) = iter.next() {
                        c_args.push("-s".to_string());
                        c_args.push(val.to_string());
                    }
                }
                _ if arg.starts_with("--seed=") => {
                    has_legacy_seed = true;
                    let val = &arg["--seed=".len()..];
                    c_args.push("-s".to_string());
                    c_args.push(val.to_string());
                    continue;
                }
                "--seed-start" => {
                    if let Some(val) = iter.next() {
                        seed_start = Some(val);
                    }
                    continue;
                }
                "--seed-end" => {
                    if let Some(val) = iter.next() {
                        seed_end = Some(val);
                    }
                    continue;
                }
                "--seed-length" => {
                    if let Some(val) = iter.next() {
                        seed_length = Some(val);
                    }
                    continue;
                }
                _ if arg.starts_with("--seed-start=") => {
                    seed_start = Some(&arg["--seed-start=".len()..]);
                    continue;
                }
                _ if arg.starts_with("--seed-end=") => {
                    seed_end = Some(&arg["--seed-end=".len()..]);
                    continue;
                }
                _ if arg.starts_with("--seed-length=") => {
                    seed_length = Some(&arg["--seed-length=".len()..]);
                    continue;
                }
                "-U" | "--no-guseed" | "--noGUseed" => c_args.push("--noGUseed".to_string()),
                "--seed-pairing" => {
                    if let Some(val) = iter.next() {
                        if val == "strict" {
                            c_args.push("--noGUseed".to_string());
                        }
                    }
                }
                _ => c_args.push(arg.to_string()),
            }
        }
        if !has_legacy_seed {
            let spec = match (seed_start, seed_end, seed_length) {
                (Some(start), Some(end), Some(len)) => Some(format!("{}:{}/{}", start, end, len)),
                (Some(start), Some(end), None) => Some(format!("{}:{}", start, end)),
                (None, None, Some(len)) => Some(len.to_string()),
                _ => None,
            };
            if let Some(s) = spec {
                c_args.push("-s".to_string());
                c_args.push(s);
            }
        }
        let c_args_ref: Vec<&str> = c_args.iter().map(|s| s.as_str()).collect();
        let c_out = self.c.search(query, &c_args_ref);
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
pub fn parse_search_args(args: &[&str]) -> risearch::config::SearchArgs {
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
    let has_pairing = args
        .iter()
        .any(|arg| *arg == "--seed-pairing" || arg.starts_with("--seed-pairing="));
    let has_no_guseed = args
        .iter()
        .any(|arg| *arg == "-U" || *arg == "--no-guseed" || *arg == "--noGUseed");
    if !has_pairing && !has_no_guseed {
        cli_args.push("--seed-pairing".into());
        cli_args.push("allow_wobble".into());
    }
    // C doesn't do shadow dedup, so disable it for parity tests
    cli_args.push("--no-dedup-shadow".into());

    #[derive(Parser)]
    struct FakeCmd {
        #[command(flatten)]
        search: risearch::cli_args::SearchArgs,
    }

    let parsed = FakeCmd::try_parse_from(&cli_args).expect("Failed to parse search args");
    let mut search: risearch::config::SearchArgs = parsed.search.into();
    let explicit_pairing = args
        .iter()
        .any(|arg| *arg == "--seed-pairing" || arg.starts_with("--seed-pairing="));
    search.seed.apply_mismatch_overrides();
    search.seed.apply_pairing_overrides(explicit_pairing);
    search
}
