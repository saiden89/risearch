//! Test runners for parity testing.
//!
//! Provides `RustRunner<S>` and `ParityRunner` for orchestrating
//! parity tests between Rust and C implementations.

use assert_cmd::cargo::cargo_bin_cmd;
use log::{debug, info, trace, warn};
use std::fs;
use std::path::{Path, PathBuf};

use crate::common::c_runner::{CRunner, NoIndex as CNoIndex};
use crate::common::comparison::{LogTag, analyze_hit_pairs};
use crate::common::record::Rec;
use crate::common::status::ParityMode;
use crate::common::{compare_results_with_recs, init_test_logging, parse_output, workspace_root};

// =============================================================================
// TYPE-STATE MARKERS
// =============================================================================

/// Marker for an unindexed Rust runner.
pub struct NoIndex;

/// Marker for an indexed Rust runner.
pub struct Indexed {
    index_path: PathBuf,
    #[allow(dead_code)]
    index_file: risearch::sa::IndexFile,
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
    #[allow(dead_code)] // Useful API - we use create_index_at instead
    pub fn create_index(self) -> RustRunner<Indexed> {
        let tmpdir = tempfile::tempdir().expect("tempdir");
        let index_path = tmpdir.path().join("target.idx");

        risearch::sa::create_suffix_array(&self.target_path, &index_path).expect("build index");
        let index_file = risearch::sa::load_index_file(&index_path).expect("load index");

        // Leak the tempdir to keep the index file alive
        std::mem::forget(tmpdir);

        RustRunner {
            target_path: self.target_path,
            state: Indexed {
                index_path,
                index_file,
            },
        }
    }

    /// Create index at a specific path (for tests that need persistent index).
    pub fn create_index_at(self, index_path: &Path) -> RustRunner<Indexed> {
        risearch::sa::create_suffix_array(&self.target_path, index_path).expect("build index");
        let index_file = risearch::sa::load_index_file(index_path).expect("load index");

        RustRunner {
            target_path: self.target_path,
            state: Indexed {
                index_path: index_path.to_path_buf(),
                index_file,
            },
        }
    }
}

impl RustRunner<Indexed> {
    /// Search using risearch as a library.
    pub fn search(&self, query_path: &Path, args: &risearch::args::SearchArgs) -> Vec<Rec> {
        use risearch::search::SaIndex;

        let index = SaIndex {
            index: &self.state.index_file,
        };
        let queries = risearch::io::read_fasta_sequences(query_path).expect("read query FASTA");
        let hits = risearch::run_search_collect(&queries, &index, args).expect("search");

        hits.iter().map(Rec::from).collect()
    }

    /// Search using CLI subprocess (for trace/debug output).
    pub fn search_cli(&self, query_path: &Path, args: &[&str]) -> std::process::Output {
        let mut cmd = cargo_bin_cmd!("risearch");

        let trace_enabled = std::env::var("RUST_LOG")
            .map(|v| v.contains("trace"))
            .unwrap_or(false);

        let mut final_args = vec!["search"];
        if trace_enabled {
            final_args.push("-vvv");
        }

        final_args.extend_from_slice(&[
            "-i",
            self.state.index_path.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "-",
        ]);
        final_args.extend_from_slice(args);

        let output = cmd.args(&final_args).output().expect("run risearch");

        if trace_enabled {
            eprint!("{}", String::from_utf8_lossy(&output.stderr));
        }

        if !output.status.success() {
            panic!("Rust search failed: {:?}", output.status);
        }
        output
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
        let rust = RustRunner::<NoIndex>::new(target).create_index_at(&rust_index_path);

        // Create C index
        let c_index_path = tmpdir.path().join("target.pksuf");
        let c = CRunner::<CNoIndex>::new(&root).create_index(target, &c_index_path);

        Self { rust, c, tmpdir }
    }

    /// Compare Rust and C results, returning the parsed outputs.
    #[allow(dead_code)] // Useful API for detailed analysis
    pub fn compare(&self, query: &Path, args: &[&str]) -> (Vec<Rec>, Vec<Rec>) {
        let search_args = parse_search_args(args);
        let rust_recs = self.rust.search(query, &search_args);

        let c_args: Vec<&str> = args
            .iter()
            .filter(|&&a| a != "--no-max-prune")
            .cloned()
            .collect();
        let c_out = self.c.search(query, &c_args);
        let (c_recs, _) = parse_output(&c_out);

        (rust_recs, c_recs)
    }

    /// Run comparison and assert parity passes.
    pub fn assert_pass(&self, query: &Path, test_name: &str, args: &[&str]) {
        let search_args = parse_search_args(args);
        let rust_recs = self.rust.search(query, &search_args);

        let c_args: Vec<&str> = args
            .iter()
            .filter(|&&a| a != "--no-max-prune")
            .cloned()
            .collect();
        let c_out = self.c.search(query, &c_args);

        compare_results_with_recs(rust_recs, &c_out, test_name, ParityMode::Relaxed);
    }

    /// Run detailed comparison with CLI output (for debugging).
    pub fn assert_pass_detailed(&self, query: &Path, test_name: &str, args: &[&str]) {
        let rust_output = self.rust.search_cli(query, args);
        let rust_out = String::from_utf8_lossy(&rust_output.stdout).to_string();

        // Log debug traces
        for line in rust_out.lines() {
            if line.contains("DEBUG_TRACE")
                || line.contains("DEBUG_MAX")
                || line.contains("DEBUG_DP_LEFT_RESULT")
                || line.contains("C_DEBUG:")
            {
                trace!("[RUST] {}", line);
            }
        }

        let c_args: Vec<&str> = args
            .iter()
            .filter(|&&a| a != "--no-max-prune")
            .cloned()
            .collect();
        let c_out = self.c.search(query, &c_args);

        let (rust_recs, _) = parse_output(&rust_out);
        let (c_recs, _) = parse_output(&c_out);

        info!(
            "{} {} - Rust={} hits, C={} hits",
            LogTag::Parity,
            test_name,
            rust_recs.len(),
            c_recs.len()
        );

        let mut c_matched = vec![false; c_recs.len()];
        let mut extras: Vec<&Rec> = Vec::new();

        for r in &rust_recs {
            let mut found = false;
            for (i, c) in c_recs.iter().enumerate() {
                if !c_matched[i] {
                    if r == c {
                        c_matched[i] = true;
                        found = true;
                        break;
                    }
                    if r.coords_match(c) {
                        let r_e = r.energy.parse::<f64>().unwrap_or(0.0);
                        let c_e = c.energy.parse::<f64>().unwrap_or(0.0);
                        if r_e < c_e - 0.001 {
                            warn!(
                                "{} Improved energy: Rust E={} vs C E={}",
                                LogTag::Parity,
                                r.energy,
                                c.energy
                            );
                        }
                        if r_e <= c_e + 0.001 {
                            c_matched[i] = true;
                            found = true;
                            break;
                        }
                    }
                }
            }
            if !found {
                extras.push(r);
            }
        }

        let missings: Vec<&Rec> = c_recs
            .iter()
            .enumerate()
            .filter(|(i, _)| !c_matched[*i])
            .map(|(_, r)| r)
            .collect();

        debug!(
            "{} Exact matches: {}, EXTRA: {}, MISSING: {}",
            LogTag::Parity,
            rust_recs.len() - extras.len(),
            extras.len(),
            missings.len()
        );

        analyze_hit_pairs(&extras, &missings);

        if !missings.is_empty() {
            panic!(
                "{} FAILED {}: {} missings ({} extras)",
                LogTag::Parity,
                test_name,
                missings.len(),
                extras.len()
            );
        }

        if !extras.is_empty() {
            warn!(
                "{} {} extras in {} (acceptable if 0 missings)",
                LogTag::Parity,
                extras.len(),
                test_name
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
pub struct SingleSeqRunner {
    query_path: PathBuf,
    #[allow(dead_code)] // Kept for potential future use
    target_path: PathBuf,
    runner: ParityRunner,
    #[allow(dead_code)]
    tmpdir: tempfile::TempDir,
}

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

    /// Run detailed comparison with CLI output.
    pub fn assert_pass_detailed(&self, test_name: &str, args: &[&str]) {
        self.runner
            .assert_pass_detailed(&self.query_path, test_name, args);
    }
}

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

/// Parse CLI-style args into SearchArgs using clap.
pub fn parse_search_args(args: &[&str]) -> risearch::args::SearchArgs {
    use clap::Parser;

    let mut cli_args = vec!["risearch"];
    cli_args.extend(args.iter().copied());

    #[derive(Parser)]
    struct FakeCmd {
        #[command(flatten)]
        search: risearch::args::SearchArgs,
    }

    let parsed = FakeCmd::try_parse_from(&cli_args).expect("Failed to parse search args");
    parsed.search
}
