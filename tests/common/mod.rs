//! Common test utilities for parity testing.
//!
//! This module provides infrastructure for comparing Rust and C risearch output.
//!
//! # Module Organization
//!
//! - `status`: `HitStatus`, `MissingReason`, `ParityMode`
//! - `table`: Table formatting utilities
//! - `comparison`: Comparison logic and `ParityComparator`
//! - `c_runner`: C binary runner with type-state pattern
//! - `runner`: Test orchestration functions

// =============================================================================
// SUBMODULES
// =============================================================================

pub mod c_runner;
pub mod comparison;
pub mod runner;
pub mod status;
pub mod table;

// =============================================================================
// RE-EXPORTS (public API)
// =============================================================================

// Used by runner.rs (internally)

// Used by c_parity.rs
pub use runner::{ParityRunner, SingleSeqRunner};

// =============================================================================
// LOGGING SETUP
// =============================================================================

use std::path::PathBuf;
use std::sync::Once;

static INIT: Once = Once::new();

pub fn init_test_logging() {
    INIT.call_once(|| {
        let _ = env_logger::builder()
            .is_test(true)
            .format(|buf, record| {
                use std::io::Write;
                let style = buf.default_level_style(record.level());
                writeln!(
                    buf,
                    "[{}{}{}] {}",
                    style.render(),
                    record.level(),
                    style.render_reset(),
                    record.args()
                )
            })
            .try_init();
    });
}

// =============================================================================
// OUTPUT LEVEL GUIDE
// =============================================================================
// INFO  - Test progress: test name, hit counts, pass/fail summary
// DEBUG - Detailed comparisons: individual hit diffs, alignment analysis, coords
// TRACE - Raw data: full record dumps, C_DEBUG traces, file writes
// WARN  - Known divergences, acceptable differences
// ERROR - Reserved for panics (test failures)
// =============================================================================

pub fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Returns (records after dedup, count before dedup)
pub fn parse_output(output: &str) -> (Vec<risearch::SearchHit>, usize) {
    let mut hits: Vec<risearch::SearchHit> = output
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty())
        .filter_map(risearch::SearchHit::from_c_output)
        .collect();
    let parsed_count = hits.len();
    // Sort by group_key, coords for deterministic comparison
    hits.sort_by(|a, b| {
        a.group_key()
            .cmp(&b.group_key())
            .then(a.q_start.cmp(&b.q_start))
            .then(a.output_t_start.cmp(&b.output_t_start))
    });
    hits.dedup_by(|a, b| {
        a.coords_match(b) && (a.energy.as_f64() - b.energy.as_f64()).abs() < 0.001
    });
    (hits, parsed_count)
}

// Comparison logic has been moved to comparison.rs (ParityComparator)
