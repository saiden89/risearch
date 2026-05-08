//! Shared test utilities for parity and Rust-only integration testing.
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

pub(crate) mod c_runner;
pub(crate) mod comparison;
pub(crate) mod runner;
pub(crate) mod search_hit;
pub(crate) mod status;
pub(crate) mod table;

// =============================================================================
// RE-EXPORTS (public API)
// =============================================================================

// Used by runner.rs (internally)

// Used by c_parity.rs (not necessarily used in every test crate)
#[allow(unused_imports)]
pub(crate) use runner::{ParityRunner, SingleSeqRunner};

// Re-export test extensions for SearchHit
pub(crate) use search_hit::{parse_bindingsite_output, SearchHitExt};

// =============================================================================
// LOGGING SETUP
// =============================================================================

use std::path::PathBuf;
use std::sync::Once;

static INIT: Once = Once::new();

pub(crate) fn init_test_logging() {
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

pub(crate) fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Returns (records after dedup, count before dedup)
pub(crate) fn parse_output(
    output: &str,
    query_registry: &risearch::QueryRegistry,
    target_store: &risearch::TargetStore,
) -> (Vec<risearch::SearchHit>, usize) {
    let mut hits: Vec<risearch::SearchHit> = output
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty())
        .filter_map(|l| parse_bindingsite_output(l, query_registry, target_store))
        .collect();
    let parsed_count = hits.len();
    // Sort by group_key, then all coordinates for deterministic dedup
    // Must sort by ALL coordinate fields to ensure true duplicates are adjacent
    hits.sort_by(|a, b| {
        (a.query_idx, a.target_idx)
            .cmp(&(b.query_idx, b.target_idx))
            .then(a.q_start.cmp(&b.q_start))
            .then(a.q_end.cmp(&b.q_end))
            .then(a.t_start.cmp(&b.t_start))
            .then(a.t_end.cmp(&b.t_end))
    });
    hits.dedup_by(|a, b| {
        a.coords_match(b) && (f64::from(a.energy) - f64::from(b.energy)).abs() < 0.001
    });
    (hits, parsed_count)
}

// Comparison logic has been moved to comparison.rs (ParityComparator)
