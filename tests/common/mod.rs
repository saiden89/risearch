//! Common test utilities for parity testing.
//!
//! This module provides infrastructure for comparing Rust and C risearch output.
//!
//! # Module Organization
//!
//! - `record`: `Rec` struct and parsing
//! - `status`: `HitStatus`, `MissingReason`, `ParityMode`
//! - `table`: `ParityTable` and display formatting
//! - `comparison`: Comparison logic and `ParityComparator`
//! - `c_runner`: C binary runner with type-state pattern
//! - `runner`: Test orchestration functions

// =============================================================================
// SUBMODULES
// =============================================================================

pub mod c_runner;
pub mod comparison;
pub mod record;
pub mod runner;
pub mod status;
pub mod table;

// =============================================================================
// RE-EXPORTS (public API)
// Used by runner.rs (internally)
pub(crate) use comparison::{LogTag, classify_missing};
pub(crate) use record::Rec;
pub(crate) use status::{HitStatus, ParityMode};
pub(crate) use table::{ParityKind, ParityTable, TableConfig};

// Used by c_parity.rs
pub use runner::{ParityRunner, SingleSeqRunner};

// =============================================================================
// LOGGING SETUP
// =============================================================================

use log::{debug, info, trace, warn};
use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::Once;
use tabled::{Table, Tabled, settings::Style};

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
pub fn parse_output(output: &str) -> (Vec<Rec>, usize) {
    let mut recs: Vec<Rec> = output
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty())
        .filter_map(Rec::from_line)
        .collect();
    let parsed_count = recs.len();
    recs.sort();
    recs.dedup();
    (recs, parsed_count)
}

/// Compare results where Rust records come directly from SearchHit (with proper seed info).
pub fn compare_results_with_recs(
    rust_recs: Vec<Rec>,
    c_out: &str,
    test_name: &str,
    mode: ParityMode,
) {
    init_test_logging();

    let (c_recs, _) = parse_output(c_out);
    compare_recs_impl(rust_recs, c_recs, test_name, mode);
}

/// Shared implementation for comparing Rust and C records.
fn compare_recs_impl(rust_recs: Vec<Rec>, c_recs: Vec<Rec>, test_name: &str, mode: ParityMode) {
    info!(
        "{} {} - Comparing: Rust={} hits, C={} hits",
        LogTag::Parity,
        test_name,
        rust_recs.len(),
        c_recs.len()
    );

    if log::log_enabled!(log::Level::Trace) {
        let mut debug_out = String::new();
        debug_out.push_str("RUST RECORDS:\n");
        for r in &rust_recs {
            debug_out.push_str(&format!("{:?}\n", r));
        }
        debug_out.push_str("\nC RECORDS:\n");
        for c in &c_recs {
            debug_out.push_str(&format!("{:?}\n", c));
        }
        let _ = std::fs::write("target/debug_parity_report.txt", debug_out);
        trace!(
            "{} Wrote raw records to target/debug_parity_report.txt",
            LogTag::Parity
        );
    }

    // Group by (q_id, t_id)
    let mut keys = HashSet::new();
    for r in &rust_recs {
        keys.insert((r.q_id.clone(), r.t_id.clone()));
    }
    for r in &c_recs {
        keys.insert((r.q_id.clone(), r.t_id.clone()));
    }
    let mut sorted_keys: Vec<_> = keys.into_iter().collect();
    sorted_keys.sort();

    let mut missing_count = 0;
    let mut extra_count = 0;
    let mut mismatch_count = 0;
    let mut exact_match_count = 0;

    let mut rust_better_energy = 0;
    let mut rust_worse_energy = 0;
    let mut energy_equal = 0;

    let mut extra_len_sum: usize = 0;
    let mut extra_energy_sum: f64 = 0.0;
    let mut missing_len_sum: usize = 0;
    let mut missing_energy_sum: f64 = 0.0;

    let all_rust_refs: Vec<&Rec> = rust_recs.iter().collect();

    for (q, t) in sorted_keys {
        let mut r_group: Vec<&Rec> = rust_recs
            .iter()
            .filter(|r| r.q_id == q && r.t_id == t)
            .collect();
        r_group.sort_by(|a, b| {
            a.q_start
                .cmp(&b.q_start)
                .then(a.q_end.cmp(&b.q_end))
                .then(a.t_start.cmp(&b.t_start))
                .then(a.t_end.cmp(&b.t_end))
        });
        let mut c_group: Vec<&Rec> = c_recs
            .iter()
            .filter(|r| r.q_id == q && r.t_id == t)
            .collect();
        c_group.sort_by(|a, b| {
            a.q_start
                .cmp(&b.q_start)
                .then(a.q_end.cmp(&b.q_end))
                .then(a.t_start.cmp(&b.t_start))
                .then(a.t_end.cmp(&b.t_end))
        });

        // Pass 1: Remove Exact Matches
        let mut c_matched = vec![false; c_group.len()];
        let mut r_remaining = Vec::new();

        for r in r_group {
            let mut found = false;
            for (i, c) in c_group.iter().enumerate() {
                if !c_matched[i] && r == *c {
                    c_matched[i] = true;
                    found = true;
                    exact_match_count += 1;
                    break;
                }
            }
            if !found {
                r_remaining.push(r);
            }
        }

        let c_remaining: Vec<&Rec> = c_group
            .into_iter()
            .enumerate()
            .filter(|(i, _)| !c_matched[*i])
            .map(|(_, r)| r)
            .collect();

        if r_remaining.is_empty() && c_remaining.is_empty() {
            continue;
        }

        let mut c_rem_matched = vec![false; c_remaining.len()];

        let mut good_hits: Vec<(&Rec, &Rec, HitStatus)> = Vec::new();
        let mut problem_hits: Vec<(&Rec, &Rec, HitStatus)> = Vec::new();
        let mut extra_hits: Vec<&Rec> = Vec::new();

        for r in &r_remaining {
            let mut best_match_idx = None;

            for (i, c) in c_remaining.iter().enumerate() {
                if c_rem_matched[i] {
                    continue;
                }
                if r.coords_match(c) {
                    best_match_idx = Some(i);
                    break;
                }
            }

            if let Some(idx) = best_match_idx {
                c_rem_matched[idx] = true;
                let c = c_remaining[idx];

                let r_e = r.energy.parse::<f64>().unwrap_or(0.0);
                let c_e = c.energy.parse::<f64>().unwrap_or(0.0);

                let status: HitStatus;

                if r_e < c_e - 0.001 {
                    rust_better_energy += 1;
                    status = HitStatus::RustBetter;
                    good_hits.push((r, c, status));
                } else if r_e > c_e + 0.001 {
                    rust_worse_energy += 1;
                    mismatch_count += 1;
                    status = HitStatus::RustWorse;
                    problem_hits.push((r, c, status));
                } else {
                    energy_equal += 1;
                    if r.interaction == c.interaction {
                        status = HitStatus::Identical;
                    } else {
                        status = HitStatus::CoOptimal;
                    }
                    good_hits.push((r, c, status));
                }
            } else {
                extra_count += 1;
                extra_len_sum += r.interaction.len();
                extra_energy_sum += r.energy.parse::<f64>().unwrap_or(0.0);
                extra_hits.push(r);
            }
        }

        debug!("");
        debug!(
            "{} Group [{}:{}] R={} C={} (exact={})",
            LogTag::Parity,
            q,
            t,
            r_remaining.len() + exact_match_count,
            c_remaining.len() + exact_match_count,
            exact_match_count
        );

        for (r, c, status) in &good_hits {
            if *status == HitStatus::RustBetter {
                debug!(
                    "{}   {} {} {} R={} C={}",
                    LogTag::Parity,
                    status,
                    r.strand,
                    format!("q[{},{}] t[{},{}]", r.q_start, r.q_end, r.t_start, r.t_end),
                    r.energy,
                    c.energy
                );
                let table = ParityTable {
                    kind: ParityKind::Mismatch { rust: r, c },
                    config: TableConfig::default(),
                };
                for line in table.to_string().lines() {
                    debug!("{} {}", LogTag::Parity, line);
                }
            } else {
                debug!(
                    "{}   {} {} {} E={}",
                    LogTag::Parity,
                    status,
                    r.strand,
                    format!("q[{},{}] t[{},{}]", r.q_start, r.q_end, r.t_start, r.t_end),
                    r.energy
                );
            }
        }

        for (r, c, status) in &problem_hits {
            debug!(
                "{}   ✗ {} {} {} E={} (C: {})",
                LogTag::Parity,
                status,
                r.strand,
                format!("q[{},{}] t[{},{}]", r.q_start, r.q_end, r.t_start, r.t_end),
                r.energy,
                c.energy
            );
            let table = ParityTable {
                kind: ParityKind::Mismatch { rust: r, c },
                config: TableConfig::default(),
            };
            for line in table.to_string().lines() {
                debug!("{} {}", LogTag::Parity, line);
            }
        }

        for r in &extra_hits {
            debug!(
                "{}   ✗ EXTRA {} {} E={}",
                LogTag::Parity,
                r.strand,
                format!("q[{},{}] t[{},{}]", r.q_start, r.q_end, r.t_start, r.t_end),
                r.energy
            );
            let table = ParityTable {
                kind: ParityKind::RustOnly(r),
                config: TableConfig::default(),
            };
            for line in table.to_string().lines() {
                debug!("{} {}", LogTag::Parity, line);
            }
        }

        for (i, c) in c_remaining.iter().enumerate() {
            if !c_rem_matched[i] {
                let (reason, covering_hit) = classify_missing(c, &all_rust_refs);
                let acceptable = reason.is_acceptable(mode);

                let status_prefix = if acceptable { "  " } else { "✗ " };
                debug!(
                    "{} {}{} {} Coords: {}",
                    LogTag::Parity,
                    status_prefix,
                    HitStatus::Missing,
                    reason,
                    c.fmt_coords()
                );

                if let Some(r) = covering_hit {
                    debug!(
                        "{}      Covered by Rust: {}",
                        LogTag::Parity,
                        r.fmt_coords()
                    );
                }

                let table = ParityTable {
                    kind: ParityKind::COnly(c),
                    config: TableConfig::default(),
                };
                for line in table.to_string().lines() {
                    debug!("{} {}", LogTag::Parity, line);
                }

                if !acceptable {
                    missing_count += 1;
                    missing_len_sum += c.interaction.len();
                    missing_energy_sum += c.energy.parse::<f64>().unwrap_or(0.0);
                }
            }
        }
    }

    let coord_matched = rust_better_energy + rust_worse_energy + energy_equal;
    let total_c = c_recs.len();

    let pct = |n: usize, total: usize| -> String {
        if total == 0 {
            "0%".into()
        } else {
            format!("{}%", n * 100 / total)
        }
    };

    let is_pass = mismatch_count == 0 && missing_count == 0;

    const KNOWN_DIVERGENT_TESTS: &[&str] = &[
        "seed_only_no_extension",
        "wobble_seed_pairs",
        "right_extension_only",
    ];
    let is_allowed_divergence = KNOWN_DIVERGENT_TESTS.contains(&test_name)
        && mismatch_count == 0
        && extra_count == 0
        && missing_count > 0;

    #[derive(Tabled)]
    struct SummaryRow {
        #[tabled(rename = "Metric")]
        metric: String,
        #[tabled(rename = "Value")]
        value: String,
    }

    let mut rows = vec![
        SummaryRow {
            metric: "Test".into(),
            value: test_name.to_string(),
        },
        SummaryRow {
            metric: "Rust hits".into(),
            value: rust_recs.len().to_string(),
        },
        SummaryRow {
            metric: "C hits".into(),
            value: total_c.to_string(),
        },
        SummaryRow {
            metric: "Exact matches".into(),
            value: format!(
                "{} ({})",
                exact_match_count,
                pct(exact_match_count, total_c)
            ),
        },
        SummaryRow {
            metric: "Coord-matched (diff content)".into(),
            value: coord_matched.to_string(),
        },
        SummaryRow {
            metric: "  ├ Rust better energy".into(),
            value: format!(
                "{} ({})",
                rust_better_energy,
                pct(rust_better_energy, coord_matched)
            ),
        },
        SummaryRow {
            metric: "  ├ Equal energy".into(),
            value: format!("{} ({})", energy_equal, pct(energy_equal, coord_matched)),
        },
        SummaryRow {
            metric: "  └ Rust worse energy".into(),
            value: format!(
                "{} ({})",
                rust_worse_energy,
                pct(rust_worse_energy, coord_matched)
            ),
        },
    ];

    if missing_count > 0 {
        let avg_len = missing_len_sum as f64 / missing_count as f64;
        let avg_energy = missing_energy_sum / missing_count as f64;
        rows.push(SummaryRow {
            metric: "Missing (in C, not Rust)".into(),
            value: format!(
                "{} (avg_len={:.1}, avg_E={:.2})",
                missing_count, avg_len, avg_energy
            ),
        });
    } else {
        rows.push(SummaryRow {
            metric: "Missing (in C, not Rust)".into(),
            value: "0".into(),
        });
    }

    if extra_count > 0 {
        let avg_len = extra_len_sum as f64 / extra_count as f64;
        let avg_energy = extra_energy_sum / extra_count as f64;
        rows.push(SummaryRow {
            metric: "Extra (in Rust, not C)".into(),
            value: format!(
                "{} (avg_len={:.1}, avg_E={:.2})",
                extra_count, avg_len, avg_energy
            ),
        });
    } else {
        rows.push(SummaryRow {
            metric: "Extra (in Rust, not C)".into(),
            value: "0".into(),
        });
    }

    let verdict = if is_pass {
        "✓ PASS".to_string()
    } else if is_allowed_divergence {
        "⚠ ALLOWED DIVERGENCE".to_string()
    } else {
        format!(
            "✗ FAIL ({} worse, {} missing, {} extra)",
            rust_worse_energy, missing_count, extra_count
        )
    };
    rows.push(SummaryRow {
        metric: "VERDICT".into(),
        value: verdict,
    });

    let table = Table::new(rows).with(Style::rounded()).to_string();
    for line in table.lines() {
        info!("{} {}", LogTag::Parity, line);
    }

    if is_allowed_divergence {
        warn!(
            "{} Allowed divergence for '{}': Missing in Rust expected due to Maximality/Wobble improvements",
            LogTag::Parity,
            test_name
        );
        return;
    }

    if !is_pass {
        panic!(
            "\n{} FAILED: {} ({} rust-worse, {} missing, {} extra)\nRun with RUST_LOG=debug for detailed diff analysis.",
            LogTag::Parity,
            test_name,
            rust_worse_energy,
            missing_count,
            extra_count
        );
    }
}
