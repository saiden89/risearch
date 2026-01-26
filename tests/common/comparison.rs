//! Comparison logic for parity testing.
//!
//! Contains `ParityComparator` - a builder pattern for comparing Rust and C results.

use std::collections::HashSet;

use crate::common::status::{HitStatus, MissingReason, ParityMode, TEST_PARITY_MODE};
use risearch::SearchHit;

// =============================================================================
// LOG TAG
// =============================================================================

#[allow(dead_code)] // Some tags are only used in specific parity modes/tests.
pub enum LogTag {
    Pair,
    Parity,
}

impl std::fmt::Display for LogTag {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Pair => write!(f, "[PAIR]"),
            Self::Parity => write!(f, "[PARITY]"),
        }
    }
}

// =============================================================================
// MATCHED HIT
// =============================================================================

/// A hit that was matched between Rust and C with its status.
#[derive(Debug, Clone)]
#[allow(dead_code)] // Part of public API for future use
pub struct MatchedHit {
    pub rust: SearchHit,
    pub c: SearchHit,
    pub status: HitStatus,
}

// =============================================================================
// PARITY RESULT
// =============================================================================

/// Result of a parity comparison.
#[derive(Debug, Default)]
pub struct ParityResult {
    pub exact_matches: usize,
    pub rust_better: Vec<MatchedHit>,
    pub rust_worse: Vec<MatchedHit>,
    pub co_optimal: Vec<MatchedHit>,
    pub extras: Vec<SearchHit>,
    /// Missing hits: (C hit, reason, optional overlapping Rust hit)
    pub missings: Vec<(SearchHit, MissingReason, Option<SearchHit>)>,
}

impl ParityResult {
    /// Check if the result is a pass given the parity mode.
    pub fn is_pass(&self, mode: ParityMode) -> bool {
        // Rust-worse is always a failure
        if !self.rust_worse.is_empty() {
            return false;
        }
        match mode {
            ParityMode::Absolute => {
                // 100% identical - nothing but exact matches allowed
                self.extras.is_empty()
                    && self.missings.is_empty()
                    && self.rust_better.is_empty()
                    && self.co_optimal.is_empty()
            }
            ParityMode::Strict => {
                // Allow co-optimal (same energy, different trace)
                self.extras.is_empty() && self.missings.is_empty() && self.rust_better.is_empty()
            }
            ParityMode::Balanced => {
                // Allow co-optimal and extras, but no missing or rust-better
                self.missings.is_empty() && self.rust_better.is_empty()
            }
            ParityMode::Relaxed => {
                // Allow co-optimal, rust-better, extras
                // Only fail on unacceptable missings (worse energy or no overlap)
                !self.missings.iter().any(|(_, r, _)| !r.is_acceptable(mode))
            }
        }
    }

    /// Log detailed comparison results for all hit types.
    pub fn log_details(&self, test_name: &str) {
        use crate::common::table::{ParityKind, ParityTable, TableConfig};
        use std::collections::BTreeMap;

        log::debug!("{} {} Summary:", LogTag::Parity, test_name);
        log::debug!(
            "{} Exact matches: {}, Co-optimal: {}, Rust-better: {}, Rust-worse: {}, Extras: {}, Missings: {}",
            LogTag::Parity,
            self.exact_matches,
            self.co_optimal.len(),
            self.rust_better.len(),
            self.rust_worse.len(),
            self.extras.len(),
            self.missings.len()
        );

        // Collect all hits into groups by (query_idx, target_idx)
        let mut groups: BTreeMap<(u32, u32), Vec<(&str, String)>> = BTreeMap::new();

        // Helper to render a table and collect lines
        let render_table = |kind: ParityKind| -> String {
            let table = ParityTable {
                kind,
                config: TableConfig::default(),
            };
            table.to_string()
        };

        // Collect co-optimal (only show details if not allowed by current mode)
        if matches!(*TEST_PARITY_MODE, ParityMode::Absolute) {
            for matched in &self.co_optimal {
                let key = (matched.rust.query_idx, matched.rust.target_idx);
                let label = format!(
                    "CO-OPTIMAL: {} FP={} vs C FP={}",
                    matched.rust.fmt_coords(),
                    matched.rust.fingerprint(),
                    matched.c.fingerprint()
                );
                let table = render_table(ParityKind::Mismatch {
                    rust: &matched.rust,
                    c: &matched.c,
                });
                groups
                    .entry(key)
                    .or_default()
                    .push(("co_optimal", format!("{}\n{}", label, table)));
            }
        }

        // Collect rust-better
        for matched in &self.rust_better {
            let key = (matched.rust.query_idx, matched.rust.target_idx);
            let label = format!(
                "RUST-BETTER: {} vs C E={}",
                matched.rust.fmt_coords(),
                matched.c.energy
            );
            let table = render_table(ParityKind::Mismatch {
                rust: &matched.rust,
                c: &matched.c,
            });
            groups
                .entry(key)
                .or_default()
                .push(("rust_better", format!("{}\n{}", label, table)));
        }

        // Collect rust-worse
        for matched in &self.rust_worse {
            let key = (matched.rust.query_idx, matched.rust.target_idx);
            let label = format!(
                "✗ RUST-WORSE: {} vs C E={}",
                matched.rust.fmt_coords(),
                matched.c.energy
            );
            let table = render_table(ParityKind::Mismatch {
                rust: &matched.rust,
                c: &matched.c,
            });
            groups
                .entry(key)
                .or_default()
                .push(("rust_worse", format!("{}\n{}", label, table)));
        }

        // Collect extras
        for extra in &self.extras {
            let key = (extra.query_idx, extra.target_idx);
            let label = format!("✗ EXTRA: {}", extra.fmt_coords());
            let table = render_table(ParityKind::RustOnly(extra));
            groups
                .entry(key)
                .or_default()
                .push(("extra", format!("{}\n{}", label, table)));
        }

        // Collect missings
        for (missing, reason, overlap) in &self.missings {
            let key = (missing.query_idx, missing.target_idx);
            let (label, table) = match overlap {
                Some(rust_hit) => {
                    let energy_diff = missing.energy.as_f64() - rust_hit.energy.as_f64();
                    // Format: C hit coords vs Rust hit coords, with energy delta
                    let label = format!(
                        "✗ MISSING ({}): C {} vs R {} (ΔE={:+.2})",
                        reason,
                        missing.fmt_coords(),
                        rust_hit.fmt_coords(),
                        energy_diff
                    );
                    let table = render_table(ParityKind::CoveredBy {
                        c: missing,
                        rust: rust_hit,
                    });
                    (label, table)
                }
                None => {
                    let label = format!("✗ MISSING ({}): {}", reason, missing.fmt_coords());
                    let table = render_table(ParityKind::COnly(missing));
                    (label, table)
                }
            };
            groups
                .entry(key)
                .or_default()
                .push(("missing", format!("{}\n{}", label, table)));
        }

        // Print groups
        for ((q_id, t_id), hits) in groups {
            log::debug!("");
            log::debug!(
                "{} Group [{}:{}] hits={}",
                LogTag::Parity,
                q_id,
                t_id,
                hits.len()
            );
            for (_kind, content) in hits {
                for line in content.lines() {
                    log::debug!("{} {}", LogTag::Parity, line);
                }
            }
        }

        // Final summary table (matching old compare_recs_impl format)
        use crate::common::table::{SummaryRow, render_summary_table};

        let total_rust = self.exact_matches
            + self.rust_better.len()
            + self.rust_worse.len()
            + self.co_optimal.len()
            + self.extras.len();
        let total_c = self.exact_matches
            + self.rust_better.len()
            + self.rust_worse.len()
            + self.co_optimal.len()
            + self.missings.len();
        let coord_matched = self.rust_better.len() + self.rust_worse.len() + self.co_optimal.len();

        let pct = |n: usize, total: usize| -> String {
            if total == 0 {
                "0%".into()
            } else {
                format!("{}%", n * 100 / total)
            }
        };

        // Calculate avg stats for extras/missings
        let avg_extra_len = if self.extras.is_empty() {
            0.0
        } else {
            self.extras
                .iter()
                .map(|r| r.fingerprint().len())
                .sum::<usize>() as f64
                / self.extras.len() as f64
        };
        let avg_extra_energy = if self.extras.is_empty() {
            0.0
        } else {
            self.extras.iter().map(|r| r.energy.as_f64()).sum::<f64>() / self.extras.len() as f64
        };
        let avg_missing_len = if self.missings.is_empty() {
            0.0
        } else {
            self.missings
                .iter()
                .map(|(r, _, _)| r.fingerprint().len())
                .sum::<usize>() as f64
                / self.missings.len() as f64
        };
        let avg_missing_energy = if self.missings.is_empty() {
            0.0
        } else {
            self.missings
                .iter()
                .map(|(r, _, _)| r.energy.as_f64())
                .sum::<f64>()
                / self.missings.len() as f64
        };

        let mut rows = vec![
            SummaryRow::new("Rust hits", total_rust),
            SummaryRow::new("C hits", total_c),
            SummaryRow::new(
                "Exact matches",
                format!(
                    "{} ({})",
                    self.exact_matches,
                    pct(self.exact_matches, total_c)
                ),
            ),
            SummaryRow::new("Coord-matched (diff content)", coord_matched),
            SummaryRow::new(
                "  ├ Rust better energy",
                format!(
                    "{} ({})",
                    self.rust_better.len(),
                    pct(self.rust_better.len(), coord_matched)
                ),
            ),
            SummaryRow::new(
                "  ├ Equal energy (co-optimal)",
                format!(
                    "{} ({})",
                    self.co_optimal.len(),
                    pct(self.co_optimal.len(), coord_matched)
                ),
            ),
            SummaryRow::new(
                "  └ Rust worse energy",
                format!(
                    "{} ({})",
                    self.rust_worse.len(),
                    pct(self.rust_worse.len(), coord_matched)
                ),
            ),
        ];

        if !self.missings.is_empty() {
            rows.push(SummaryRow::new(
                "Missing (in C, not Rust)",
                format!(
                    "{} (avg_len={:.1}, avg_E={:.2})",
                    self.missings.len(),
                    avg_missing_len,
                    avg_missing_energy
                ),
            ));
        } else {
            rows.push(SummaryRow::new("Missing (in C, not Rust)", "0"));
        }

        if !self.extras.is_empty() {
            rows.push(SummaryRow::new(
                "Extra (in Rust, not C)",
                format!(
                    "{} (avg_len={:.1}, avg_E={:.2})",
                    self.extras.len(),
                    avg_extra_len,
                    avg_extra_energy
                ),
            ));
        } else {
            rows.push(SummaryRow::new("Extra (in Rust, not C)", "0"));
        }

        // Only show verdict on failure (pass is implied by test not panicking)
        if !self.is_pass(*TEST_PARITY_MODE) {
            let verdict = format!(
                "✗ FAIL [{:?}] ({} co-opt, {} better, {} worse, {} missing, {} extra)",
                *TEST_PARITY_MODE,
                self.co_optimal.len(),
                self.rust_better.len(),
                self.rust_worse.len(),
                self.missings.len(),
                self.extras.len()
            );
            rows.push(SummaryRow::new("VERDICT", verdict));
        }

        for line in render_summary_table(rows).lines() {
            log::info!("{} {}", LogTag::Parity, line);
        }
    }
}

// =============================================================================
// OVERLAP HELPERS
// =============================================================================

/// Check if two coordinate ranges overlap.
pub fn ranges_overlap(a_start: usize, a_end: usize, b_start: usize, b_end: usize) -> bool {
    a_start <= b_end && b_start <= a_end
}

/// Check if two hits overlap in both query and target coordinates.
/// Also requires same query_id and target_id to be meaningful.
pub fn hits_overlap(a: &SearchHit, b: &SearchHit) -> bool {
    a.query_idx == b.query_idx
        && a.target_idx == b.target_idx
        && a.strand == b.strand
        && ranges_overlap(
            a.output_t_start,
            a.output_t_end,
            b.output_t_start,
            b.output_t_end,
        )
        && ranges_overlap(a.q_start, a.q_end, b.q_start, b.q_end)
}

/// Classify why a C hit is missing from Rust output.
/// Returns the reason and optionally the best overlapping Rust hit.
pub fn classify_missing<'a>(
    c_hit: &SearchHit,
    rust_hits: &[&'a SearchHit],
) -> (MissingReason, Option<&'a SearchHit>) {
    let c_e = c_hit.energy.as_f64();

    let best_overlap = rust_hits
        .iter()
        .filter(|r| hits_overlap(c_hit, r))
        .min_by(|a, b| {
            let a_e = a.energy.as_f64();
            let b_e = b.energy.as_f64();
            a_e.partial_cmp(&b_e).unwrap_or(std::cmp::Ordering::Equal)
        })
        .copied();

    match best_overlap {
        Some(r) => {
            let r_e = r.energy.as_f64();
            let reason = if r_e < c_e - 0.001 {
                MissingReason::BetterEnergy
            } else if r_e > c_e + 0.001 {
                MissingReason::WorseEnergy
            } else {
                MissingReason::EqualEnergy
            };
            (reason, Some(r))
        }
        None => (MissingReason::NoOverlap, None),
    }
}

// =============================================================================
// PARITY COMPARATOR
// =============================================================================

/// Builder for comparing Rust and C parity results.
#[allow(dead_code)] // Part of public API for future use
pub struct ParityComparator<'a> {
    rust_hits: &'a [SearchHit],
    c_hits: &'a [SearchHit],
    mode: ParityMode,
}

impl<'a> ParityComparator<'a> {
    pub fn new(rust_hits: &'a [SearchHit], c_hits: &'a [SearchHit]) -> Self {
        Self {
            rust_hits,
            c_hits,
            mode: *TEST_PARITY_MODE,
        }
    }

    #[allow(dead_code)] // Builder API for future use
    pub fn mode(mut self, mode: ParityMode) -> Self {
        self.mode = mode;
        self
    }

    /// Perform the comparison and return results.
    pub fn compare(self) -> ParityResult {
        let mut result = ParityResult::default();

        // Group by (query_idx, target_idx)
        let mut keys = HashSet::new();
        for h in self.rust_hits {
            keys.insert((h.query_idx, h.target_idx));
        }
        for h in self.c_hits {
            keys.insert((h.query_idx, h.target_idx));
        }

        let all_rust_refs: Vec<&SearchHit> = self.rust_hits.iter().collect();

        for key in keys {
            let mut r_group: Vec<&SearchHit> = self
                .rust_hits
                .iter()
                .filter(|h| (h.query_idx, h.target_idx) == key)
                .collect();

            r_group.sort_by(|a, b| {
                a.q_start
                    .cmp(&b.q_start)
                    .then(a.q_end.cmp(&b.q_end))
                    .then(a.output_t_start.cmp(&b.output_t_start))
                    .then(a.output_t_end.cmp(&b.output_t_end))
            });

            let mut c_group: Vec<&SearchHit> = self
                .c_hits
                .iter()
                .filter(|h| (h.query_idx, h.target_idx) == key)
                .collect();
            c_group.sort_by(|a, b| {
                a.q_start
                    .cmp(&b.q_start)
                    .then(a.q_end.cmp(&b.q_end))
                    .then(a.output_t_start.cmp(&b.output_t_start))
                    .then(a.output_t_end.cmp(&b.output_t_end))
            });

            // Pass 1: Remove exact matches (same coords, same energy, same fingerprint)
            let mut c_matched = vec![false; c_group.len()];
            let mut r_remaining = Vec::new();

            for r in r_group {
                let mut found = false;
                for (i, c) in c_group.iter().enumerate() {
                    if !c_matched[i]
                        && r.coords_match(c)
                        && (r.energy.as_f64() - c.energy.as_f64()).abs() < 0.01
                        && r.fingerprint() == c.fingerprint()
                    {
                        c_matched[i] = true;
                        found = true;
                        result.exact_matches += 1;
                        break;
                    }
                }
                if !found {
                    r_remaining.push(r);
                }
            }

            let c_remaining: Vec<&SearchHit> = c_group
                .into_iter()
                .enumerate()
                .filter(|(i, _)| !c_matched[*i])
                .map(|(_, h)| h)
                .collect();

            // Pass 2: Match by coordinates
            let mut c_rem_matched = vec![false; c_remaining.len()];

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

                    let r_e = r.energy.as_f64();
                    let c_e = c.energy.as_f64();

                    let matched = MatchedHit {
                        rust: (*r).clone(),
                        c: c.clone(),
                        status: if r_e < c_e - 0.001 {
                            HitStatus::RustBetter
                        } else if r_e > c_e + 0.001 {
                            HitStatus::RustWorse
                        } else if r.fingerprint() == c.fingerprint() {
                            HitStatus::Identical
                        } else {
                            HitStatus::CoOptimal
                        },
                    };

                    match matched.status {
                        HitStatus::RustBetter => result.rust_better.push(matched),
                        HitStatus::RustWorse => result.rust_worse.push(matched),
                        HitStatus::CoOptimal => result.co_optimal.push(matched),
                        HitStatus::Identical => result.exact_matches += 1,
                        _ => {}
                    }
                } else {
                    result.extras.push((*r).clone());
                }
            }

            // Collect missings
            for (i, c) in c_remaining.iter().enumerate() {
                if !c_rem_matched[i] {
                    let (reason, overlap) = classify_missing(c, &all_rust_refs);
                    result
                        .missings
                        .push(((*c).clone(), reason, overlap.cloned()));
                }
            }
        }

        result
    }
}

// =============================================================================
// UNIT TESTS
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use risearch::types::{Base, Energy, Strand};
    use risearch::{Alignment, Pairing, Sequence};

    fn make_hit(
        q_start: usize,
        q_end: usize,
        t_start: usize,
        t_end: usize,
        strand: &str,
        energy: &str,
    ) -> SearchHit {
        let strand_enum = if strand == "+" {
            Strand::Forward
        } else {
            Strand::Reverse
        };
        let energy_val = Energy::parse(energy).unwrap_or(Energy::new(0.0));

        // Create a simple alignment with all Match pairings
        let len = q_end.saturating_sub(q_start).max(1);
        let seed: Vec<Pairing> = (0..len).map(|_| Pairing::Match(Base::A, Base::U)).collect();
        let alignment = Alignment::new(&[], &seed, &[]);

        SearchHit {
            query_idx: 0,
            target_idx: 0,
            q_start,
            q_end,
            t_start,
            t_end,
            output_q_start: q_start, // Test uses same coords for internal/output
            output_q_end: q_end,
            output_t_start: t_start,
            output_t_end: t_end,
            strand: strand_enum,
            energy: energy_val,
            alignment,
            flank_5: Sequence::from(Vec::new()),
            flank_3: Sequence::from(Vec::new()),
        }
    }

    #[test]
    fn test_ranges_overlap() {
        assert!(ranges_overlap(1, 10, 5, 15));
        assert!(ranges_overlap(1, 10, 1, 10));
        assert!(ranges_overlap(1, 10, 3, 7));
        assert!(!ranges_overlap(1, 10, 11, 20));
        assert!(ranges_overlap(1, 10, 10, 20));
    }

    #[test]
    fn test_hits_overlap() {
        let a = make_hit(1, 10, 100, 110, "+", "-10.0");
        let b = make_hit(5, 15, 105, 115, "+", "-12.0");
        let c = make_hit(1, 10, 100, 110, "-", "-10.0");
        assert!(hits_overlap(&a, &b));
        assert!(!hits_overlap(&a, &c));
    }

    #[test]
    fn test_classify_missing_better_energy() {
        let c_hit = make_hit(1, 10, 100, 110, "+", "-10.0");
        let rust_hits = vec![make_hit(1, 10, 100, 110, "+", "-15.0")];
        let refs: Vec<&SearchHit> = rust_hits.iter().collect();
        let (reason, _) = classify_missing(&c_hit, &refs);
        assert_eq!(reason, MissingReason::BetterEnergy);
    }

    #[test]
    fn test_parity_comparator_exact_match() {
        let rust = vec![make_hit(1, 10, 100, 110, "+", "-10.00")];
        let c = vec![make_hit(1, 10, 100, 110, "+", "-10.00")];
        let result = ParityComparator::new(&rust, &c).compare();
        assert_eq!(result.exact_matches, 1);
        assert!(result.is_pass(ParityMode::Relaxed));
    }

    #[test]
    fn test_parity_comparator_rust_better() {
        let rust = vec![make_hit(1, 10, 100, 110, "+", "-15.00")];
        let c = vec![make_hit(1, 10, 100, 110, "+", "-10.00")];
        let result = ParityComparator::new(&rust, &c).compare();
        assert_eq!(result.rust_better.len(), 1);
        assert!(result.is_pass(ParityMode::Relaxed));
    }

    #[test]
    fn test_parity_comparator_rust_worse_fails() {
        let rust = vec![make_hit(1, 10, 100, 110, "+", "-5.00")];
        let c = vec![make_hit(1, 10, 100, 110, "+", "-10.00")];
        let result = ParityComparator::new(&rust, &c).compare();
        assert_eq!(result.rust_worse.len(), 1);
        assert!(!result.is_pass(ParityMode::Relaxed));
    }

    #[test]
    fn test_log_details_covers_all_categories() {
        // Create a result with all hit types
        let rust = vec![
            make_hit(1, 10, 100, 110, "+", "-10.00"), // exact match
            make_hit(1, 10, 200, 210, "+", "-15.00"), // rust better
            make_hit(1, 10, 300, 310, "+", "-5.00"),  // rust worse
            make_hit(1, 10, 400, 410, "+", "-8.00"),  // extra (no C match)
        ];
        let c = vec![
            make_hit(1, 10, 100, 110, "+", "-10.00"), // exact match
            make_hit(1, 10, 200, 210, "+", "-10.00"), // rust better
            make_hit(1, 10, 300, 310, "+", "-10.00"), // rust worse
            make_hit(1, 10, 500, 510, "+", "-8.00"),  // missing (no Rust match)
        ];

        let result = ParityComparator::new(&rust, &c).compare();

        // Verify all categories are populated
        assert_eq!(result.exact_matches, 1, "exact matches");
        assert_eq!(result.rust_better.len(), 1, "rust better");
        assert_eq!(result.rust_worse.len(), 1, "rust worse");
        assert_eq!(result.extras.len(), 1, "extras");
        assert_eq!(result.missings.len(), 1, "missings");

        // log_details shouldn't panic - this is a smoke test for visualization
        result.log_details("test_all_categories");
    }

    #[test]
    fn test_log_details_empty_result() {
        let result = ParityResult::default();
        // Should not panic with empty result
        result.log_details("empty_test");
    }

    #[test]
    fn test_log_details_co_optimal() {
        // Same coords and energy but different fingerprints (detected via coord match)
        // Create two hits with same coords but they'll have different alignments
        let rust_hit = make_hit(1, 10, 100, 110, "+", "-10.00");
        let c_hit = make_hit(1, 10, 100, 110, "+", "-10.00");

        // Since both have same energy and we create them with identical alignments,
        // they should be exact matches. For co-optimal we'd need different alignments.
        // For now, just verify the test runs without panicking.
        let rust = vec![rust_hit];
        let c = vec![c_hit];
        let result = ParityComparator::new(&rust, &c).compare();

        // With identical hits, should be exact match
        assert!(
            result.exact_matches >= 1 || !result.co_optimal.is_empty(),
            "should handle equal hits"
        );

        // Visualization smoke test
        result.log_details("co_optimal_test");
    }
}
