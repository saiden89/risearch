//! Comparison logic for parity testing.
//!
//! Contains `ParityComparator` - a builder pattern for comparing Rust and C results.

use log::{debug, info};
use std::collections::HashSet;

use crate::common::record::Rec;
use crate::common::status::{HitStatus, MissingReason, ParityMode};
use crate::common::table::{ParityKind, ParityTable, TableConfig};

// =============================================================================
// LOG TAG
// =============================================================================

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
#[derive(Debug)]
#[allow(dead_code)] // Part of public API for future use
pub struct MatchedHit {
    pub rust: Rec,
    pub c: Rec,
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
    pub extras: Vec<Rec>,
    pub missings: Vec<(Rec, MissingReason)>,
}

impl ParityResult {
    /// Check if the result is a pass given the parity mode.
    pub fn is_pass(&self, mode: ParityMode) -> bool {
        if !self.rust_worse.is_empty() {
            return false;
        }
        match mode {
            ParityMode::Strict => {
                self.extras.is_empty()
                    && self.missings.is_empty()
                    && self.rust_better.is_empty()
                    && self.co_optimal.is_empty()
            }
            ParityMode::Relaxed => {
                // In relaxed mode, only unacceptable missings count as failures
                !self.missings.iter().any(|(_, r)| !r.is_acceptable(mode))
            }
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
pub fn hits_overlap(a: &Rec, b: &Rec) -> bool {
    a.strand == b.strand
        && ranges_overlap(a.t_start, a.t_end, b.t_start, b.t_end)
        && ranges_overlap(a.q_start, a.q_end, b.q_start, b.q_end)
}

/// Classify why a C hit is missing from Rust output.
/// Returns the reason and optionally the best overlapping Rust hit.
pub fn classify_missing<'a>(
    c_hit: &Rec,
    rust_hits: &[&'a Rec],
) -> (MissingReason, Option<&'a Rec>) {
    let c_e = c_hit.energy.parse::<f64>().unwrap_or(0.0);

    let best_overlap = rust_hits
        .iter()
        .filter(|r| hits_overlap(c_hit, r))
        .min_by(|a, b| {
            let a_e = a.energy.parse::<f64>().unwrap_or(0.0);
            let b_e = b.energy.parse::<f64>().unwrap_or(0.0);
            a_e.partial_cmp(&b_e).unwrap_or(std::cmp::Ordering::Equal)
        })
        .copied();

    match best_overlap {
        Some(r) => {
            let r_e = r.energy.parse::<f64>().unwrap_or(0.0);
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
    rust_recs: &'a [Rec],
    c_recs: &'a [Rec],
    mode: ParityMode,
}

impl<'a> ParityComparator<'a> {
    pub fn new(rust_recs: &'a [Rec], c_recs: &'a [Rec]) -> Self {
        Self {
            rust_recs,
            c_recs,
            mode: ParityMode::default(),
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

        // Group by (q_id, t_id)
        let mut keys = HashSet::new();
        for r in self.rust_recs {
            keys.insert((r.q_id.clone(), r.t_id.clone()));
        }
        for r in self.c_recs {
            keys.insert((r.q_id.clone(), r.t_id.clone()));
        }

        let all_rust_refs: Vec<&Rec> = self.rust_recs.iter().collect();

        for (q, t) in keys {
            let mut r_group: Vec<&Rec> = self
                .rust_recs
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

            let mut c_group: Vec<&Rec> = self
                .c_recs
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

            // Pass 1: Remove exact matches
            let mut c_matched = vec![false; c_group.len()];
            let mut r_remaining = Vec::new();

            for r in r_group {
                let mut found = false;
                for (i, c) in c_group.iter().enumerate() {
                    if !c_matched[i] && r == *c {
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

            let c_remaining: Vec<&Rec> = c_group
                .into_iter()
                .enumerate()
                .filter(|(i, _)| !c_matched[*i])
                .map(|(_, r)| r)
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

                    let r_e = r.energy.parse::<f64>().unwrap_or(0.0);
                    let c_e = c.energy.parse::<f64>().unwrap_or(0.0);

                    let matched = MatchedHit {
                        rust: (*r).clone(),
                        c: c.clone(),
                        status: if r_e < c_e - 0.001 {
                            HitStatus::RustBetter
                        } else if r_e > c_e + 0.001 {
                            HitStatus::RustWorse
                        } else if r.interaction == c.interaction {
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
                    let (reason, _) = classify_missing(c, &all_rust_refs);
                    result.missings.push(((*c).clone(), reason));
                }
            }
        }

        result
    }
}

// =============================================================================
// ANALYZE HIT PAIRS
// =============================================================================

/// Detailed analysis for focused tests with few hits.
pub fn analyze_hit_pairs(extras: &[&Rec], missings: &[&Rec]) {
    if extras.is_empty() || missings.is_empty() {
        return;
    }

    info!(
        "{} Paired hit analysis ({} extras, {} missings)",
        LogTag::Pair,
        extras.len(),
        missings.len()
    );

    for extra in extras {
        let best_match = missings.iter().min_by_key(|missing| {
            let q_start_diff = (extra.q_start as i32 - missing.q_start as i32).abs();
            let q_end_diff = (extra.q_end as i32 - missing.q_end as i32).abs();
            let t_start_diff = (extra.t_start as i32 - missing.t_start as i32).abs();
            let t_end_diff = (extra.t_end as i32 - missing.t_end as i32).abs();
            q_start_diff + q_end_diff + t_start_diff + t_end_diff
        });

        if let Some(m) = best_match {
            let score = (extra.q_start as i32 - m.q_start as i32).abs()
                + (extra.q_end as i32 - m.q_end as i32).abs()
                + (extra.t_start as i32 - m.t_start as i32).abs()
                + (extra.t_end as i32 - m.t_end as i32).abs();

            debug!("{} Likely pair (distance={})", LogTag::Pair, score);

            if extra.interaction == m.interaction {
                debug!("{} [NOTE] Interactions Identical!", LogTag::Pair);
            }

            let table = ParityTable {
                kind: ParityKind::Mismatch { rust: extra, c: m },
                config: TableConfig::default(),
            };

            for line in table.to_string().lines() {
                debug!("{} {}", LogTag::Pair, line);
            }
        }
    }
}

// =============================================================================
// UNIT TESTS
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_rec(
        q_start: usize,
        q_end: usize,
        t_start: usize,
        t_end: usize,
        strand: &str,
        energy: &str,
    ) -> Rec {
        Rec {
            q_id: "query".into(),
            t_id: "target".into(),
            q_start,
            q_end,
            t_start,
            t_end,
            strand: strand.into(),
            energy: energy.into(),
            interaction: "PPPPP".into(),
            target_seq: "AUGCG".into(),
            query_seq: String::new(),
            flank_5: String::new(),
            flank_3: String::new(),
            seed_start: None,
            seed_end: None,
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
        let a = make_rec(1, 10, 100, 110, "+", "-10.0");
        let b = make_rec(5, 15, 105, 115, "+", "-12.0");
        let c = make_rec(1, 10, 100, 110, "-", "-10.0");
        assert!(hits_overlap(&a, &b));
        assert!(!hits_overlap(&a, &c));
    }

    #[test]
    fn test_classify_missing_better_energy() {
        let c_hit = make_rec(1, 10, 100, 110, "+", "-10.0");
        let rust_hits = vec![make_rec(1, 10, 100, 110, "+", "-15.0")];
        let refs: Vec<&Rec> = rust_hits.iter().collect();
        let (reason, _) = classify_missing(&c_hit, &refs);
        assert_eq!(reason, MissingReason::BetterEnergy);
    }

    #[test]
    fn test_parity_comparator_exact_match() {
        let rust = vec![make_rec(1, 10, 100, 110, "+", "-10.00")];
        let c = vec![make_rec(1, 10, 100, 110, "+", "-10.00")];
        let result = ParityComparator::new(&rust, &c).compare();
        assert_eq!(result.exact_matches, 1);
        assert!(result.is_pass(ParityMode::Relaxed));
    }

    #[test]
    fn test_parity_comparator_rust_better() {
        let rust = vec![make_rec(1, 10, 100, 110, "+", "-15.00")];
        let c = vec![make_rec(1, 10, 100, 110, "+", "-10.00")];
        let result = ParityComparator::new(&rust, &c).compare();
        assert_eq!(result.rust_better.len(), 1);
        assert!(result.is_pass(ParityMode::Relaxed));
    }

    #[test]
    fn test_parity_comparator_rust_worse_fails() {
        let rust = vec![make_rec(1, 10, 100, 110, "+", "-5.00")];
        let c = vec![make_rec(1, 10, 100, 110, "+", "-10.00")];
        let result = ParityComparator::new(&rust, &c).compare();
        assert_eq!(result.rust_worse.len(), 1);
        assert!(!result.is_pass(ParityMode::Relaxed));
    }
}
