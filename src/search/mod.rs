use crate::config::SearchArgs;
use crate::dp;
use crate::dsm::EnergyModel;
use crate::seed::SeedHit;
use crate::seed::SeedMatch;
use crate::seq::Sequence;
use crate::registry::QueryRegistry;
use crate::types::{Alignment, Energy, Strand};
use crate::sa::TargetRegistry;

use std::collections::HashMap;

mod core;
mod stream;

pub use core::run_search;
pub use stream::run_search_streaming;

const MAX_DP_EXT: usize = 50;

/// High-level algorithm stages for structured logging
#[derive(Debug, Clone, Copy)]
enum SearchStage {
    Extend, // DP extension, maximality checks
}

impl std::fmt::Display for SearchStage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Extend => write!(f, "[EXTEND]"),
        }
    }
}

/// Reasons why a seed, candidate, or hit was filtered out
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FilterReason {
    // Seed generation
    SeedContainsN,

    // Candidate validation
    SeedOutOfBounds,

    // Maximality check
    MaximalityLeft,
    MaximalityRight,

    // Energy filtering
    EnergyAboveThreshold,
}

/// Statistics for search filtering
#[derive(Debug, Default)]
pub struct SearchStats {
    pub seeds_tried: usize,
    pub candidates_processed: usize,
    pub hits_final: usize,
    pub filtered: HashMap<FilterReason, usize>,
}

impl SearchStats {
    pub fn record_filter(&mut self, reason: FilterReason) {
        *self.filtered.entry(reason).or_insert(0) += 1;
    }
}

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SearchHit {
    pub query_idx: u32,
    pub target_idx: u32,

    pub q_start: usize,        // 0-based internal
    pub q_end: usize,          // 0-based internal
    pub t_start: usize,        // 0-based internal
    pub t_end: usize,          // 0-based internal
    pub output_q_start: usize, // 1-based for output/comparison
    pub output_q_end: usize,   // 1-based for output/comparison
    pub output_t_start: usize, // 1-based, strand-aware
    pub output_t_end: usize,   // 1-based, strand-aware
    pub strand: Strand,
    pub energy: Energy,
    pub alignment: Alignment,
    pub flank_5: Sequence,
    pub flank_3: Sequence,
}

impl SearchHit {
    /// Parse a SearchHit from C risearch output line.
    ///
    /// C output format (tab-separated):
    /// `q_id, q_start, q_end, t_id, t_start, t_end, strand, energy, interaction, target_seq, [flank_5, flank_3]`
    ///
    /// Handles C quirks:
    /// - Seed markers 'y' and 'x' in interaction/target strings
    /// - Missing optional columns (flanks)
    /// - 1-based coordinates
    pub fn from_c_output(
        line: &str,
        query_registry: &QueryRegistry,
        target_registry: &TargetRegistry,
    ) -> Option<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            return None;
        }
        let query_idx = query_registry.index_of(fields[0])?;
        let target_idx = target_registry.index_of(fields[3])?;

        // Parse and strip seed markers from interaction
        let (interaction, seed_start, seed_end) = Self::strip_c_markers(fields[8]);
        let target_seq = Self::strip_markers_simple(fields[9]);

        // Parse coordinates (C uses 1-based)
        let q_start: usize = fields[1].parse().ok()?;
        let q_end: usize = fields[2].parse().ok()?;
        let t_start: usize = fields[4].parse().ok()?;
        let t_end: usize = fields[5].parse().ok()?;

        // Parse strand
        let strand = fields[6].chars().next().unwrap_or('+').into();

        // Parse energy
        let energy = Energy::parse(fields[7])?;

        // Create alignment from fingerprint and target sequence
        let alignment = Alignment::from_c_output(&interaction, &target_seq, seed_start, seed_end);

        // Optional flanks
        let flank_5 = match fields.get(10) {
            Some(s) => {
                let clean = Self::strip_markers_simple(s);
                Sequence::normalize("flank_5", clean.as_bytes()).ok()?.0
            }
            None => Sequence::from(Vec::new()),
        };
        let flank_3 = match fields.get(11) {
            Some(s) => {
                let clean = Self::strip_markers_simple(s);
                Sequence::normalize("flank_3", clean.as_bytes()).ok()?.0
            }
            None => Sequence::from(Vec::new()),
        };

        Some(SearchHit {
            query_idx,
            target_idx,
            q_start: q_start.saturating_sub(1), // 0-based internal
            q_end: q_end.saturating_sub(1),     // 0-based internal
            t_start: t_start.saturating_sub(1), // 0-based internal
            t_end: t_end.saturating_sub(1),     // 0-based internal
            output_q_start: q_start,            // 1-based for comparison
            output_q_end: q_end,                // 1-based for comparison
            output_t_start: t_start,            // 1-based
            output_t_end: t_end,                // 1-based
            strand,
            energy,
            alignment,
            flank_5,
            flank_3,
        })
    }

    /// Strip y/x seed markers and extract seed range positions.
    fn strip_c_markers(s: &str) -> (String, Option<usize>, Option<usize>) {
        let y_pos = s.find('y');
        let x_pos = s.find('x');

        let (seed_start, seed_end) = match (y_pos, x_pos) {
            (Some(y), Some(x)) if x > y => (Some(y), Some(x - 1)),
            _ => (None, None),
        };

        let clean: String = s.chars().filter(|&c| c != 'y' && c != 'x').collect();
        (clean, seed_start, seed_end)
    }

    /// Strip y/x markers without tracking positions.
    fn strip_markers_simple(s: &str) -> String {
        s.chars().filter(|&c| c != 'y' && c != 'x').collect()
    }

    // =========================================================================
    // COMPARISON HELPERS (for parity testing)
    // =========================================================================

    /// Returns true if coordinates and strand match another hit.
    /// Uses output coordinates (1-based) for comparison.
    pub fn coords_match(&self, other: &Self) -> bool {
        self.output_q_start == other.output_q_start
            && self.output_q_end == other.output_q_end
            && self.output_t_start == other.output_t_start
            && self.output_t_end == other.output_t_end
            && self.strand == other.strand
    }

    /// Group key for matching hits (query_idx:target_idx -> names).
    pub fn group_key(&self, query_registry: &QueryRegistry, target_registry: &TargetRegistry) -> String {
        format!(
            "{}:{}",
            query_registry.get_name(self.query_idx),
            target_registry.get_name(self.target_idx)
        )
    }

    /// Fingerprint string for comparison.
    pub fn fingerprint(&self) -> String {
        self.alignment.fingerprint()
    }

    /// Target sequence for comparison.
    pub fn target_seq(&self) -> String {
        self.alignment.target_sequence()
    }

    /// Query sequence for comparison (may have N placeholders if from C output).
    pub fn query_seq(&self) -> String {
        self.alignment.query_sequence()
    }

    /// Seed start position within interaction.
    pub fn seed_start(&self) -> Option<usize> {
        Some(self.alignment.left_extension().len())
    }

    /// Seed end position within interaction.
    pub fn seed_end(&self) -> Option<usize> {
        let start = self.alignment.left_extension().len();
        let seed_len = self.alignment.seed().len();
        Some(start + seed_len)
    }

    /// Format for debug output (uses 1-based output coordinates).
    pub fn fmt_coords(&self) -> String {
        format!(
            "q=[{},{}] t=[{},{}] S={} E={}",
            self.output_q_start,
            self.output_q_end,
            self.output_t_start,
            self.output_t_end,
            self.strand,
            self.energy
        )
    }
}

// Reimplementing mapping locally for safety and speed

pub struct SearchContext<'a, 'e> {
    pub index: &'a TargetRegistry,
    pub args: &'a SearchArgs,
    pub extender: &'e mut dp::DpExtender,
    pub stats: SearchStats,
    pub energy: EnergyModel,
}

impl<'a, 'e> SearchContext<'a, 'e> {
    pub fn with_extender(
        index: &'a TargetRegistry,
        args: &'a SearchArgs,
        extender: &'e mut dp::DpExtender,
    ) -> Self {
        Self {
            index,
            args,
            extender,
            stats: SearchStats::default(),
            energy: EnergyModel::T04,
        }
    }
}

// Thread-local storage for reusable buffers - avoids allocation per query.
thread_local! {
    static THREAD_EXTENDER: std::cell::RefCell<dp::DpExtender> =
        std::cell::RefCell::new(dp::DpExtender::new());
    static THREAD_SEEDS: std::cell::RefCell<Vec<SeedHit>> =
        std::cell::RefCell::new(Vec::with_capacity(100_000));
    static THREAD_MATCHES: std::cell::RefCell<Vec<SeedMatch>> =
        std::cell::RefCell::new(Vec::with_capacity(1024));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_search_hit_from_c_output() {
        // Sample C output line (C uses 1-based coords)
        let line =
            "hsa-miR-1\t1\t10\tENSG00000001\t100\t110\t+\t-15.50\tyPPPUPPPx\tacguacgu\tAA\tCC";

        let query_registry = QueryRegistry::from_names(vec!["hsa-miR-1".to_string()]);
        let target_registry = TargetRegistry::from_names(vec!["ENSG00000001".to_string()]);
        let hit =
            SearchHit::from_c_output(line, &query_registry, &target_registry).expect("should parse");

        assert_eq!(hit.query_idx, 0);
        assert_eq!(hit.target_idx, 0);
        // Internal coords are 0-based (converted from C's 1-based)
        assert_eq!(hit.q_start, 0);
        assert_eq!(hit.q_end, 9);
        assert_eq!(hit.t_start, 99);
        assert_eq!(hit.t_end, 109);
        // Output coords stay 1-based for display
        assert_eq!(hit.output_q_start, 1);
        assert_eq!(hit.output_q_end, 10);
        assert_eq!(hit.output_t_start, 100);
        assert_eq!(hit.output_t_end, 110);
        assert_eq!(hit.strand, Strand::Forward);
        assert!((hit.energy.as_f64() - (-15.50)).abs() < 0.01);
        assert_eq!(hit.alignment.fingerprint(), "PPPUPPP"); // markers stripped
        let (f5, _) = Sequence::normalize("f5", b"AA").unwrap();
        let (f3, _) = Sequence::normalize("f3", b"CC").unwrap();
        assert_eq!(hit.flank_5, f5);
        assert_eq!(hit.flank_3, f3);
    }

    #[test]
    fn test_search_hit_from_c_output_minimal() {
        // Minimal 10 columns (no flanks), C uses 1-based coords
        let line = "q1\t1\t5\tt1\t10\t15\t-\t-8.00\tPPPPP\tacgua";

        let query_registry = QueryRegistry::from_names(vec!["q1".to_string()]);
        let target_registry = TargetRegistry::from_names(vec!["t1".to_string()]);
        let hit =
            SearchHit::from_c_output(line, &query_registry, &target_registry).expect("should parse minimal");
        assert_eq!(hit.query_idx, 0);
        // Internal coords are 0-based
        assert_eq!(hit.q_start, 0);
        assert_eq!(hit.q_end, 4);
        assert_eq!(hit.t_start, 9);
        assert_eq!(hit.t_end, 14);
        // Output coords stay 1-based
        assert_eq!(hit.output_q_start, 1);
        assert_eq!(hit.output_q_end, 5);
        assert_eq!(hit.output_t_start, 10);
        assert_eq!(hit.output_t_end, 15);
        assert_eq!(hit.strand, Strand::Reverse);
        assert!(hit.flank_5.is_empty());
        assert!(hit.flank_3.is_empty());
    }

    #[test]
    fn test_search_hit_from_c_output_invalid() {
        let query_registry = QueryRegistry::from_names(vec![]);
        let target_registry = TargetRegistry::from_names(vec![]);
        assert!(SearchHit::from_c_output("too\tfew\tcolumns", &query_registry, &target_registry).is_none());
        assert!(SearchHit::from_c_output("", &query_registry, &target_registry).is_none());
    }
}
