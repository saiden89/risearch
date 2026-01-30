//! Test-only extensions for `SearchHit`.
//!
//! This module provides:
//! - `SearchHitExt` trait: Comparison helpers for parity testing
//! - `parse_c_output()`: Parse C risearch output into `SearchHit`

use risearch::SearchHit;
use risearch::alignment::Alignment;
use risearch::registry::{QueryRegistry, TargetRegistry};
use risearch::seq::Sequence;
use risearch::types::{Energy, Strand};

// =============================================================================
// EXTENSION TRAIT: Parity-testing helpers
// =============================================================================

/// Extension trait providing parity-testing methods for `SearchHit`.
pub(crate) trait SearchHitExt {
    /// Returns true if coordinates and strand match another hit.
    fn coords_match(&self, other: &Self) -> bool;

    /// Group key for matching hits (query_idx:target_idx -> names).
    #[allow(dead_code)]
    fn group_key(&self, query_registry: &QueryRegistry, target_registry: &TargetRegistry)
    -> String;

    /// Fingerprint string for comparison. None if no alignment data.
    fn fingerprint(&self) -> Option<String>;

    /// Target sequence for comparison. None if no alignment data.
    fn target_seq(&self) -> Option<String>;

    /// Query sequence for comparison (may have N placeholders if from C output).
    fn query_seq(&self) -> Option<String>;

    /// Seed start position within interaction.
    fn seed_start(&self) -> Option<usize>;

    /// Seed end position within interaction.
    fn seed_end(&self) -> Option<usize>;

    /// Format for debug output (uses 1-based output coordinates).
    fn fmt_coords(&self) -> String;
}

impl SearchHitExt for SearchHit {
    fn coords_match(&self, other: &Self) -> bool {
        self.output_q_start == other.output_q_start
            && self.output_q_end == other.output_q_end
            && self.output_t_start == other.output_t_start
            && self.output_t_end == other.output_t_end
            && self.strand == other.strand
    }

    fn group_key(
        &self,
        query_registry: &QueryRegistry,
        target_registry: &TargetRegistry,
    ) -> String {
        format!(
            "{}:{}",
            query_registry.get_name(self.query_idx),
            target_registry.get_name(self.target_idx)
        )
    }

    fn fingerprint(&self) -> Option<String> {
        self.alignment.as_ref().map(|a| a.fingerprint())
    }

    fn target_seq(&self) -> Option<String> {
        self.alignment.as_ref().map(|a| a.target_sequence())
    }

    fn query_seq(&self) -> Option<String> {
        self.alignment.as_ref().map(|a| a.query_sequence())
    }

    fn seed_start(&self) -> Option<usize> {
        self.alignment.as_ref().map(|a| a.left_extension().len())
    }

    fn seed_end(&self) -> Option<usize> {
        self.alignment.as_ref().map(|a| {
            let start = a.left_extension().len();
            let seed_len = a.seed().len();
            start + seed_len
        })
    }

    fn fmt_coords(&self) -> String {
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

// =============================================================================
// C OUTPUT PARSING
// =============================================================================

/// Parse a SearchHit from C risearch output line.
///
/// C output format (tab-separated):
/// `q_id, q_start, q_end, t_id, t_start, t_end, strand, energy, interaction, target_seq, [flank_5, flank_3]`
///
/// Handles C quirks:
/// - Seed markers 'y' and 'x' in interaction/target strings
/// - Missing optional columns (flanks)
/// - 1-based coordinates
pub(crate) fn parse_c_output(
    line: &str,
    query_registry: &QueryRegistry,
    target_registry: &TargetRegistry,
) -> Option<SearchHit> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 10 {
        return None;
    }
    let query_idx = query_registry.index_of(fields[0])?;
    let target_idx = target_registry.index_of(fields[3])?;

    // Parse and strip seed markers from interaction
    let (interaction, seed_start, seed_end) = strip_c_markers(fields[8]);
    let target_seq = strip_markers_simple(fields[9]);

    // Parse coordinates (C uses 1-based)
    let q_start: usize = fields[1].parse().ok()?;
    let q_end: usize = fields[2].parse().ok()?;
    let t_start: usize = fields[4].parse().ok()?;
    let t_end: usize = fields[5].parse().ok()?;

    // Parse strand
    let strand: Strand = fields[6].chars().next().unwrap_or('+').into();

    // Parse energy
    let energy = Energy::parse(fields[7])?;

    // Create alignment from fingerprint and target sequence
    let alignment = Alignment::from_c_output(&interaction, &target_seq, seed_start, seed_end);

    // Optional flanks
    let flank_5 = match fields.get(10) {
        Some(s) => {
            let clean = strip_markers_simple(s);
            Sequence::normalize("flank_5", clean.as_bytes()).ok()?.0
        }
        None => Sequence::from(Vec::new()),
    };
    let flank_3 = match fields.get(11) {
        Some(s) => {
            let clean = strip_markers_simple(s);
            Sequence::normalize("flank_3", clean.as_bytes()).ok()?.0
        }
        None => Sequence::from(Vec::new()),
    };

    Some(SearchHit {
        query_idx,
        target_idx,
        q_start: q_start.saturating_sub(1),
        q_end: q_end.saturating_sub(1),
        t_start: t_start.saturating_sub(1),
        t_end: t_end.saturating_sub(1),
        output_q_start: q_start,
        output_q_end: q_end,
        output_t_start: t_start,
        output_t_end: t_end,
        strand,
        energy,
        alignment: Some(alignment),
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
