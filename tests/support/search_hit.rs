//! Test-only extensions for `SearchHit`.
//!
//! This module provides:
//! - `SearchHitExt` trait: Comparison helpers for parity testing
//! - `parse_bindingsite_output()`: Parse bindingsite output into `SearchHit`

use risearch::alignment::Alignment;
use risearch::index::store::TargetRegistry;
use risearch::registry::QueryRegistry;
use risearch::types::{Energy, Strand};
use risearch::SearchHit;

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

    /// Seed start position within interaction.
    fn seed_start(&self) -> Option<usize>;

    /// Seed end position within interaction.
    fn seed_end(&self) -> Option<usize>;

    /// Format for debug output (uses 1-based output coordinates).
    fn fmt_coords(&self) -> String;
}

impl SearchHitExt for SearchHit {
    fn coords_match(&self, other: &Self) -> bool {
        (self.q_start + 1) == (other.q_start + 1)
            && (self.q_end + 1) == (other.q_end + 1)
            && (self.t_start + 1) == (other.t_start + 1)
            && (self.t_end + 1) == (other.t_end + 1)
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

    fn seed_start(&self) -> Option<usize> {
        self.seed_start
    }

    fn seed_end(&self) -> Option<usize> {
        self.seed_end
    }

    fn fmt_coords(&self) -> String {
        format!(
            "q=[{},{}] t=[{},{}] S={} E={:.2}",
            self.q_start + 1,
            self.q_end + 1,
            self.t_start + 1,
            self.t_end + 1,
            self.strand,
            self.energy
        )
    }
}

// =============================================================================
// BINDINGSITE OUTPUT PARSING
// =============================================================================

/// Parse a SearchHit from Rust or C bindingsite output line.
///
/// Binding-site output format (tab-separated):
/// `q_id, q_start, q_end, t_id, t_start, t_end, strand, energy, interaction, target_seq, [flank_5, flank_3]`
///
/// Handles C quirks when present:
/// - Seed markers 'y' and 'x' in interaction/target strings
/// - Missing optional columns (flanks)
/// - 1-based coordinates
pub(crate) fn parse_bindingsite_output(
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
    let energy = Energy::from_kcal(fields[7].parse::<f64>().ok()?);

    // Create alignment from C interaction/target columns.
    let alignment = Alignment::from_c_output(&interaction, &target_seq);
    let (seed_start, seed_end) = match (seed_start, seed_end) {
        (Some(s), Some(e)) => {
            let s = s.min(interaction.len());
            let e = e.min(interaction.len()).max(s);
            (Some(s), Some(e))
        }
        _ => (None, None),
    };

    Some(SearchHit {
        query_idx,
        target_idx,
        q_start: q_start.saturating_sub(1),
        q_end: q_end.saturating_sub(1),
        t_start: t_start.saturating_sub(1),
        t_end: t_end.saturating_sub(1),
        strand,
        energy,
        seed_start,
        seed_end,
        alignment: Some(alignment),
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
