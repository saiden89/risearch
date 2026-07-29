//! Test-only extensions for `SearchHit`.
//!
//! This module provides:
//! - `SearchHitExt` trait: Comparison helpers for parity testing
//! - `parse_bindingsite_output()`: Parse bindingsite output into `SearchHit`

use anyhow::{anyhow, bail, Context, Result};
use risearch::alignment::Alignment;
use risearch::index::store::TargetRegistry;
use risearch::registry::QueryRegistry;
use risearch::types::{Energy, Strand};
use risearch::{PairClass, SearchHit};

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
        self.alignment
            .as_ref()
            .and_then(|a| a.seed())
            .map(|s| s.start)
    }

    fn seed_end(&self) -> Option<usize> {
        self.alignment
            .as_ref()
            .and_then(|a| a.seed())
            .map(|s| s.end)
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
) -> Result<SearchHit> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 10 {
        bail!("expected at least 10 columns, got {}", fields.len());
    }
    let query_idx = (0..query_registry.len())
        .find(|&i| name_token(query_registry.get_name(i)) == name_token(fields[0]))
        .with_context(|| format!("unknown query name {:?}", fields[0]))?;
    let target_idx = (0..target_registry.len())
        .find(|&i| name_token(target_registry.get_name(i)) == name_token(fields[3]))
        .with_context(|| format!("unknown target name {:?}", fields[3]))?;

    // Parse and strip seed markers from interaction
    let (interaction, seed_start, seed_end) = strip_c_markers(fields[8]);

    // Parse coordinates (C uses 1-based)
    let q_start: usize = column(&fields, 1, "q_start")?;
    let q_end: usize = column(&fields, 2, "q_end")?;
    let t_start: usize = column(&fields, 4, "t_start")?;
    let t_end: usize = column(&fields, 5, "t_end")?;

    // Parse strand
    let strand = Strand::try_from(fields[6].chars().next().context("empty strand column")?)
        .map_err(|e| anyhow!(e))?;

    // Parse energy
    let energy = Energy::from_kcal(column(&fields, 7, "energy")?);

    let classes: Vec<PairClass> = interaction
        .chars()
        .map(|c| match c {
            'P' => PairClass::Canonical,
            'W' => PairClass::Wobble,
            'U' => PairClass::Mismatch,
            'T' => PairClass::TargetBulge,
            'Q' => PairClass::QueryBulge,
            _ => PairClass::Mismatch,
        })
        .collect();

    let mut hit = SearchHit {
        query_idx,
        target_idx,
        q_start: q_start.saturating_sub(1),
        q_end: q_end.saturating_sub(1),
        t_start: t_start.saturating_sub(1),
        t_end: t_end.saturating_sub(1),
        strand,
        energy,
        alignment: None,
    };

    // Resolve C's classes into columns against our own registries, reusing the
    // hit's slice accessors rather than restating the reverse-strand remap.
    let seed = match (seed_start, seed_end) {
        (Some(s), Some(e)) => {
            let s = s.min(classes.len());
            Some(s..e.min(classes.len()).max(s))
        }
        _ => None,
    };
    let q_seq = query_registry.entries()[query_idx].sequence();
    hit.alignment = Some(Box::new(Alignment::from_classes(
        &classes,
        seed,
        hit.query(q_seq),
        hit.target(target_registry.view()),
    )));

    Ok(hit)
}

/// C truncates FASTA ids at the first whitespace (`fasta.c` strtok); Rust keeps the
/// whole header, so the two sides' name columns only agree on the first token.
fn name_token(name: &str) -> &str {
    name.split_ascii_whitespace().next().unwrap_or("")
}

fn column<T: std::str::FromStr>(fields: &[&str], idx: usize, name: &str) -> Result<T>
where
    T::Err: std::fmt::Display,
{
    fields[idx]
        .parse()
        .map_err(|e| anyhow!("invalid {name} {:?}: {e}", fields[idx]))
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
