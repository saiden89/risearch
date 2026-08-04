//! Test-only extensions for `SearchHit`.
//!
//! This module provides:
//! - `SearchHitExt` trait: Comparison helpers for parity testing
//! - `parse_bindingsite_output()`: Parse bindingsite output into `SearchHit`

use anyhow::{anyhow, bail, Context, Result};
use risearch::alignment::{fingerprint_symbols, AlignColumn};
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
            query_registry.get_name(usize::try_from(self.query_idx).unwrap()),
            target_registry.get_name(usize::try_from(self.target_idx).unwrap())
        )
    }

    fn fingerprint(&self) -> Option<String> {
        self.alignment
            .as_ref()
            .map(|columns| fingerprint_symbols(columns).collect())
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
    let interaction = strip_c_markers(fields[8]);

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
        query_idx: u32::try_from(query_idx).context("query index does not fit u32")?,
        target_idx: u32::try_from(target_idx).context("target index does not fit u32")?,
        q_start: q_start.saturating_sub(1),
        q_end: q_end.saturating_sub(1),
        t_start: t_start.saturating_sub(1),
        t_end: t_end.saturating_sub(1),
        strand,
        energy,
        alignment: None,
    };

    // Resolve C's reported classes against our own registries, reusing the hit's
    // slice accessors rather than restating the reverse-strand remap. Keep the
    // reported class even if it differs from what Rust would derive.
    let q_seq = query_registry.entries()[query_idx].sequence();
    let mut query = hit.query(q_seq).iter().copied();
    let mut target = hit.target(target_registry.view()).iter().copied();
    let mut columns = Vec::with_capacity(classes.len());
    for class in classes {
        let column = match class {
            PairClass::TargetBulge => AlignColumn::target_only(
                target
                    .next()
                    .context("C alignment consumes too many target bases")?,
            ),
            PairClass::QueryBulge => AlignColumn::query_only(
                query
                    .next()
                    .context("C alignment consumes too many query bases")?,
            ),
            _ => AlignColumn::reported_pair(
                class,
                query
                    .next()
                    .context("C alignment consumes too many query bases")?,
                target
                    .next()
                    .context("C alignment consumes too many target bases")?,
            ),
        };
        columns.push(column);
    }
    hit.alignment = Some(columns.into_boxed_slice());

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

/// Strip the seed markers embedded by the C output format.
fn strip_c_markers(s: &str) -> String {
    s.chars().filter(|&c| c != 'y' && c != 'x').collect()
}
