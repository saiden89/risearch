use anyhow::{anyhow, bail, Context, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::ops::Range;
use std::path::Path;

use crate::config::SeedConfig;
use crate::fastx::read_fasta_sequences;
use crate::index::io::validate_readable_file;
use crate::index::sa::SuffixArray;
use crate::seq::{SeqView, Sequence};
use crate::types::Base;

pub trait RegistryEntry {
    fn name(&self) -> &str;
}

#[derive(Serialize, Deserialize)]
pub struct Registry<T> {
    entries: Vec<T>,
}

impl<T> Registry<T> {
    pub fn new(entries: Vec<T>) -> Self {
        Self { entries }
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn get(&self, idx: u32) -> &T {
        &self.entries[idx as usize]
    }

    pub fn entries(&self) -> &[T] {
        &self.entries
    }

    pub fn into_entries(self) -> Vec<T> {
        self.entries
    }
}

impl<T: RegistryEntry> Registry<T> {
    /// Get name by index (unchecked for hot path).
    ///
    /// # Safety
    /// Caller must ensure idx is valid (< number of entries).
    /// In practice, idx comes from hit.query_idx which is always valid.
    #[inline(always)]
    pub fn get_name(&self, idx: u32) -> &str {
        debug_assert!(
            (idx as usize) < self.entries.len(),
            "Registry index out of bounds: {} >= {}",
            idx,
            self.entries.len()
        );
        // SAFETY: idx comes from generated hits, guaranteed to be valid
        unsafe { self.entries.get_unchecked(idx as usize).name() }
    }

    pub fn index_of(&self, name: &str) -> Option<u32> {
        self.entries
            .iter()
            .position(|e| e.name() == name)
            .map(|i| i as u32)
    }

    pub fn iter(&self) -> impl Iterator<Item = (u32, &T)> {
        self.entries.iter().enumerate().map(|(i, e)| (i as u32, e))
    }
}

/// Query prepared once at load time for the current seed configuration.
pub struct Query {
    /// Sequence identifier
    name: String,
    /// Full forward sequence used for extension and output
    sequence: Sequence,
    /// Seed-interval slice used by the seeding suffix array
    seed_sequence: Sequence,
    /// Suffix array built over `seed_sequence`
    sa: SuffixArray,
    /// Pre-computed normalized seed interval bounds on the full query
    pub(crate) seed_interval: Range<usize>,
    /// Prefix sum of N positions for O(1) N-checking
    n_prefix: Vec<u32>,
    /// Fast path when query has no Ns
    has_n_any: bool,
}

impl Query {
    fn from_parts(name: String, sequence: Sequence, config: &SeedConfig) -> Result<Self> {
        let q_len = sequence.len();

        // Compute N-prefix for O(1) N-checking
        let mut n_prefix = Vec::with_capacity(q_len + 1);
        n_prefix.push(0);
        let mut n_total = 0;
        for &base in sequence.iter() {
            if base == Base::N {
                n_total += 1;
            }
            n_prefix.push(n_total);
        }
        let has_n_any = n_total != 0;

        // Compute seed interval once and fail early at boundary if invalid.
        let (start1, end1, _min_seed_len) = config
            .seed
            .normalize(q_len)
            .map_err(|err| anyhow!("Invalid seed spec for query '{}': {}", name, err))?;
        let seed_interval = (start1 - 1)..end1;
        let seed_sequence =
            Sequence::from(sequence[seed_interval.start..seed_interval.end].to_vec());
        let sa = SuffixArray::try_from(&seed_sequence[..])
            .map_err(|err| anyhow!("Failed to build seed SA for query '{}': {}", name, err))?;

        Ok(Self {
            name,
            sequence,
            seed_sequence,
            sa,
            seed_interval,
            n_prefix,
            has_n_any,
        })
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[inline(always)]
    pub fn sequence(&self) -> SeqView<'_> {
        self.sequence.as_view()
    }

    #[inline(always)]
    pub fn sa(&self) -> &[u64] {
        self.sa.as_ref()
    }

    #[inline]
    pub fn n_prefix(&self) -> &[u32] {
        &self.n_prefix
    }

    #[inline]
    pub fn has_n_any(&self) -> bool {
        self.has_n_any
    }

    #[inline]
    pub fn seed_interval(&self) -> Range<usize> {
        self.seed_interval.clone()
    }

    #[inline(always)]
    pub fn seed_sequence(&self) -> SeqView<'_> {
        self.seed_sequence.as_view()
    }
}

impl RegistryEntry for Query {
    fn name(&self) -> &str {
        Query::name(self)
    }
}

pub type QueryRegistry = Registry<Query>;

fn read_and_validate_sequences(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    validate_readable_file(path)?;

    let sequences = read_fasta_sequences(path).context("Failed to read FASTA sequences")?;
    if sequences.is_empty() {
        bail!("No sequences found in input file: {}", path.display());
    }

    let mut seen = HashSet::with_capacity(sequences.len());
    for (id, _) in &sequences {
        if id.trim().is_empty() {
            bail!("Encountered empty FASTA record id in {}", path.display());
        }
        if !seen.insert(id.clone()) {
            bail!("Duplicate FASTA record id '{}' in {}", id, path.display());
        }
    }

    Ok(sequences)
}

fn normalize_record(id: String, seq: Vec<u8>) -> Result<Option<(String, Sequence)>> {
    let (seq_norm, stats) = Sequence::normalize(&id, &seq)
        .with_context(|| format!("Failed to normalize sequence '{}'", id))?;

    if seq_norm.is_empty() {
        log::warn!(
            "Skipping empty sequence after normalization: '{}' (removed_gaps={}, converted_to_n={})",
            id,
            stats.removed_gaps,
            stats.converted_to_n
        );
        return Ok(None);
    }

    if stats.removed_gaps > 0 || stats.converted_to_n > 0 {
        log::debug!(
            "Normalized sequence '{}': removed_gaps={}, converted_to_n={}",
            id,
            stats.removed_gaps,
            stats.converted_to_n
        );
    }

    Ok(Some((id, seq_norm)))
}

impl QueryRegistry {
    pub fn from_fasta(path: &Path, config: &SeedConfig) -> Result<Self> {
        let sequences = read_and_validate_sequences(path)?;

        let maybe_entries: Vec<Option<Query>> = sequences
            .into_par_iter()
            .map(|(id, seq)| -> Result<Option<Query>> {
                let Some((name, sequence)) = normalize_record(id, seq)? else {
                    return Ok(None);
                };

                Ok(Some(Query::from_parts(name, sequence, config)?))
            })
            .collect::<Result<Vec<_>>>()?;

        let entries: Vec<Query> = maybe_entries.into_iter().flatten().collect();
        if entries.is_empty() {
            bail!(
                "All sequences were empty after normalization in {}",
                path.display()
            );
        }

        Ok(Self::new(entries))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{MismatchSpec, SeedSpec};

    fn make_query_data(sequence: Sequence, seed: SeedSpec) -> Query {
        let cfg = SeedConfig::with_wobble(seed, MismatchSpec::exact(), false);
        Query::from_parts("q".into(), sequence, &cfg).expect("query data")
    }

    #[test]
    fn seed_sa_is_built_on_interval_slice() {
        let sequence = Sequence::from(vec![Base::A, Base::U, Base::G, Base::C, Base::A, Base::U]);
        let seed = SeedSpec::IntervalWithLength {
            start: 2,
            end: 5,
            length: 2,
        };
        let query = make_query_data(sequence, seed.clone());

        assert_eq!(query.seed_interval, 1..5);
        assert_eq!(
            query.seed_sequence().as_slice(),
            &query.sequence().as_slice()[1..5]
        );

        let expected_sa =
            SuffixArray::try_from(query.seed_sequence().as_slice()).expect("slice SA");
        assert_eq!(query.sa(), expected_sa.as_ref());
    }

    #[test]
    fn seed_sequence_tracks_only_valid_interval_starts() {
        let sequence = Sequence::from(vec![Base::A, Base::G, Base::C, Base::U, Base::A]);
        let seed = SeedSpec::LengthOnly(3);
        let query = make_query_data(sequence, seed.clone());

        assert_eq!(query.seed_interval, 0..5);
        assert_eq!(query.seed_sequence().len(), 5);
        assert_eq!(
            query.seed_sequence().as_slice(),
            query.sequence().as_slice()
        );
    }
}
