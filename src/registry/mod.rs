use anyhow::{anyhow, bail, Context, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::path::Path;

use crate::config::SeedConfig;
use crate::fastx::read_fasta_sequences;
use crate::index::io::validate_readable_file;
use crate::index::sa::SuffixArray;
use crate::seq::Sequence;
use crate::types::{Base, Interval};

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

/// Target data computed once at load/index-build time.
#[derive(Serialize, Deserialize)]
pub struct TargetData {
    /// Sequence identifier.
    pub name: String,
    /// Suffix array for the forward strand.
    pub forward_sa: SuffixArray,
    /// Suffix array for the reverse strand.
    pub reverse_sa: SuffixArray,
    /// Forward sequence.
    pub sequence: Sequence,
    /// Reverse-complement sequence.
    pub sequence_rc: Sequence,
}

/// Query data computed once at load time.
///
/// Owns all sequence data, suffix arrays, and N-position metadata.
/// For config-dependent views (interval bounds), use `QueryView`.
pub struct QueryData {
    /// Sequence identifier
    name: String,
    /// Forward sequence
    sequence: Sequence,
    /// Reverse complement sequence
    sequence_rc: Sequence,
    /// Suffix array for reverse complement
    reverse_sa: SuffixArray,
    /// Pre-computed seed interval bounds
    seed_interval: Interval,
    /// Minimum seed length
    min_seed_len: usize,
    /// Prefix sum of N positions for O(1) N-checking
    n_prefix: Vec<u32>,
    /// Fast path when query has no Ns
    has_n_any: bool,
}

impl QueryData {
    fn from_parts(
        name: String,
        sequence: Sequence,
        sequence_rc: Sequence,
        reverse_sa: SuffixArray,
        config: &SeedConfig,
    ) -> Result<Self> {
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
        let (start1, end1, min_seed_len) = config
            .seed
            .normalize(q_len)
            .map_err(|err| anyhow!("Invalid seed spec for query '{}': {}", name, err))?;
        let seed_interval = Interval::new(start1 - 1, end1);

        Ok(Self {
            name,
            sequence,
            sequence_rc,
            reverse_sa,
            seed_interval,
            min_seed_len,
            n_prefix,
            has_n_any,
        })
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[inline]
    pub fn sequence(&self) -> &Sequence {
        &self.sequence
    }

    #[inline]
    pub fn sequence_rc(&self) -> &Sequence {
        &self.sequence_rc
    }

    #[inline]
    pub fn reverse_sa(&self) -> &SuffixArray {
        &self.reverse_sa
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
    pub fn seed_interval(&self) -> Interval {
        self.seed_interval
    }

    #[inline]
    pub fn min_seed_len(&self) -> usize {
        self.min_seed_len
    }
}

impl RegistryEntry for QueryData {
    fn name(&self) -> &str {
        QueryData::name(self)
    }
}

pub type QueryRegistry = Registry<QueryData>;
pub type TargetRegistry = Registry<TargetData>;

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

fn normalize_record(id: String, seq: Vec<u8>) -> Result<Option<(String, Sequence, Sequence)>> {
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

    if seq_norm.len() > u32::MAX as usize {
        bail!(
            "Sequence '{}' too long for u32 suffix array ({} bases > {} max). \
             Consider chunking the sequence or using a future u64-enabled build.",
            id,
            seq_norm.len(),
            u32::MAX
        );
    }

    let seq_rc = seq_norm.reverse_complement();
    Ok(Some((id, seq_norm, seq_rc)))
}

impl QueryRegistry {
    pub fn from_fasta(path: &Path, config: &SeedConfig) -> Result<Self> {
        let sequences = read_and_validate_sequences(path)?;

        let maybe_entries: Vec<Option<QueryData>> = sequences
            .into_par_iter()
            .map(|(id, seq)| -> Result<Option<QueryData>> {
                let Some((name, sequence, sequence_rc)) = normalize_record(id, seq)? else {
                    return Ok(None);
                };

                let reverse_sa = SuffixArray::try_from(&sequence_rc)?;

                Ok(Some(QueryData::from_parts(
                    name,
                    sequence,
                    sequence_rc,
                    reverse_sa,
                    config,
                )?))
            })
            .collect::<Result<Vec<_>>>()?;

        let entries: Vec<QueryData> = maybe_entries.into_iter().flatten().collect();
        if entries.is_empty() {
            bail!(
                "All sequences were empty after normalization in {}",
                path.display()
            );
        }

        Ok(Self::new(entries))
    }
}

impl RegistryEntry for TargetData {
    fn name(&self) -> &str {
        &self.name
    }
}

impl TargetRegistry {
    pub fn from_fasta(path: &Path) -> Result<Self> {
        let sequences = read_and_validate_sequences(path)?;

        let maybe_entries: Vec<Option<TargetData>> = sequences
            .into_par_iter()
            .map(|(id, seq)| -> Result<Option<TargetData>> {
                let Some((name, sequence, sequence_rc)) = normalize_record(id, seq)? else {
                    return Ok(None);
                };

                let forward_sa = SuffixArray::try_from(&sequence)?;
                let reverse_sa = SuffixArray::try_from(&sequence_rc)?;

                Ok(Some(TargetData {
                    name,
                    forward_sa,
                    reverse_sa,
                    sequence,
                    sequence_rc,
                }))
            })
            .collect::<Result<Vec<_>>>()?;

        let entries: Vec<TargetData> = maybe_entries.into_iter().flatten().collect();
        if entries.is_empty() {
            bail!(
                "All sequences were empty after normalization in {}",
                path.display()
            );
        }

        Ok(Self::new(entries))
    }

    pub fn load(path: &Path) -> Result<Self> {
        crate::index::io::load_index_file(path)
    }

    pub fn save(&self, path: &Path) -> Result<()> {
        crate::index::io::write_index_file(self, path)
    }

    pub fn get_sequence(&self, seq_idx: usize) -> &Sequence {
        &self.entries[seq_idx].sequence
    }

    pub fn get_sequence_rc(&self, seq_idx: usize) -> &Sequence {
        &self.entries[seq_idx].sequence_rc
    }

    pub fn get_sequence_len(&self, seq_idx: usize) -> usize {
        self.entries[seq_idx].sequence.len()
    }
}
