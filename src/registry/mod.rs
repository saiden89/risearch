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
    /// Combined sequence: forward ++ [Gap] ++ reverse-complement
    pub combined_seq: Sequence,
    /// Suffix array built on combined_seq
    pub combined_sa: SuffixArray,
    /// Length of the original forward sequence
    pub seq_len: usize,
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
    /// Suffix array for forward sequence (bit-packed u64)
    sa: SuffixArray,
    /// Pre-computed seed interval bounds
    seed_interval: Interval,
    /// Minimum seed length
    min_seed_len: usize,
    /// Prefix sum of N positions for O(1) N-checking
    n_prefix: Vec<u32>,
    /// Fast path when query has no Ns
    has_n_any: bool,
    /// Reversed sequence (NOT reverse-complement)
    sequence_rc: Sequence,
    /// Suffix array for reversed sequence
    reverse_sa: SuffixArray,
}

impl QueryData {
    fn from_parts(
        name: String,
        sequence: Sequence,
        sa: SuffixArray,
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

        let reverse_seq_bases: Vec<Base> = sequence.as_ref().iter().rev().copied().collect();
        let sequence_rc = Sequence::from(reverse_seq_bases);
        let reverse_sa = SuffixArray::try_from(&sequence_rc)
            .map_err(|e| anyhow!("Failed to build reverse SA for query '{}': {}", name, e))?;

        Ok(Self {
            name,
            sequence,
            sa,
            seed_interval,
            min_seed_len,
            n_prefix,
            has_n_any,
            sequence_rc,
            reverse_sa,
        })
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[inline(always)]
    pub fn sequence(&self) -> &[Base] {
        self.sequence.as_ref()
    }

    #[inline(always)]
    pub fn sa(&self) -> &[u64] {
        self.sa.as_ref()
    }

    #[inline(always)]
    pub fn sequence_rc(&self) -> &[Base] {
        self.sequence_rc.as_ref()
    }

    #[inline(always)]
    pub fn reverse_sa(&self) -> &[u64] {
        self.reverse_sa.as_ref()
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

        let maybe_entries: Vec<Option<QueryData>> = sequences
            .into_par_iter()
            .map(|(id, seq)| -> Result<Option<QueryData>> {
                let Some((name, sequence)) = normalize_record(id, seq)? else {
                    return Ok(None);
                };

                // Build SA on forward sequence
                let sa = SuffixArray::try_from(&sequence)?;

                Ok(Some(QueryData::from_parts(name, sequence, sa, config)?))
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
                let Some((name, sequence)) = normalize_record(id, seq)? else {
                    return Ok(None);
                };

                // Build combined sequence: fwd ++ [Gap] ++ rc
                let seq_len = sequence.len();
                let sequence_rc = sequence.reverse_complement();
                let mut combined_bases: Vec<Base> = Vec::with_capacity(2 * seq_len + 2);
                combined_bases.extend_from_slice(&sequence);
                combined_bases.push(Base::Gap);
                combined_bases.extend_from_slice(&sequence_rc);
                let combined_seq = Sequence::from(combined_bases.clone());
                let combined_sa = SuffixArray::try_from(&combined_seq)?;

                Ok(Some(TargetData {
                    name,
                    combined_seq,
                    combined_sa,
                    seq_len,
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

    pub fn get_sequence(&self, seq_idx: usize) -> &[Base] {
        let t = &self.entries[seq_idx];
        &t.combined_seq[..t.seq_len]
    }

    pub fn get_sequence_rc(&self, seq_idx: usize) -> &[Base] {
        let t = &self.entries[seq_idx];
        &t.combined_seq[t.seq_len + 1..2 * t.seq_len + 1]
    }

    pub fn get_sequence_len(&self, seq_idx: usize) -> usize {
        self.entries[seq_idx].seq_len
    }
}
