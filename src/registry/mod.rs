use serde::{Deserialize, Serialize};

use crate::config::SeedConfig;
use crate::index::sa::SequenceIndex;
use crate::sa::SuffixArray;
use crate::seq::Sequence;
use crate::types::Base;

#[derive(Clone, Copy, Debug)]
pub struct SeedInterval {
    /// Start position (0-based, inclusive)
    pub start: usize,
    /// End position (0-based, exclusive)
    pub end: usize,
    /// Minimum seed length
    pub min_len: usize,
}

impl SeedInterval {
    /// Compute interval from SeedConfig and query length.
    ///
    /// Converts from 1-based (SeedSpec) to 0-based indexing.
    pub fn from_config(config: &SeedConfig, query_len: usize) -> Self {
        match config.seed.normalize(query_len) {
            Ok((start1, end1, min_len)) => Self {
                start: start1 - 1, // Convert to 0-based
                end: end1,         // end1 is already exclusive in 0-based terms
                min_len,
            },
            Err(_) => {
                // Fallback for invalid spec: use entire sequence
                Self {
                    start: 0,
                    end: query_len,
                    min_len: query_len,
                }
            }
        }
    }
}

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
    seed_interval: SeedInterval,
    /// Prefix sum of N positions for O(1) N-checking
    n_prefix: Vec<u32>,
    /// Fast path when query has no Ns
    has_n_any: bool,
}

impl QueryData {
    /// Build QueryData from a SequenceIndex, computing N-prefix and seed interval.
    pub fn from_index(index: SequenceIndex, config: &SeedConfig) -> Self {
        let q_len = index.sequence.len();

        // Compute N-prefix for O(1) N-checking
        let mut n_prefix = Vec::with_capacity(q_len + 1);
        n_prefix.push(0);
        let mut n_total = 0;
        for &base in index.sequence.iter() {
            if base == Base::N {
                n_total += 1;
            }
            n_prefix.push(n_total);
        }
        let has_n_any = n_total != 0;

        // Compute seed interval once
        let seed_interval = SeedInterval::from_config(config, q_len);

        Self {
            name: index.name,
            sequence: index.sequence,
            sequence_rc: index.sequence_rc,
            reverse_sa: index.reverse_sa,
            seed_interval,
            n_prefix,
            has_n_any,
        }
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
    pub fn seed_interval(&self) -> SeedInterval {
        self.seed_interval
    }
}

impl RegistryEntry for QueryData {
    fn name(&self) -> &str {
        QueryData::name(self)
    }
}

pub type QueryRegistry = Registry<QueryData>;
pub type TargetRegistry = Registry<SequenceIndex>;

impl QueryRegistry {
    pub fn from_indices(indices: Vec<SequenceIndex>, config: &SeedConfig) -> Self {
        Self::new(
            indices
                .into_iter()
                .map(|idx| QueryData::from_index(idx, config))
                .collect(),
        )
    }

    pub fn from_names(names: Vec<String>, config: &SeedConfig) -> Self {
        let entries = names
            .into_iter()
            .map(|name| {
                QueryData::from_index(
                    SequenceIndex {
                        name,
                        forward_sa: SuffixArray::from(Vec::new()),
                        reverse_sa: SuffixArray::from(Vec::new()),
                        sequence: Sequence::from(Vec::new()),
                        sequence_rc: Sequence::from(Vec::new()),
                    },
                    config,
                )
            })
            .collect();
        Self { entries }
    }
}

impl RegistryEntry for SequenceIndex {
    fn name(&self) -> &str {
        &self.name
    }
}

impl TargetRegistry {
    pub fn get_sequence(&self, seq_idx: usize) -> &Sequence {
        &self.entries[seq_idx].sequence
    }

    pub fn get_sequence_rc(&self, seq_idx: usize) -> &Sequence {
        &self.entries[seq_idx].sequence_rc
    }

    pub fn get_sequence_len(&self, seq_idx: usize) -> usize {
        self.entries[seq_idx].sequence.len()
    }

    pub fn from_names(names: Vec<String>) -> Self {
        let entries = names
            .into_iter()
            .map(|name| SequenceIndex {
                name,
                forward_sa: SuffixArray::from(Vec::new()),
                reverse_sa: SuffixArray::from(Vec::new()),
                sequence: Sequence::from(Vec::new()),
                sequence_rc: Sequence::from(Vec::new()),
            })
            .collect();
        Self { entries }
    }
}
