use anyhow::{anyhow, bail, Context, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::ops::Range;
use std::path::Path;

use crate::config::SeedConfig;
use crate::fastx::{read_and_validate_fasta, normalize_record};
use crate::index::io::validate_readable_file;
use crate::index::sa::SuffixArray;
use crate::index::store::SA_CHAR_PADDING;
use crate::seed::SeedView;
use crate::seq::Sequence;
use crate::types::Base;

pub trait RegistryEntry {
    fn name(&self) -> &str;
}

impl RegistryEntry for String {
    fn name(&self) -> &str {
        self.as_str()
    }
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

    pub fn get(&self, idx: usize) -> &T {
        &self.entries[idx]
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
    pub fn get_name(&self, idx: usize) -> &str {
        debug_assert!(
            idx < self.entries.len(),
            "Registry index out of bounds: {} >= {}",
            idx,
            self.entries.len()
        );
        // SAFETY: idx comes from generated hits, guaranteed to be valid
        unsafe { self.entries.get_unchecked(idx).name() }
    }

    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.entries.iter().position(|e| e.name() == name)
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, &T)> {
        self.entries.iter().enumerate()
    }
}

/// Query prepared once at load time for the current seed configuration.
pub struct Query {
    /// Sequence identifier
    name: String,
    /// Full forward sequence used for extension and output
    sequence: Sequence,
    /// Seed-interval slice used by the combined seeding suffix array
    seed_sequence: Sequence,
    /// Pre-computed normalized seed interval bounds on the full query
    pub(crate) seed_interval: Range<usize>,
    /// Minimum admissible seed length for this query after normalization.
    pub(crate) min_seed_len: usize,
    /// Maximum admissible seed length for this query after normalization.
    pub(crate) max_seed_len: usize,
    /// Prefix sum of N positions for O(1) N-checking
    n_prefix: Vec<u32>,
    /// Fast path when query has no Ns
    has_n_any: bool,
}

impl Query {
    fn from_parts(id: String, sequence: Sequence, config: &SeedConfig) -> Result<Self> {
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
            .resolve(q_len)
            .map_err(|err| anyhow!("Invalid seed spec for query '{}': {}", id, err))?;
        let seed_interval = (start1 - 1)..end1;
        let max_seed_len = seed_interval.end.saturating_sub(seed_interval.start);
        let seed_sequence =
            Sequence::from(sequence[seed_interval.start..seed_interval.end].to_vec());

        Ok(Self {
            name: id,
            sequence,
            seed_sequence,
            seed_interval,
            min_seed_len,
            max_seed_len,
            n_prefix,
            has_n_any,
        })
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[inline(always)]
    pub fn sequence(&self) -> &[Base] {
        &self.sequence
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
    pub fn seed_sequence(&self) -> &[Base] {
        &self.seed_sequence
    }

    /// Map a local position within the seed view to a global query start coordinate.
    ///
    /// Returns `None` if the seed overflows the query's seed-search interval,
    /// violates the query-specific length constraints, or contains an 'N' base.
    #[inline]
    pub fn map_seed_pos(&self, local_pos: usize, seed_len: usize) -> Option<usize> {
        if local_pos + seed_len > self.seed_sequence.len() {
            return None;
        }

        if seed_len < self.min_seed_len || seed_len > self.max_seed_len {
            return None;
        }

        let query_start = self.seed_interval.start + local_pos;

        if self.has_n_any && self.n_prefix[query_start + seed_len] != self.n_prefix[query_start] {
            return None;
        }

        Some(query_start)
    }
}

impl RegistryEntry for Query {
    fn name(&self) -> &str {
        Query::name(self)
    }
}

/// Registry of queries with a combined suffix array for efficient seed search.
///
/// Individual `Query` entries hold per-query metadata (sequence, seed interval,
/// N-prefix). The combined SA/sequence enables a single recursive traversal
/// across all queries × all targets, avoiding redundant target-side partitioning.
pub struct QueryRegistry {
    inner: Registry<Query>,
    /// Concatenated seed sequences with Gap separators + SA_CHAR_PADDING sentinels.
    combined_seed_seq: Vec<Base>,
    /// Global SA over combined_seed_seq + SA_CHAR_PADDING sentinel zeros.
    combined_sa: Vec<u64>,
    /// Number of real SA entries (excluding padding).
    len: usize,
    /// Start offset of each query's seed sequence in combined_seed_seq.
    offsets: Vec<usize>,
    /// Length of each query's seed sequence.
    seed_seq_lens: Vec<usize>,
}

impl QueryRegistry {
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    pub fn get(&self, idx: usize) -> &Query {
        self.inner.get(idx)
    }

    pub fn entries(&self) -> &[Query] {
        self.inner.entries()
    }

    pub fn get_name(&self, idx: usize) -> &str {
        self.inner.get_name(idx)
    }

    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.inner.index_of(name)
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, &Query)> {
        self.inner.iter()
    }

    pub fn view(&self) -> SeedView<'_> {
        SeedView {
            combined_seq: &self.combined_seed_seq,
            combined_sa: &self.combined_sa,
            sa_real_len: self.len,
            offsets: &self.offsets,
            seq_lens: &self.seed_seq_lens,
        }
    }
}

fn read_and_validate_sequences(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    validate_readable_file(path)?;

    let sequences = read_and_validate_fasta(path)?;
    if sequences.is_empty() {
        bail!("No sequences found in input file: {}", path.display());
    }

    Ok(sequences)
}

impl QueryRegistry {
    /// Load queries from a FASTA file.
    ///
    /// Sequences are validated (duplicate IDs are rejected).
    /// Query SA construction is parallelised via rayon.
    pub fn from_fasta(path: &Path, config: &SeedConfig) -> Result<Self> {
        let seqs = read_and_validate_sequences(path)
            .with_context(|| format!("in {}", path.display()))?;

        let maybe_entries: Vec<Option<Query>> = seqs
            .into_par_iter()
            .map(|(id, seq)| -> Result<Option<Query>> {
                let Some(sequence) = normalize_record(&id, &seq)? else {
                    return Ok(None);
                };
                Ok(Some(Query::from_parts(id, sequence, config)?))
            })
            .collect::<Result<Vec<_>>>()?;

        let entries: Vec<Query> = maybe_entries.into_iter().flatten().collect();
        if entries.is_empty() {
            bail!("All sequences were empty after normalization");
        }

        // Build combined seed sequence and SA (mirrors TargetStore::build_from_fasta).
        let mut combined_seed_seq: Vec<Base> = Vec::new();
        let mut offsets: Vec<usize> = Vec::with_capacity(entries.len());
        let mut seed_seq_lens: Vec<usize> = Vec::with_capacity(entries.len());

        for query in &entries {
            offsets.push(combined_seed_seq.len());
            seed_seq_lens.push(query.seed_sequence().len());
            combined_seed_seq.extend_from_slice(query.seed_sequence());
            combined_seed_seq.push(Base::Gap);
        }

        let combined_sa_raw = SuffixArray::try_from(combined_seed_seq.as_slice())
            .context("Failed to build combined query SA")?;
        let mut combined_sa = filter_seedable_query_suffixes(
            combined_sa_raw.into_inner(),
            &entries,
            &offsets,
            &seed_seq_lens,
        );
        let sa_real_len = combined_sa.len();

        // Pad both for branchless SA character lookup.
        combined_seed_seq.resize(combined_seed_seq.len() + SA_CHAR_PADDING, Base::Gap);
        combined_sa.resize(combined_sa.len() + SA_CHAR_PADDING, 0u64);

        Ok(Self {
            inner: Registry::new(entries),
            combined_seed_seq,
            combined_sa,
            len: sa_real_len,
            offsets,
            seed_seq_lens,
        })
    }
}

fn filter_seedable_query_suffixes(
    sa: Vec<u64>,
    entries: &[Query],
    offsets: &[usize],
    seed_seq_lens: &[usize],
) -> Vec<u64> {
    sa.into_iter()
        .filter(|&suffix_pos| {
            let suffix_pos = suffix_pos as usize;
            let query_idx = offsets.partition_point(|&offset| offset <= suffix_pos) - 1;
            let min_seed_len = entries[query_idx].min_seed_len;
            let seed_seq_len = seed_seq_lens[query_idx];
            let local_pos = suffix_pos - offsets[query_idx];
            min_seed_len * 2 <= seed_seq_len + 1 || local_pos + min_seed_len <= seed_seq_len
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn seed_config() -> SeedConfig {
        SeedConfig {
            seed_start: None,
            seed_end: None,
            seed_length: Some(4),
            seed_wobble: true,
            max_mismatches: 0,
            min_prefix_matches: 1,
            min_suffix_matches: 0,
        }
    }

    fn temp_fasta(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }


    fn make_query_data(sequence: Sequence, seed_start: Option<i64>, seed_end: Option<i64>, seed_length: Option<i64>) -> Query {
        let cfg = SeedConfig {
            seed_start,
            seed_end,
            seed_length,
            seed_wobble: false,
            max_mismatches: 0,
            min_prefix_matches: 1,
            min_suffix_matches: 0,
        };
        Query::from_parts("q".into(), sequence, &cfg).expect("query data")
    }

    #[test]
    fn seed_sequence_matches_interval_slice() {
        let sequence = Sequence::from(vec![Base::A, Base::U, Base::G, Base::C, Base::A, Base::U]);
        let query = make_query_data(sequence, Some(2), Some(5), Some(2));

        assert_eq!(query.seed_interval, 1..5);
        assert_eq!(query.min_seed_len, 2);
        assert_eq!(query.max_seed_len, 4);
        assert_eq!(query.seed_sequence(), &query.sequence()[1..5]);
    }

    #[test]
    fn seed_sequence_tracks_only_valid_interval_starts() {
        let sequence = Sequence::from(vec![Base::A, Base::G, Base::C, Base::U, Base::A]);
        let query = make_query_data(sequence, None, None, Some(3));

        assert_eq!(query.seed_interval, 0..5);
        assert_eq!(query.min_seed_len, 3);
        assert_eq!(query.max_seed_len, 5);
        assert_eq!(query.seed_sequence().len(), 5);
        assert_eq!(query.seed_sequence(), query.sequence());
    }

    #[test]
    fn query_sa_keeps_only_suffixes_that_can_reach_min_seed_len() {
        let f = temp_fasta(">q1\nACGUACGUACGUACGUACGUAC\n>q2\nUGCAUGCAUGCAUGCAUGCAUG\n");
        let cfg = SeedConfig {
            seed_start: None,
            seed_end: None,
            seed_length: Some(22),
            seed_wobble: false,
            max_mismatches: 0,
            min_prefix_matches: 1,
            min_suffix_matches: 0,
        };

        let registry = QueryRegistry::from_fasta(f.path(), &cfg).unwrap();
        let view = registry.view();

        assert_eq!(view.sa_real_len, 2);
        let mut starts: Vec<usize> = view.combined_sa[..view.sa_real_len]
            .iter()
            .map(|&pos| pos as usize)
            .collect();
        starts.sort_unstable();
        assert_eq!(starts, vec![0, 23]);
    }
}
