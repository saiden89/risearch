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
use crate::index::store::SA_CHAR_PADDING;
use crate::seq::Sequence;
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
        let (start1, end1, min_seed_len) = config
            .seed
            .normalize(q_len)
            .map_err(|err| anyhow!("Invalid seed spec for query '{}': {}", name, err))?;
        let seed_interval = (start1 - 1)..end1;
        let max_seed_len = seed_interval.end.saturating_sub(seed_interval.start);
        let seed_sequence =
            Sequence::from(sequence[seed_interval.start..seed_interval.end].to_vec());

        Ok(Self {
            name,
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
}

impl RegistryEntry for Query {
    fn name(&self) -> &str {
        Query::name(self)
    }
}

/// Combined query data for seed search — mirrors `TargetView` on the target side.
#[derive(Clone, Copy)]
pub struct QueryView<'a> {
    pub combined_seed_seq: &'a [Base],
    pub combined_sa: &'a [u64],
    pub len: usize,
    pub offsets: &'a [usize],
    pub seed_seq_lens: &'a [usize],
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

    pub fn view(&self) -> QueryView<'_> {
        QueryView {
            combined_seed_seq: &self.combined_seed_seq,
            combined_sa: &self.combined_sa,
            len: self.len,
            offsets: &self.offsets,
            seed_seq_lens: &self.seed_seq_lens,
        }
    }
}

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
        Self::from_fastas(&[path], config).with_context(|| format!("in {}", path.display()))
    }

    /// Load queries from one or more FASTA files.
    ///
    /// Sequences are validated per-file (duplicate IDs within a file are
    /// rejected) and then checked for cross-file duplicates before building
    /// the registry. Query SA construction is parallelised via rayon over
    /// the merged sequence list.
    pub fn from_fastas(paths: &[&Path], config: &SeedConfig) -> Result<Self> {
        if paths.is_empty() {
            bail!("No query files provided");
        }

        // Read and validate each file; collect all sequences into one vec.
        let mut all_sequences: Vec<(String, Vec<u8>)> = Vec::new();
        for path in paths {
            let seqs = read_and_validate_sequences(path)?;
            all_sequences.extend(seqs);
        }

        // Check for duplicate IDs across all files.
        let mut seen = HashSet::with_capacity(all_sequences.len());
        for (id, _) in &all_sequences {
            if !seen.insert(id.as_str()) {
                bail!("Duplicate FASTA record id '{}' across input files", id);
            }
        }

        let maybe_entries: Vec<Option<Query>> = all_sequences
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
        let sa_real_len = combined_sa_raw.len();

        // Pad both for branchless SA character lookup.
        combined_seed_seq.resize(combined_seed_seq.len() + SA_CHAR_PADDING, Base::Gap);
        let mut combined_sa = combined_sa_raw.into_inner();
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{MismatchSpec, SeedSpec};
    use std::io::Write;

    fn seed_config() -> SeedConfig {
        SeedConfig::with_wobble(SeedSpec::LengthOnly(4), MismatchSpec::exact(), true)
    }

    fn temp_fasta(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    #[test]
    fn from_fastas_single_file_matches_from_fasta() {
        let f = temp_fasta(">seq1\nACGUACGU\n>seq2\nUUUUAAAA\n");
        let cfg = seed_config();

        let r1 = QueryRegistry::from_fasta(f.path(), &cfg).unwrap();
        let r2 = QueryRegistry::from_fastas(&[f.path()], &cfg).unwrap();

        assert_eq!(r1.len(), r2.len());
        for i in 0..r1.len() {
            assert_eq!(r1.get_name(i), r2.get_name(i));
        }
    }

    #[test]
    fn from_fastas_merges_two_files() {
        let f1 = temp_fasta(">seq1\nACGUACGU\n");
        let f2 = temp_fasta(">seq2\nUUUUAAAA\n");
        let cfg = seed_config();

        let registry = QueryRegistry::from_fastas(&[f1.path(), f2.path()], &cfg).unwrap();

        assert_eq!(registry.len(), 2);
        let names: Vec<&str> = (0..2).map(|i| registry.get_name(i)).collect();
        assert!(names.contains(&"seq1") && names.contains(&"seq2"));
    }

    #[test]
    fn from_fastas_preserves_order_across_files() {
        let f1 = temp_fasta(">alpha\nACGUACGU\n");
        let f2 = temp_fasta(">beta\nUUUUAAAA\n>gamma\nGGGGCCCC\n");
        let cfg = seed_config();

        let registry = QueryRegistry::from_fastas(&[f1.path(), f2.path()], &cfg).unwrap();

        assert_eq!(registry.len(), 3);
        assert_eq!(registry.get_name(0), "alpha");
        assert_eq!(registry.get_name(1), "beta");
        assert_eq!(registry.get_name(2), "gamma");
    }

    #[test]
    fn from_fastas_rejects_cross_file_duplicates() {
        let f1 = temp_fasta(">seq1\nACGUACGU\n");
        let f2 = temp_fasta(">seq1\nUUUUAAAA\n");
        let cfg = seed_config();

        let result = QueryRegistry::from_fastas(&[f1.path(), f2.path()], &cfg);
        assert!(result.is_err());
        let msg = result.err().unwrap().to_string();
        assert!(
            msg.contains("Duplicate"),
            "expected duplicate error, got: {msg}"
        );
    }

    #[test]
    fn from_fastas_empty_paths_returns_error() {
        let cfg = seed_config();
        assert!(QueryRegistry::from_fastas(&[], &cfg).is_err());
    }

    fn make_query_data(sequence: Sequence, seed: SeedSpec) -> Query {
        let cfg = SeedConfig::with_wobble(seed, MismatchSpec::exact(), false);
        Query::from_parts("q".into(), sequence, &cfg).expect("query data")
    }

    #[test]
    fn seed_sequence_matches_interval_slice() {
        let sequence = Sequence::from(vec![Base::A, Base::U, Base::G, Base::C, Base::A, Base::U]);
        let seed = SeedSpec::Interval {
            start: 2,
            end: 5,
            length: Some(2),
        };
        let query = make_query_data(sequence, seed.clone());

        assert_eq!(query.seed_interval, 1..5);
        assert_eq!(query.min_seed_len, 2);
        assert_eq!(query.max_seed_len, 4);
        assert_eq!(query.seed_sequence(), &query.sequence()[1..5]);
    }

    #[test]
    fn seed_sequence_tracks_only_valid_interval_starts() {
        let sequence = Sequence::from(vec![Base::A, Base::G, Base::C, Base::U, Base::A]);
        let seed = SeedSpec::LengthOnly(3);
        let query = make_query_data(sequence, seed.clone());

        assert_eq!(query.seed_interval, 0..5);
        assert_eq!(query.min_seed_len, 3);
        assert_eq!(query.max_seed_len, 5);
        assert_eq!(query.seed_sequence().len(), 5);
        assert_eq!(query.seed_sequence(), query.sequence());
    }
}
