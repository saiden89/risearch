mod parallel_sa;
mod engine;

pub use engine::collect;

use crate::types::{Base, Strand};

/// Zero-copy view of a combined sequence + suffix array used by the seed search.
///
/// Both the target index and the query registry expose this view. The target
/// view covers the full fwd+rc concatenated sequences; the query view covers
/// the concatenated seed-interval subsequences.
#[derive(Clone, Copy)]
pub struct SeedView<'a> {
    pub combined_seq: &'a [Base],
    pub combined_sa: &'a [u64],
    pub sa_real_len: usize,
    /// Start offset of each entry's block in `combined_seq`.
    pub offsets: &'a [usize],
    /// Sequence length for each entry (forward only).
    pub seq_lens: &'a [usize],
}

impl<'a> SeedView<'a> {
    /// Read the base at `combined_seq[sa[sa_idx] + offset]` — branchless hot-path lookup.
    #[inline(always)]
    pub(crate) unsafe fn sa_base_unchecked(&self, sa_idx: usize, offset: usize) -> Base {
        // SAFETY: SA indices come from valid suffix array construction; SA_CHAR_PADDING
        // sentinels ensure combined_seq[suffix_pos + offset] is always in bounds.
        let suffix_pos = unsafe { *self.combined_sa.get_unchecked(sa_idx) as usize };
        unsafe { *self.combined_seq.get_unchecked(suffix_pos + offset) }
    }
}

/// A candidate seed match found during suffix array search.
#[derive(Debug, Clone)]
pub struct SeedHit {
    /// Index of query in the registry
    pub query_idx: usize,
    /// Start position in query sequence (0-based)
    pub query_start: usize,
    /// Index of target sequence in the index
    pub target_idx: usize,
    /// Start position in the strand-selected target view used by seeding/DP (0-based)
    pub target_start: usize,
    /// Length of the seed match
    pub len: usize,
    /// Strand of the match (forward or reverse)
    pub strand: Strand,
}
