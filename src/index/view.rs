use crate::types::Base;

/// Zero-copy view of a combined sequence + suffix array.
///
/// Both the target index and the query registry expose this view. The target
/// view covers the concatenated forward and reverse physical duplex blocks;
/// the query view covers the concatenated seed-interval subsequences.
#[derive(Clone, Copy)]
pub struct RegistryView<'a> {
    pub combined_seq: &'a [Base],
    pub combined_sa: &'a [u64],
    pub sa_real_len: usize,
    /// Start offset of each entry's block in `combined_seq`.
    pub offsets: &'a [usize],
    /// Sequence length for each entry (forward only).
    pub seq_lens: &'a [usize],
}

impl<'a> RegistryView<'a> {
    /// Read the base at `combined_seq[sa[sa_idx] + offset]` — branchless hot-path lookup.
    #[inline(always)]
    pub(crate) unsafe fn sa_base_unchecked(&self, sa_idx: usize, offset: usize) -> Base {
        // SAFETY: SA indices come from valid suffix array construction; SA_CHAR_PADDING
        // sentinels ensure combined_seq[suffix_pos + offset] is always in bounds.
        let suffix_pos = unsafe { *self.combined_sa.get_unchecked(sa_idx) as usize };
        unsafe { *self.combined_seq.get_unchecked(suffix_pos + offset) }
    }
}
