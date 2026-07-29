//! Borrowed target-index views and target-coordinate mapping.

use std::ops::Range;

use crate::index::sa::SuffixIndexView;
use crate::types::{Base, Strand};

/// Zero-copy target index used by seed search.
///
/// The paired suffix index covers the concatenation
/// `R(T) + Gap + C(T) + Gap` for every target `T`. The offset directory maps
/// positions in that global sequence back to a target, physical strand, and
/// strand-local coordinate. A block ends at the next target offset, or at the
/// end of the indexed sequence for the final target.
///
/// All fields are borrowed from the same validated [`crate::TargetRegistry`].
#[derive(Clone, Copy)]
pub struct TargetView<'a> {
    suffixes: SuffixIndexView<'a>,
    /// Start offset of each target block in the indexed sequence.
    offsets: &'a [usize],
}

impl<'a> TargetView<'a> {
    /// Construct a target view from a suffix index and its target directory.
    ///
    /// `offsets[i]` must point to the start of target `i`'s
    /// `R(T) + Gap + C(T) + Gap` block. The offsets must be strictly increasing
    /// and collectively cover the indexed sequence.
    pub(crate) fn new(suffixes: SuffixIndexView<'a>, offsets: &'a [usize]) -> Self {
        Self { suffixes, offsets }
    }

    /// Return the paired sequence and suffix positions used for traversal.
    #[inline]
    pub(crate) fn suffixes(self) -> SuffixIndexView<'a> {
        self.suffixes
    }

    /// Return a range of target positions in lexicographic suffix order.
    ///
    /// # Panics
    ///
    /// Panics if `interval` is not a valid range within the suffix index.
    #[inline]
    pub(crate) fn suffix_positions(self, interval: Range<usize>) -> &'a [u64] {
        self.suffixes.suffix_positions(interval)
    }

    /// Return the selected physical target strand in duplex-column order.
    ///
    /// For input `T` written 5' to 3', Forward is `R(T)` and Reverse is `C(T)`.
    ///
    /// # Panics
    ///
    /// Panics if `target_idx` is outside the target directory.
    #[inline]
    pub fn target(&self, target_idx: usize, strand: Strand) -> &'a [Base] {
        let seq_len = self.seq_len(target_idx);
        let offset = self.offsets[target_idx]
            + match strand {
                Strand::Forward => 0,
                Strand::Reverse => seq_len + 1,
            };
        &self.suffixes.sequence()[offset..offset + seq_len]
    }

    /// Map a global sequence position to a physical target coordinate.
    ///
    /// The returned start is relative to the physical strand selected by
    /// [`Strand`]. Returns `None` when `global_pos` is outside a physical strand,
    /// points at a separator [`Base::Gap`], or leaves fewer than `seed_len`
    /// bases before the next separator.
    #[inline]
    pub(crate) fn map_seed_pos(
        self,
        global_pos: usize,
        seed_len: usize,
    ) -> Option<(usize, Strand, usize)> {
        let target_idx = self
            .offsets
            .partition_point(|&offset| offset <= global_pos)
            .checked_sub(1)?;
        let local_pos = global_pos - self.offsets[target_idx];
        let seq_len = self.seq_len(target_idx);

        let (strand, start) = if local_pos < seq_len {
            (Strand::Forward, local_pos)
        } else {
            (Strand::Reverse, local_pos.checked_sub(seq_len + 1)?)
        };
        (start < seq_len && seed_len <= seq_len - start).then_some((target_idx, strand, start))
    }

    /// Derive one physical strand's length from its enclosing target block.
    #[inline]
    fn seq_len(self, target_idx: usize) -> usize {
        let block_start = self.offsets[target_idx];
        let block_end = self
            .offsets
            .get(target_idx + 1)
            .copied()
            .unwrap_or(self.suffixes.sequence().len());
        (block_end - block_start - 2) / 2
    }
}
