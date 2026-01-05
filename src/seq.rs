//! Sequence abstraction for clean base access.
//!
//! This module provides a `Seq` wrapper that encapsulates sequence data and strand,
//! providing clean methods for base access without scattered `Base::from_byte()` calls.

use crate::types::{Base, Strand};

/// A sequence wrapper that provides clean base access.
///
/// Wraps raw sequence bytes with strand information, providing methods
/// for type-safe base access that eliminate scattered `Base::from_byte()` calls.
#[derive(Debug, Clone, Copy)]
pub struct Seq<'a> {
    data: &'a [u8],
    strand: Strand,
}

impl<'a> Seq<'a> {
    /// Create a new Seq from raw bytes and strand information.
    #[inline]
    pub fn new(data: &'a [u8], strand: Strand) -> Self {
        Self { data, strand }
    }

    /// Create a forward strand Seq (most common case).
    #[inline]
    pub fn forward(data: &'a [u8]) -> Self {
        Self::new(data, Strand::Forward)
    }

    /// Length of the underlying sequence.
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Check if sequence is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get base at position as Base enum.
    ///
    /// # Panics
    /// Panics if pos >= len().
    #[inline]
    pub fn base(&self, pos: usize) -> Base {
        Base::from_byte(self.data[pos])
    }

    /// Get base at position, or Gap if out of bounds.
    ///
    /// Useful in DP where extensions can go beyond sequence boundaries.
    #[inline]
    pub fn base_or_gap(&self, pos: usize) -> Base {
        if pos >= self.data.len() {
            Base::Gap
        } else {
            Base::from_byte(self.data[pos])
        }
    }

    /// Get DSM index at position (0-5 range).
    ///
    /// This is `base(pos).idx()` - a common pattern in energy calculations.
    #[inline]
    pub fn dsm_idx(&self, pos: usize) -> usize {
        self.base(pos).idx()
    }

    /// Get raw byte at position.
    #[inline]
    pub fn byte(&self, pos: usize) -> u8 {
        self.data[pos]
    }

    /// Get the strand of this sequence.
    #[inline]
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Get base going LEFT (5' direction) from anchor.
    ///
    /// Used in dp_left where query goes q_start-i.
    /// Returns Base::Gap if offset would go before position 0.
    #[inline]
    pub fn left(&self, anchor: usize, offset: usize) -> Base {
        if offset > anchor {
            Base::Gap
        } else {
            self.base(anchor - offset)
        }
    }

    /// Get base going RIGHT (3' direction) from anchor.
    ///
    /// Used in dp_right where query goes q_end+i.
    /// Returns Base::Gap if position would exceed sequence length.
    #[inline]
    pub fn right(&self, anchor: usize, offset: usize) -> Base {
        self.base_or_gap(anchor + offset)
    }

    /// Access underlying bytes slice.
    #[inline]
    pub fn as_bytes(&self) -> &[u8] {
        self.data
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seq_base_access() {
        let data = b"ACGU";
        let seq = Seq::forward(data);

        assert_eq!(seq.base(0), Base::A);
        assert_eq!(seq.base(1), Base::C);
        assert_eq!(seq.base(2), Base::G);
        assert_eq!(seq.base(3), Base::U);
    }

    #[test]
    fn test_seq_base_or_gap() {
        let data = b"ACG";
        let seq = Seq::forward(data);

        assert_eq!(seq.base_or_gap(0), Base::A);
        assert_eq!(seq.base_or_gap(2), Base::G);
        assert_eq!(seq.base_or_gap(3), Base::Gap); // Out of bounds
        assert_eq!(seq.base_or_gap(100), Base::Gap);
    }

    #[test]
    fn test_seq_left_right() {
        let data = b"ACGU";
        let seq = Seq::forward(data);

        // left from position 2
        assert_eq!(seq.left(2, 0), Base::G); // position 2
        assert_eq!(seq.left(2, 1), Base::C); // position 1
        assert_eq!(seq.left(2, 2), Base::A); // position 0
        assert_eq!(seq.left(2, 3), Base::Gap); // would be -1

        // right from position 1
        assert_eq!(seq.right(1, 0), Base::C); // position 1
        assert_eq!(seq.right(1, 1), Base::G); // position 2
        assert_eq!(seq.right(1, 2), Base::U); // position 3
        assert_eq!(seq.right(1, 3), Base::Gap); // position 4, OOB
    }

    #[test]
    fn test_seq_dsm_idx() {
        let data = b"ACGU";
        let seq = Seq::forward(data);

        assert_eq!(seq.dsm_idx(0), Base::A.idx());
        assert_eq!(seq.dsm_idx(1), Base::C.idx());
    }
}
