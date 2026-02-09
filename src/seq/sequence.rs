//! Sequence newtype wrapping Vec<Base>
//!
//! This module provides the core `Sequence` type that wraps `Vec<Base>`,
//! ensuring sequences are normalized once at the input boundary and then
//! work with Base enums throughout the codebase.

use crate::types::Base;
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::ops::{Deref, Index, Range, RangeFrom, RangeFull, RangeTo};

use super::normalize::{normalize_rna_sequence, NormalizationStats};

/// A normalized RNA sequence stored as Vec<Base>.
///
/// Sequences are normalized once at input (FASTA parsing) and stored as Base enums.
/// This eliminates repeated ASCII→Base conversions during search and alignment.
#[derive(
    Clone, Debug, PartialEq, Eq, Serialize, Deserialize, wincode::SchemaWrite, wincode::SchemaRead,
)]
pub struct Sequence(Vec<Base>);

impl Sequence {
    /// Normalize raw FASTA bytes into a Sequence.
    ///
    /// This is the only place where ASCII→Base conversion happens.
    /// Input bytes are validated and converted to uppercase Base values,
    /// with gaps removed and ambiguous bases mapped to N.
    ///
    /// # Arguments
    /// * `id` - Sequence identifier (for error messages)
    /// * `raw` - Raw FASTA sequence bytes
    ///
    /// # Returns
    /// A normalized Sequence and statistics about the normalization process
    pub fn normalize(id: &str, raw: &[u8]) -> Result<(Self, NormalizationStats)> {
        let (bases, stats) = normalize_rna_sequence(id, raw)?;
        Ok((Self(bases), stats))
    }

    /// Length of the sequence
    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Check if sequence is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Compute reverse complement of this sequence.
    ///
    /// Uses Base::complement() for each base, then reverses.
    pub fn reverse_complement(&self) -> Sequence {
        let rc_bases: Vec<Base> = self.0.iter().rev().map(|&base| base.complement()).collect();

        Self(rc_bases)
    }

    /// Convert to bytes for suffix array construction.
    ///
    /// This returns discriminant values (0-5), NOT ASCII bytes.
    /// This is done once at index build time for libsais.
    ///
    /// # Performance
    /// O(n) operation, should only be called once per sequence during index building.
    pub fn to_bytes(&self) -> Vec<u8> {
        self.0.iter().map(|&b| b as u8).collect()
    }

    /// Convert to ASCII bytes for backward compatibility.
    ///
    /// Returns lowercase ASCII representation (a, g, c, t, n, -).
    /// This is a temporary method for Phase 3 compatibility.
    /// Will be removed once all code is updated to work with &[Base].
    ///
    /// # Performance
    /// O(n) operation, allocates a new Vec.
    pub fn to_ascii_bytes(&self) -> Vec<u8> {
        self.0.iter().map(|&b| b.to_byte()).collect()
    }

    /// Iterator over bases
    pub fn iter(&self) -> impl Iterator<Item = &Base> {
        self.0.iter()
    }

    /// Unchecked access to base at index (for DP hot path).
    ///
    /// # Safety
    /// Caller must ensure `index < self.len()`.
    /// Use only in performance-critical code where bounds are pre-validated.
    #[inline(always)]
    pub unsafe fn get_unchecked(&self, index: usize) -> Base {
        // SAFETY: Caller guarantees index is in bounds
        unsafe { *self.0.get_unchecked(index) }
    }
}

// Implement Deref to allow treating Sequence as &[Base]
impl Deref for Sequence {
    type Target = [Base];

    fn deref(&self) -> &[Base] {
        &self.0
    }
}

// Implement Index to allow seq[i] access
impl Index<usize> for Sequence {
    type Output = Base;

    fn index(&self, index: usize) -> &Base {
        &self.0[index]
    }
}

// Implement Index for ranges
impl Index<Range<usize>> for Sequence {
    type Output = [Base];

    fn index(&self, index: Range<usize>) -> &[Base] {
        &self.0[index]
    }
}

impl Index<RangeFrom<usize>> for Sequence {
    type Output = [Base];

    fn index(&self, index: RangeFrom<usize>) -> &[Base] {
        &self.0[index]
    }
}

impl Index<RangeTo<usize>> for Sequence {
    type Output = [Base];

    fn index(&self, index: RangeTo<usize>) -> &[Base] {
        &self.0[index]
    }
}

impl Index<RangeFull> for Sequence {
    type Output = [Base];

    fn index(&self, index: RangeFull) -> &[Base] {
        &self.0[index]
    }
}

impl From<Vec<Base>> for Sequence {
    fn from(bases: Vec<Base>) -> Self {
        Self(bases)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_simple() {
        let (seq, stats) = Sequence::normalize("test", b"ACGU").unwrap();
        assert_eq!(seq.len(), 4);
        assert_eq!(seq[0], Base::A);
        assert_eq!(seq[1], Base::C);
        assert_eq!(seq[2], Base::G);
        assert_eq!(seq[3], Base::U);
        assert_eq!(stats.removed_gaps, 0);
        assert_eq!(stats.converted_to_n, 0);
    }

    #[test]
    fn test_normalize_with_gaps() {
        let (seq, stats) = Sequence::normalize("test", b"AC-GU").unwrap();
        assert_eq!(seq.len(), 4); // Gap removed
        assert_eq!(seq[0], Base::A);
        assert_eq!(seq[1], Base::C);
        assert_eq!(seq[2], Base::G);
        assert_eq!(seq[3], Base::U);
        assert_eq!(stats.removed_gaps, 1);
    }

    #[test]
    fn test_normalize_ambiguous() {
        let (seq, stats) = Sequence::normalize("test", b"ACRGU").unwrap();
        assert_eq!(seq.len(), 5);
        assert_eq!(seq[0], Base::A);
        assert_eq!(seq[1], Base::C);
        assert_eq!(seq[2], Base::N); // R -> N
        assert_eq!(seq[3], Base::G);
        assert_eq!(seq[4], Base::U);
        assert_eq!(stats.converted_to_n, 1);
    }

    #[test]
    fn test_reverse_complement() {
        let (seq, _) = Sequence::normalize("test", b"ACGU").unwrap();
        let rc = seq.reverse_complement();

        assert_eq!(rc.len(), 4);
        assert_eq!(rc[0], Base::A); // U -> A
        assert_eq!(rc[1], Base::C); // G -> C
        assert_eq!(rc[2], Base::G); // C -> G
        assert_eq!(rc[3], Base::U); // A -> U
    }

    #[test]
    fn test_to_bytes() {
        let (seq, _) = Sequence::normalize("test", b"ACGU").unwrap();
        let bytes = seq.to_bytes();

        // Should be discriminant values, not ASCII
        assert_eq!(bytes, vec![1, 3, 2, 4]); // A=1, C=3, G=2, U=4
    }

    #[test]
    fn test_deref() {
        let (seq, _) = Sequence::normalize("test", b"ACGU").unwrap();
        let slice: &[Base] = &seq; // Deref coercion

        assert_eq!(slice.len(), 4);
        assert_eq!(slice[0], Base::A);
    }

    #[test]
    fn test_iter() {
        let (seq, _) = Sequence::normalize("test", b"ACG").unwrap();
        let bases: Vec<Base> = seq.iter().copied().collect();

        assert_eq!(bases, vec![Base::A, Base::C, Base::G]);
    }
}
