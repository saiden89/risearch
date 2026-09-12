//! Sequence newtype wrapping Vec<Base>
//!
//! This module provides the core `Sequence` type that wraps `Vec<Base>`,
//! ensuring sequences are normalized once at the input boundary and then
//! work with Base enums throughout the codebase.

use crate::error::Result;
use crate::types::Base;
use std::ops::{Deref, Index, Range, RangeFrom, RangeFull, RangeTo};

use super::normalize::{normalize_rna_sequence, NormalizationStats};

/// A normalized RNA sequence stored as Vec<Base>.
///
/// Sequences are normalized once at input (FASTA parsing) and stored as Base enums.
/// This eliminates repeated ASCII→Base conversions during search and alignment.
#[derive(Clone, Debug, PartialEq, Eq)]
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
    use rstest::rstest;

    #[test]
    fn normalize_folds_case_and_maps_t_to_u() {
        let (seq, stats) = Sequence::normalize("test", b"aCgUtN").unwrap();
        assert_eq!(
            &*seq,
            &[Base::A, Base::C, Base::G, Base::U, Base::U, Base::N]
        );
        assert_eq!(stats.removed_gaps, 0);
        assert_eq!(stats.converted_to_n, 0);
    }

    #[rstest]
    #[case::gaps_between_bases(b"AC-.GU", &[Base::A, Base::C, Base::G, Base::U])]
    #[case::only_gaps(b"-.", &[])]
    fn normalize_strips_gaps_and_counts_them(#[case] raw: &[u8], #[case] expected: &[Base]) {
        let (seq, stats) = Sequence::normalize("test", raw).unwrap();
        assert_eq!(&*seq, expected);
        assert_eq!(stats.removed_gaps, 2);
    }

    // Alphabetic ambiguity is deliberately accepted and accounted for;
    // punctuation/non-ASCII bytes are rejected by a separate contract.
    #[rstest]
    #[case::single_code(b"ACRGU", &[Base::A, Base::C, Base::N, Base::G, Base::U], 1)]
    #[case::every_code(b"RYSWKMBDHVXZJryswkmbdhvxzj", &[Base::N; 26], 26)]
    fn normalize_converts_ambiguous_to_n(
        #[case] raw: &[u8],
        #[case] expected: &[Base],
        #[case] converted: usize,
    ) {
        let (seq, stats) = Sequence::normalize("ambiguities", raw).unwrap();
        assert_eq!(&*seq, expected);
        assert_eq!(stats.converted_to_n, converted);
    }

    #[rstest]
    #[case::punctuation(b'?')]
    #[case::digit(b'1')]
    #[case::nul(0)]
    #[case::high_byte(0xff)]
    fn normalize_rejects_non_base_bytes_with_record_context(#[case] byte: u8) {
        let error = Sequence::normalize("bad-record", &[byte]).unwrap_err();
        assert!(matches!(error, crate::Error::Input(_)));
        let message = error.to_string();
        assert!(message.contains("bad-record"));
        assert!(message.contains(&format!("0x{byte:02X}")));
    }

    #[test]
    fn deref_coercion_to_slice() {
        let (seq, _) = Sequence::normalize("test", b"ACGU").unwrap();
        let slice: &[Base] = &seq; // Deref coercion

        assert_eq!(slice.len(), 4);
        assert_eq!(slice[0], Base::A);
    }

    #[test]
    fn indexing_returns_the_underlying_bases() {
        let (seq, _) = Sequence::normalize("test", b"ACGU").unwrap();

        assert_eq!(seq.len(), 4);
        assert!(!seq.is_empty());
        assert!(Sequence::from(Vec::new()).is_empty());
        assert_eq!(seq[0], Base::A);
        assert_eq!(seq[1], Base::C);
        assert_eq!(seq[2], Base::G);
        assert_eq!(seq[3], Base::U);
        assert_eq!(&seq[1..3], &[Base::C, Base::G]);
        assert_eq!(&seq[2..], &[Base::G, Base::U]);
        assert_eq!(&seq[..2], &[Base::A, Base::C]);
        assert_eq!(&seq[..], &[Base::A, Base::C, Base::G, Base::U]);
    }

    #[test]
    fn iter_yields_all_bases() {
        let (seq, _) = Sequence::normalize("test", b"ACG").unwrap();
        let bases: Vec<Base> = seq.iter().copied().collect();

        assert_eq!(bases, vec![Base::A, Base::C, Base::G]);
    }
}
