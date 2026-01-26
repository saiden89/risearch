use crate::seq::Sequence;

/// Compute RNA reverse complement for a normalized Sequence.
#[inline]
pub fn reverse_complement_rna(data: &Sequence) -> Sequence {
    data.reverse_complement()
}

/// Compute DNA reverse complement for a normalized Sequence.
///
/// Note: sequences are normalized to RNA (U instead of T), so this mirrors RNA behavior.
#[inline]
pub fn reverse_complement_dna(data: &Sequence) -> Sequence {
    data.reverse_complement()
}
