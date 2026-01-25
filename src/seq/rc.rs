use crate::types::{RC_DNA_TABLE, RC_RNA_TABLE};

/// Compute RNA reverse complement using LUT (single lookup per base).
/// Returns uppercase: A↔U, G↔C.
#[inline]
pub fn reverse_complement_rna(data: &[u8]) -> Vec<u8> {
    data.iter()
        .rev()
        .map(|&b| RC_RNA_TABLE[b as usize])
        .collect()
}

/// Compute DNA reverse complement using LUT (single lookup per base).
/// Returns lowercase: A↔T, G↔C.
#[inline]
pub fn reverse_complement_dna(data: &[u8]) -> Vec<u8> {
    data.iter()
        .rev()
        .map(|&b| RC_DNA_TABLE[b as usize])
        .collect()
}
