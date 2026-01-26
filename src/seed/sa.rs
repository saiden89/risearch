use libsais::SuffixArrayConstruction;

use crate::types::Base;

/// Compute the complement of a sequence.
///
/// Takes a slice of Base enums and returns a Vec<Base> with complements.
/// A<->U, G<->C, others unchanged.
pub fn complement_sequence(seq: &[Base]) -> Vec<Base> {
    seq.iter().map(|&b| b.complement()).collect()
}

/// Build a suffix array for the given sequence.
///
/// Converts the sequence to discriminant bytes (0-5) for libsais,
/// which sorts by these values: Gap=0 < A=1 < G=2 < C=3 < U=4 < N=5.
///
/// Returns the SA as u32 values.
pub fn build_suffix_array(seq: &[Base]) -> Vec<u32> {
    // Convert Base slice to discriminant bytes for libsais
    let bytes: Vec<u8> = seq.iter().map(|&b| b as u8).collect();

    SuffixArrayConstruction::for_text(&bytes)
        .in_owned_buffer()
        .single_threaded()
        .run()
        .expect("SA construction should not fail for valid sequences")
        .into_vec()
        .into_iter()
        .map(|x: i64| x as u32)
        .collect()
}
