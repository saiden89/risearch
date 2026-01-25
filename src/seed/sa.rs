use libsais::SuffixArrayConstruction;

use crate::types::COMPLEMENT;

/// Compute the complement of an RNA sequence (A<->U, C<->G).
pub fn complement_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| COMPLEMENT[b as usize]).collect()
}

/// Build a suffix array for the given text. Returns the SA as u32 values.
pub fn build_suffix_array(seq: &[u8]) -> Vec<u32> {
    SuffixArrayConstruction::for_text(seq)
        .in_owned_buffer()
        .single_threaded()
        .run()
        .expect("SA construction should not fail for valid sequences")
        .into_vec()
        .into_iter()
        .map(|x: i64| x as u32)
        .collect()
}
