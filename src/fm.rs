use crate::io::read_fasta_sequences;
use genedex::{FmIndexConfig, alphabet};
use needletail::Sequence;
use savefile::prelude::*;
use std::{error::Error, path::Path};

#[derive(Savefile)]
pub struct FmIndexFile {
    pub fm_index: FmIndex<i64>,
    pub sequences: Vec<Vec<u8>>,
    pub ids: Vec<String>,
}

use genedex::FmIndex;
use serde::{Deserialize, Serialize};

pub fn create_fm_index(
    input_file: impl AsRef<Path>,
    output_file: impl AsRef<Path>,
) -> Result<(), Box<dyn Error>> {
    let sequences = read_fasta_sequences(input_file)?;

    // Build forward index
    let seq_refs: Vec<&[u8]> = sequences.iter().map(|(_, seq)| seq.as_slice()).collect();
    let dna_n_alphabet = alphabet::ascii_dna_with_n();
    let fwd_index = FmIndexConfig::<i64>::new().construct_index(&seq_refs, dna_n_alphabet.clone());

    // Unzip sequences into ids and raw sequences, avoiding clone
    let (ids, raw_sequences): (Vec<String>, Vec<Vec<u8>>) = sequences.into_iter().unzip();

    let index_file = FmIndexFile {
        fm_index: fwd_index,
        sequences: raw_sequences,
        ids,
    };

    // Save index using savefile
    let output_path = output_file.as_ref();
    save_file(output_path, 0, &index_file)?;

    Ok(())
}

pub fn load_fm_index(input_file: impl AsRef<Path>) -> Result<FmIndexFile, Box<dyn Error>> {
    let output_path = input_file.as_ref();
    let index: FmIndexFile = load_file(output_path, 0)?;
    Ok(index)
}
