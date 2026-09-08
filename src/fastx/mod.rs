//! FASTA/FASTQ reading and normalization into [`Sequence`](crate::Sequence).

use needletail::parse_fastx_file;
use std::collections::HashSet;
use std::path::Path;

use crate::error::{Error, Result};

/// Type alias for FASTA records to reduce type complexity
pub(crate) type FastaRecords = Vec<(String, Vec<u8>)>;

/// Read FASTA/FASTQ and return Vec of (id, sequence) tuples.
/// Validates that records have non-empty IDs and that there are no duplicates.
pub(crate) fn read_and_validate_fasta(filename: impl AsRef<Path>) -> Result<FastaRecords> {
    let filename_ref = filename.as_ref();

    let md = fs_err::metadata(filename_ref)?;
    if !md.is_file() {
        return Err(Error::Input(format!(
            "Input path is not a file: {}",
            filename_ref.display()
        )));
    }
    if md.len() == 0 {
        return Ok(Vec::new());
    }

    let mut reader = parse_fastx_file(filename_ref).map_err(|err| {
        Error::Input(format!(
            "Failed to open FASTA/FASTQ file {}: {err}",
            filename_ref.display()
        ))
    })?;

    let mut records = Vec::new();
    let mut seen = HashSet::new();

    while let Some(record) = reader.next() {
        let rec = record.map_err(|err| {
            Error::Input(format!(
                "Failed to parse FASTA/FASTQ record from {}: {err}",
                filename_ref.display()
            ))
        })?;

        let id = String::from_utf8_lossy(rec.id()).into_owned();
        if id.trim().is_empty() {
            return Err(Error::Input(format!(
                "Encountered empty FASTA record id in {}",
                filename_ref.display()
            )));
        }
        if !seen.insert(id.clone()) {
            return Err(Error::Input(format!(
                "Duplicate FASTA record id '{}' in {}",
                id,
                filename_ref.display()
            )));
        }

        let seq = rec.seq().to_vec();
        records.push((id, seq));
    }

    Ok(records)
}

/// Read a FASTA/FASTQ file into normalized sequences.
///
/// Records that normalize to nothing are dropped; an input where every record
/// does is an error.
pub fn read_sequences(filename: impl AsRef<Path>) -> Result<Vec<(String, crate::seq::Sequence)>> {
    let filename = filename.as_ref();
    let records = read_and_validate_fasta(filename)?;

    let mut sequences = Vec::with_capacity(records.len());
    for (id, raw_seq) in records {
        if let Some(sequence) = normalize_record(&id, &raw_seq)? {
            sequences.push((id, sequence));
        }
    }

    if sequences.is_empty() {
        return Err(Error::Input(format!(
            "All sequences were empty after normalization in {}",
            filename.display()
        )));
    }

    Ok(sequences)
}

/// Normalizes a raw sequence and handles logging for gaps/N-conversions.
/// Returns `Ok(None)` if the sequence is entirely empty after normalization.
pub(crate) fn normalize_record(id: &str, raw_seq: &[u8]) -> Result<Option<crate::seq::Sequence>> {
    let (sequence, stats) = crate::seq::Sequence::normalize(id, raw_seq)?;

    if sequence.is_empty() {
        log::warn!(
            "Skipping empty sequence after normalization: '{}' (removed_gaps={}, converted_to_n={})",
            id,
            stats.removed_gaps,
            stats.converted_to_n
        );
        return Ok(None);
    }

    if stats.removed_gaps > 0 || stats.converted_to_n > 0 {
        log::debug!(
            "Normalized sequence '{}': removed_gaps={}, converted_to_n={}",
            id,
            stats.removed_gaps,
            stats.converted_to_n
        );
    }

    Ok(Some(sequence))
}
