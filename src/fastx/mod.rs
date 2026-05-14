use anyhow::{bail, Context, Result};
use needletail::parse_fastx_file;
use std::collections::HashSet;
use std::path::Path;

/// Type alias for FASTA records to reduce type complexity
pub type FastaRecords = Vec<(String, Vec<u8>)>;

/// Read FASTA/FASTQ and return Vec of (id, sequence) tuples.
/// Validates that records have non-empty IDs and that there are no duplicates.
pub fn read_and_validate_fasta(filename: impl AsRef<Path>) -> Result<FastaRecords> {
    let filename_ref = filename.as_ref();

    let md = fs_err::metadata(filename_ref).with_context(|| {
        format!(
            "Failed to access FASTA/FASTQ file: {}",
            filename_ref.display()
        )
    })?;
    if md.is_file() && md.len() == 0 {
        return Ok(Vec::new());
    }

    let mut reader = parse_fastx_file(filename_ref).with_context(|| {
        format!(
            "Failed to open FASTA/FASTQ file: {}",
            filename_ref.display()
        )
    })?;

    let mut records = Vec::new();
    let mut seen = HashSet::new();

    while let Some(record) = reader.next() {
        let rec = record.with_context(|| {
            format!(
                "Failed to parse FASTA/FASTQ record from {}",
                filename_ref.display()
            )
        })?;

        let id = String::from_utf8_lossy(rec.id()).into_owned();
        if id.trim().is_empty() {
            bail!(
                "Encountered empty FASTA record id in {}",
                filename_ref.display()
            );
        }
        if !seen.insert(id.clone()) {
            bail!(
                "Duplicate FASTA record id '{}' in {}",
                id,
                filename_ref.display()
            );
        }

        let seq = rec.seq().to_vec();
        records.push((id, seq));
    }

    Ok(records)
}

/// Normalizes a raw sequence and handles logging for gaps/N-conversions.
/// Returns `Ok(None)` if the sequence is entirely empty after normalization.
pub fn normalize_record(id: &str, raw_seq: &[u8]) -> Result<Option<crate::seq::Sequence>> {
    let (sequence, stats) = crate::seq::Sequence::normalize(id, raw_seq)
        .with_context(|| format!("Failed to normalize sequence '{}'", id))?;

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
