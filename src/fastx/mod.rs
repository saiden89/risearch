use anyhow::{Context, Result};
use needletail::parse_fastx_file;
use std::path::Path;

/// Type alias for FASTA records to reduce type complexity
pub type FastaRecords = Vec<(String, Vec<u8>)>;

/// Read FASTA/FASTQ and return Vec of (id, sequence) tuples.
pub fn read_fasta_sequences(filename: impl AsRef<Path>) -> Result<FastaRecords> {
    let filename_ref = filename.as_ref();

    // needletail reports an error for an empty file; for CLI UX it's nicer to
    // treat it as "no records" and let callers decide how to handle that.
    let md = std::fs::metadata(filename_ref).with_context(|| {
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

    while let Some(record) = reader.next() {
        let rec = record.with_context(|| {
            format!(
                "Failed to parse FASTA/FASTQ record from {}",
                filename_ref.display()
            )
        })?;
        let id = String::from_utf8_lossy(rec.id()).into_owned();
        let seq = rec.seq().to_vec();
        records.push((id, seq));
    }

    Ok(records)
}
