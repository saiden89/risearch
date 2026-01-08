use std::{collections::HashSet, path::Path};

use anyhow::{Context, Result, anyhow, bail};
use libsais::SuffixArrayConstruction;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::io::read_fasta_sequences;
use needletail::Sequence;

/// Structure representing the suffix array index for a single sequence
#[derive(Serialize, Deserialize)]
pub struct SequenceIndex {
    /// Sequence identifier.
    pub name: String,
    /// Suffix array for the forward strand.
    pub forward_sa: Vec<i64>,
    /// Suffix array for the reverse strand.
    pub reverse_sa: Vec<i64>,
    /// Normalized RNA sequence (lowercase, gaps removed, ambiguous bases as 'n')
    pub sequence: Vec<u8>,
    /// Pre-computed reverse complement (avoids allocation on every access)
    pub sequence_rc: Vec<u8>,
}

/// Structure representing the entire index file containing multiple sequences.
#[derive(Serialize, Deserialize)]
pub struct SaIndexFile {
    /// List of sequence indices in the index file.
    pub sequences: Vec<SequenceIndex>,
}

#[derive(Default, Debug, Clone, Copy)]
struct NormalizationStats {
    removed_gaps: usize,
    converted_to_n: usize,
}

fn validate_readable_file(path: &Path) -> Result<()> {
    let md = std::fs::metadata(path)
        .with_context(|| format!("Failed to access input path: {}", path.display()))?;
    if !md.is_file() {
        bail!("Input path is not a file: {}", path.display());
    }
    Ok(())
}

fn validate_output_path(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        let md = std::fs::metadata(parent)
            .with_context(|| format!("Output directory does not exist: {}", parent.display()))?;
        if !md.is_dir() {
            bail!("Output parent is not a directory: {}", parent.display());
        }
    }
    Ok(())
}

// TODO: this should adheto to Rust patterns, deuplicated code
fn normalize_rna_sequence(id: &str, seq: &[u8]) -> Result<(Vec<u8>, NormalizationStats)> {
    let mut out = Vec::with_capacity(seq.len());
    let mut stats = NormalizationStats::default();

    for &b in seq {
        let c = b.to_ascii_lowercase();
        match c {
            b'a' | b'c' | b'g' => out.push(c),
            b'u' | b't' => out.push(b't'),
            b'n' => out.push(b'n'),
            b'-' | b'.' => {
                stats.removed_gaps += 1;
            }
            _ if c.is_ascii_alphabetic() => {
                // Preserve behavior similar to legacy RIsearch2: map ambiguous bases to N.
                out.push(b'n');
                stats.converted_to_n += 1;
            }
            _ => {
                bail!(
                    "Invalid character in sequence '{}': byte=0x{:02X} ('{}')",
                    id,
                    b,
                    b as char
                );
            }
        }
    }

    Ok((out, stats))
}

// TODO: This probablu belongs in some trait
pub fn create_suffix_array(
    input_file: impl AsRef<Path>,
    output_file: impl AsRef<Path>,
) -> Result<()> {
    validate_readable_file(input_file.as_ref())?;
    validate_output_path(output_file.as_ref())?;

    let index = process_sequences(&input_file).context("Failed to process input sequences")?;
    write_index_file(&index, &output_file).context("Failed to write index file")?;

    Ok(())
}

pub fn process_sequences(filename: impl AsRef<Path>) -> Result<SaIndexFile> {
    validate_readable_file(filename.as_ref())?;

    let sequences = read_fasta_sequences(&filename).context("Failed to read FASTA sequences")?;

    if sequences.is_empty() {
        bail!(
            "No sequences found in input file: {}",
            filename.as_ref().display()
        );
    }

    let mut seen = HashSet::with_capacity(sequences.len());
    for (id, _) in &sequences {
        if id.trim().is_empty() {
            bail!(
                "Encountered empty FASTA record id in {}",
                filename.as_ref().display()
            );
        }
        if !seen.insert(id.clone()) {
            bail!(
                "Duplicate FASTA record id '{}' in {}",
                id,
                filename.as_ref().display()
            );
        }
    }

    // Compute suffix arrays in parallel and collect results
    let maybe_indices: Vec<Option<SequenceIndex>> = sequences
        .into_par_iter()
        .map(|(id, seq)| -> Result<Option<SequenceIndex>> {
            let (seq_norm, stats) = normalize_rna_sequence(&id, &seq)
                .with_context(|| format!("Failed to normalize sequence '{}'", id))?;

            if seq_norm.is_empty() {
                // Header-only records or all-gaps should not crash indexing.
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

            let seq_rc = seq_norm.as_slice().reverse_complement();
            let sa_fwd = SuffixArrayConstruction::for_text(&seq_norm)
                .in_owned_buffer()
                .single_threaded()
                .run()
                .map_err(|e| anyhow!("Suffix array construction failed for '{}': {e:?}", id))?
                .into_vec();
            let sa_rev = SuffixArrayConstruction::for_text(&seq_rc)
                .in_owned_buffer()
                .single_threaded()
                .run()
                .map_err(|e| anyhow!("Suffix array construction failed for '{}' (revcomp): {e:?}", id))?
                .into_vec();

            Ok(Some(SequenceIndex {
                name: id,
                forward_sa: sa_fwd,
                reverse_sa: sa_rev,
                sequence: seq_norm,
                sequence_rc: seq_rc,
            }))
        })
        .collect::<Result<Vec<_>>>()?;

    let sequence_indices: Vec<SequenceIndex> = maybe_indices.into_iter().flatten().collect();

    if sequence_indices.is_empty() {
        bail!(
            "All sequences were empty after normalization in {}",
            filename.as_ref().display()
        );
    }

    Ok(SaIndexFile {
        sequences: sequence_indices,
    })
}

// TODO: trait for indices

pub fn write_index_file(index: &SaIndexFile, output_file: impl AsRef<Path>) -> Result<()> {
    let encoded = bincode::serialize(index).context("Failed to serialize index")?;
    std::fs::write(&output_file, encoded).context("Failed to write index file")?;
    Ok(())
}

// TODO: trait for indices
pub fn load_index_file(input_file: impl AsRef<Path>) -> Result<SaIndexFile> {
    let data = std::fs::read(&input_file).context("Failed to read index file")?;
    let index: SaIndexFile = bincode::deserialize(&data).context("Failed to deserialize index")?;
    Ok(index)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::{NamedTempFile, TempDir};

    /// Helper to create a temporary FASTA file with given content
    fn create_temp_fasta(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().expect("Failed to create temp file");
        file.write_all(content.as_bytes())
            .expect("Failed to write to temp file");
        file.flush().expect("Failed to flush temp file");
        file
    }

    #[test]
    fn test_process_sequences_rejects_empty_fasta() {
        let file = create_temp_fasta("");
        let err = process_sequences(file.path())
            .err()
            .expect("expected an error for empty FASTA");
        assert!(
            err.to_string().contains("No sequences found"),
            "unexpected error: {}",
            err
        );
    }

    #[test]
    fn test_process_sequences_rejects_duplicate_ids() {
        let file = create_temp_fasta(">dup\nACGU\n>dup\nACGU\n");
        let err = process_sequences(file.path())
            .err()
            .expect("expected an error for duplicate ids");
        assert!(
            err.to_string().contains("Duplicate FASTA record id"),
            "unexpected error: {}",
            err
        );
    }

    #[test]
    fn test_process_sequences_normalizes_t_and_ambiguous_to_n() {
        let file = create_temp_fasta(">s1\nACGTYR\n");
        let idx = process_sequences(file.path()).expect("indexing failed");
        assert_eq!(idx.sequences.len(), 1);
        assert_eq!(idx.sequences[0].name, "s1");
        assert_eq!(idx.sequences[0].sequence, b"acgtnn".to_vec());
    }

    #[test]
    fn test_process_sequences_skips_empty_after_normalization() {
        // First record is gaps only and becomes empty; second is valid.
        let file = create_temp_fasta(">empty\n--..\n>ok\nACGU\n");
        let idx = process_sequences(file.path()).expect("indexing failed");
        assert_eq!(idx.sequences.len(), 1);
        assert_eq!(idx.sequences[0].name, "ok");
        assert_eq!(idx.sequences[0].sequence, b"acgt".to_vec()); // U normalized to T
    }

    #[test]
    fn test_process_sequences_rejects_non_alphabetic() {
        let file = create_temp_fasta(">bad\nACG*U\n");
        let err = process_sequences(file.path())
            .err()
            .expect("expected an error for non-alphabetic characters");
        assert!(
            err.to_string()
                .contains("Failed to normalize sequence 'bad'"),
            "unexpected error: {}",
            err
        );
    }

    // ==================== SequenceIndex Tests ====================

    #[test]
    fn test_sequence_index_struct_fields() {
        let index = SequenceIndex {
            name: "test_seq".to_string(),
            forward_sa: vec![0, 1, 2],
            reverse_sa: vec![2, 1, 0],
            sequence: b"acgt".to_vec(),
            sequence_rc: b"acgt".to_vec(), // RC of acgt is acgt
        };

        assert_eq!(index.name, "test_seq");
        assert_eq!(index.forward_sa, vec![0, 1, 2]);
        assert_eq!(index.reverse_sa, vec![2, 1, 0]);
        assert_eq!(index.sequence, b"acgt".to_vec());
    }

    #[test]
    fn test_sequence_index_serialization() {
        let index = SequenceIndex {
            name: "seq1".to_string(),
            forward_sa: vec![3, 0, 1, 2],
            reverse_sa: vec![0, 3, 2, 1],
            sequence: b"acgt".to_vec(),
            sequence_rc: b"acgt".to_vec(),
        };

        let encoded = bincode::serialize(&index).expect("Serialization failed");
        let decoded: SequenceIndex =
            bincode::deserialize(&encoded).expect("Deserialization failed");

        assert_eq!(decoded.name, index.name);
        assert_eq!(decoded.forward_sa, index.forward_sa);
        assert_eq!(decoded.reverse_sa, index.reverse_sa);
        assert_eq!(decoded.sequence, index.sequence);
    }

    // ==================== IndexFile Tests ====================

    #[test]
    fn test_index_file_empty() {
        let index_file = SaIndexFile {
            sequences: Vec::new(),
        };
        assert!(index_file.sequences.is_empty());
    }

    #[test]
    fn test_index_file_multiple_sequences() {
        let index_file = SaIndexFile {
            sequences: vec![
                SequenceIndex {
                    name: "seq1".to_string(),
                    forward_sa: vec![0],
                    reverse_sa: vec![0],
                    sequence: b"a".to_vec(),
                    sequence_rc: b"t".to_vec(), // RC of a is t
                },
                SequenceIndex {
                    name: "seq2".to_string(),
                    forward_sa: vec![0, 1],
                    reverse_sa: vec![1, 0],
                    sequence: b"ac".to_vec(),
                    sequence_rc: b"gt".to_vec(), // RC of ac is gt
                },
            ],
        };
        assert_eq!(index_file.sequences.len(), 2);
        assert_eq!(index_file.sequences[0].name, "seq1");
        assert_eq!(index_file.sequences[1].name, "seq2");
    }

    #[test]
    fn test_index_file_serialization() {
        let index_file = SaIndexFile {
            sequences: vec![SequenceIndex {
                name: "test".to_string(),
                forward_sa: vec![0, 1],
                reverse_sa: vec![1, 0],
                sequence: b"at".to_vec(),
                sequence_rc: b"at".to_vec(), // RC of at is at
            }],
        };

        let encoded = bincode::serialize(&index_file).expect("Serialization failed");
        let decoded: SaIndexFile = bincode::deserialize(&encoded).expect("Deserialization failed");

        assert_eq!(decoded.sequences.len(), 1);
        assert_eq!(decoded.sequences[0].name, "test");
    }

    // ==================== process_sequences Tests ====================

    #[test]
    fn test_process_sequences_single_sequence() {
        let fasta_content = ">seq1\nACGT\n";
        let temp_file = create_temp_fasta(fasta_content);

        let result = process_sequences(temp_file.path()).expect("Processing failed");

        assert_eq!(result.sequences.len(), 1);
        assert_eq!(result.sequences[0].name, "seq1");
        // Sequence should be lowercase
        assert_eq!(result.sequences[0].sequence, b"acgt".to_vec());
        // Suffix arrays should have same length as sequence
        assert_eq!(result.sequences[0].forward_sa.len(), 4);
        assert_eq!(result.sequences[0].reverse_sa.len(), 4);
    }

    #[test]
    fn test_process_sequences_multiple_sequences() {
        let fasta_content = ">seq1\nACGT\n>seq2\nTGCA\n>seq3\nAAAA\n";
        let temp_file = create_temp_fasta(fasta_content);

        let result = process_sequences(temp_file.path()).expect("Processing failed");

        assert_eq!(result.sequences.len(), 3);
        // Note: parallel processing may change order, so check by finding sequences
        let names: Vec<&str> = result.sequences.iter().map(|s| s.name.as_str()).collect();
        assert!(names.contains(&"seq1"));
        assert!(names.contains(&"seq2"));
        assert!(names.contains(&"seq3"));
    }

    #[test]
    fn test_process_sequences_lowercase_conversion() {
        let fasta_content = ">test\nAcGtAcGt\n";
        let temp_file = create_temp_fasta(fasta_content);

        let result = process_sequences(temp_file.path()).expect("Processing failed");

        assert_eq!(result.sequences[0].sequence, b"acgtacgt".to_vec());
    }

    #[test]
    fn test_process_sequences_already_lowercase() {
        let fasta_content = ">test\nacgt\n";
        let temp_file = create_temp_fasta(fasta_content);

        let result = process_sequences(temp_file.path()).expect("Processing failed");

        assert_eq!(result.sequences[0].sequence, b"acgt".to_vec());
    }

    #[test]
    fn test_process_sequences_multiline_sequence() {
        let fasta_content = ">seq1\nACGT\nTGCA\nAAAA\n";
        let temp_file = create_temp_fasta(fasta_content);

        let result = process_sequences(temp_file.path()).expect("Processing failed");

        assert_eq!(result.sequences.len(), 1);
        assert_eq!(result.sequences[0].sequence, b"acgttgcaaaaa".to_vec());
    }

    #[test]
    fn test_process_sequences_nonexistent_file() {
        let result = process_sequences("/nonexistent/path/file.fa");
        assert!(result.is_err());
    }

    #[test]
    fn test_process_sequences_suffix_array_validity() {
        // Test that suffix arrays are valid (contain all indices 0..n-1)
        let fasta_content = ">test\nACGT\n";
        let temp_file = create_temp_fasta(fasta_content);

        let result = process_sequences(temp_file.path()).expect("Processing failed");
        let seq_idx = &result.sequences[0];

        // Check forward SA contains all indices
        let mut forward_sorted: Vec<i64> = seq_idx.forward_sa.clone();
        forward_sorted.sort();
        assert_eq!(forward_sorted, vec![0, 1, 2, 3]);

        // Check reverse SA contains all indices
        let mut reverse_sorted: Vec<i64> = seq_idx.reverse_sa.clone();
        reverse_sorted.sort();
        assert_eq!(reverse_sorted, vec![0, 1, 2, 3]);
    }

    // ==================== write_index_file Tests ====================

    #[test]
    fn test_write_index_file_success() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let output_path = temp_dir.path().join("test.idx");

        let index = SaIndexFile {
            sequences: vec![SequenceIndex {
                name: "test".to_string(),
                forward_sa: vec![0, 1, 2],
                reverse_sa: vec![2, 1, 0],
                sequence: b"acg".to_vec(),
                sequence_rc: b"cgt".to_vec(), // RC of acg is cgt
            }],
        };

        let result = write_index_file(&index, &output_path);
        assert!(result.is_ok());
        assert!(output_path.exists());
    }

    #[test]
    fn test_write_index_file_invalid_path() {
        let index = SaIndexFile {
            sequences: Vec::new(),
        };

        let result = write_index_file(&index, "/nonexistent/directory/file.idx");
        assert!(result.is_err());
    }

    #[test]
    fn test_write_index_file_empty_index() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let output_path = temp_dir.path().join("empty.idx");

        let index = SaIndexFile {
            sequences: Vec::new(),
        };

        let result = write_index_file(&index, &output_path);
        assert!(result.is_ok());
        assert!(output_path.exists());

        // File should still be readable
        let data = std::fs::read(&output_path).expect("Failed to read file");
        assert!(!data.is_empty()); // bincode header still present
    }

    // ==================== load_index_file Tests ====================

    #[test]
    fn test_load_index_file_success() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("test.idx");

        let original = SaIndexFile {
            sequences: vec![SequenceIndex {
                name: "loaded_seq".to_string(),
                forward_sa: vec![3, 0, 1, 2],
                reverse_sa: vec![0, 3, 2, 1],
                sequence: b"acgt".to_vec(),
                sequence_rc: b"acgt".to_vec(), // RC of acgt is acgt
            }],
        };

        write_index_file(&original, &file_path).expect("Write failed");

        let loaded = load_index_file(&file_path).expect("Load failed");

        assert_eq!(loaded.sequences.len(), 1);
        assert_eq!(loaded.sequences[0].name, "loaded_seq");
        assert_eq!(loaded.sequences[0].forward_sa, vec![3, 0, 1, 2]);
        assert_eq!(loaded.sequences[0].reverse_sa, vec![0, 3, 2, 1]);
        assert_eq!(loaded.sequences[0].sequence, b"acgt".to_vec());
    }

    #[test]
    fn test_load_index_file_nonexistent() {
        let result = load_index_file("/nonexistent/path/file.idx");
        assert!(result.is_err());
    }

    #[test]
    fn test_load_index_file_invalid_data() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("invalid.idx");

        // Write invalid binary data
        std::fs::write(&file_path, b"not valid bincode data").expect("Write failed");

        let result = load_index_file(&file_path);
        assert!(result.is_err());
    }

    #[test]
    fn test_load_index_file_empty_file() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("empty.idx");

        std::fs::write(&file_path, b"").expect("Write failed");

        let result = load_index_file(&file_path);
        assert!(result.is_err());
    }

    // ==================== create_suffix_array Tests ====================

    #[test]
    fn test_create_suffix_array_success() {
        let fasta_content = ">seq1\nACGTACGT\n";
        let temp_fasta = create_temp_fasta(fasta_content);
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let output_path = temp_dir.path().join("output.idx");

        let result = create_suffix_array(temp_fasta.path(), &output_path);
        assert!(result.is_ok());
        assert!(output_path.exists());

        // Verify the created file can be loaded
        let loaded = load_index_file(&output_path).expect("Failed to load created index");
        assert_eq!(loaded.sequences.len(), 1);
        assert_eq!(loaded.sequences[0].name, "seq1");
    }

    #[test]
    fn test_create_suffix_array_invalid_input() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let output_path = temp_dir.path().join("output.idx");

        let result = create_suffix_array("/nonexistent/input.fa", &output_path);
        assert!(result.is_err());
    }

    #[test]
    fn test_create_suffix_array_invalid_output() {
        let fasta_content = ">seq1\nACGT\n";
        let temp_fasta = create_temp_fasta(fasta_content);

        let result = create_suffix_array(temp_fasta.path(), "/nonexistent/dir/output.idx");
        assert!(result.is_err());
    }

    #[test]
    fn test_create_suffix_array_multiple_sequences() {
        let fasta_content = ">seq1\nACGT\n>seq2\nTGCA\n";
        let temp_fasta = create_temp_fasta(fasta_content);
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let output_path = temp_dir.path().join("multi.idx");

        let result = create_suffix_array(temp_fasta.path(), &output_path);
        assert!(result.is_ok());

        let loaded = load_index_file(&output_path).expect("Failed to load");
        assert_eq!(loaded.sequences.len(), 2);
    }

    // ==================== Round-trip Tests ====================

    #[test]
    fn test_full_roundtrip() {
        let fasta_content = ">myseq description here\nACGTACGTAAAACCCCGGGGTTTT\n";
        let temp_fasta = create_temp_fasta(fasta_content);
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let index_path = temp_dir.path().join("roundtrip.idx");

        // Create index from FASTA
        create_suffix_array(temp_fasta.path(), &index_path).expect("Create failed");

        // Load it back
        let loaded = load_index_file(&index_path).expect("Load failed");

        // Verify content
        assert_eq!(loaded.sequences.len(), 1);
        assert_eq!(loaded.sequences[0].name, "myseq description here");
        assert_eq!(
            loaded.sequences[0].sequence,
            b"acgtacgtaaaaccccggggtttt".to_vec()
        );
        assert_eq!(loaded.sequences[0].forward_sa.len(), 24);
        assert_eq!(loaded.sequences[0].reverse_sa.len(), 24);
    }

    #[test]
    fn test_write_load_roundtrip_preserves_data() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("roundtrip.idx");

        let original = SaIndexFile {
            sequences: vec![
                SequenceIndex {
                    name: "first".to_string(),
                    forward_sa: vec![5, 4, 3, 2, 1, 0],
                    reverse_sa: vec![0, 1, 2, 3, 4, 5],
                    sequence: b"aaaaaa".to_vec(),
                    sequence_rc: b"tttttt".to_vec(), // RC of all a's is all t's
                },
                SequenceIndex {
                    name: "second".to_string(),
                    forward_sa: vec![0],
                    reverse_sa: vec![0],
                    sequence: b"c".to_vec(),
                    sequence_rc: b"g".to_vec(), // RC of c is g
                },
            ],
        };

        write_index_file(&original, &file_path).expect("Write failed");
        let loaded = load_index_file(&file_path).expect("Load failed");

        assert_eq!(loaded.sequences.len(), original.sequences.len());
        for (orig, load) in original.sequences.iter().zip(loaded.sequences.iter()) {
            assert_eq!(orig.name, load.name);
            assert_eq!(orig.forward_sa, load.forward_sa);
            assert_eq!(orig.reverse_sa, load.reverse_sa);
            assert_eq!(orig.sequence, load.sequence);
        }
    }
}
