use std::{collections::HashSet, path::Path};

use anyhow::{Context, Result, anyhow, bail};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::fastx::read_fasta_sequences;
use crate::registry::TargetRegistry;
use crate::sa::SuffixArray;
use crate::seq::Sequence;

use super::io::{validate_output_path, validate_readable_file, write_index_file};

/// Structure representing the suffix array index for a single sequence
#[derive(Serialize, Deserialize)]
pub struct SequenceIndex {
    /// Sequence identifier.
    pub name: String,
    /// Suffix array for the forward strand
    pub forward_sa: SuffixArray,
    /// Suffix array for the reverse strand
    pub reverse_sa: SuffixArray,
    /// Forward sequence
    pub sequence: Sequence,
    /// Reverse complement sequence
    pub sequence_rc: Sequence,
}

pub fn create_suffix_array(
    input_file: impl AsRef<Path>,
    output_file: impl AsRef<Path>,
) -> Result<()> {
    validate_readable_file(input_file.as_ref())?;
    validate_output_path(output_file.as_ref())?;

    let sa = process_sequences(&input_file).context("Failed to process input sequences")?;
    write_index_file(&sa, &output_file).context("Failed to write index file")?;

    Ok(())
}

pub fn process_sequences(filename: impl AsRef<Path>) -> Result<TargetRegistry> {
    validate_readable_file(filename.as_ref())?;

    let sequences = read_fasta_sequences(&filename).context("Failed to read FASTA sequences")?;
    // TODO decouple processing from validation
    if sequences.is_empty() {
        bail!(
            "No sequences found in input file: {}",
            filename.as_ref().display()
        );
    }

    // TODO: decouple validation from processing
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
            let (seq_norm, stats) = Sequence::normalize(&id, &seq)
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

            // TODO: Add feature flag for u64 SA to support sequences > 4.2B bases (e.g. some plant genomes)
            if seq_norm.len() > u32::MAX as usize {
                bail!(
                    "Sequence '{}' too long for u32 suffix array ({} bases > {} max). \
                     Consider chunking the sequence or using a future u64-enabled build.",
                    id,
                    seq_norm.len(),
                    u32::MAX
                );
            }

            let seq_rc = seq_norm.reverse_complement();

            let sa_fwd = SuffixArray::try_build(&seq_norm.to_bytes())
                .map_err(|e| anyhow!("Suffix array construction failed for '{}': {e:?}", id))?;
            let sa_rev = SuffixArray::try_build(&seq_rc.to_bytes())
                .map_err(|e| anyhow!("Suffix array construction failed for '{}' (revcomp): {e:?}", id))?;

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

    Ok(TargetRegistry::new(sequence_indices))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::io::load_index_file;
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
        assert_eq!(idx.entries().len(), 1);
        assert_eq!(idx.entries()[0].name, "s1");
        // Compare with expected sequence using Sequence::normalize
        let (expected, _) = Sequence::normalize("test", b"ACGTNN").unwrap();
        assert_eq!(idx.entries()[0].sequence, expected);
    }

    #[test]
    fn test_process_sequences_skips_empty_after_normalization() {
        // First record is gaps only and becomes empty; second is valid.
        let file = create_temp_fasta(">empty\n--..\n>ok\nACGU\n");
        let idx = process_sequences(file.path()).expect("indexing failed");
        assert_eq!(idx.entries().len(), 1);
        assert_eq!(idx.entries()[0].name, "ok");
        // Compare with expected sequence using Sequence::normalize
        let (expected, _) = Sequence::normalize("test", b"ACGU").unwrap();
        assert_eq!(idx.entries()[0].sequence, expected);
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
        let (seq, _) = Sequence::normalize("test", b"ACGT").unwrap();
        let (seq_rc, _) = Sequence::normalize("test", b"ACGT").unwrap();

        let index = SequenceIndex {
            name: "test_seq".to_string(),
            forward_sa: SuffixArray::from(vec![0, 1, 2]),
            reverse_sa: SuffixArray::from(vec![2, 1, 0]),
            sequence: seq.clone(),
            sequence_rc: seq_rc,
        };

        assert_eq!(index.name, "test_seq");
        assert_eq!(&index.forward_sa[..], &[0, 1, 2]);
        assert_eq!(&index.reverse_sa[..], &[2, 1, 0]);
        assert_eq!(index.sequence, seq);
    }

    #[test]
    fn test_sequence_index_serialization() {
        let (seq, _) = Sequence::normalize("test", b"ACGT").unwrap();
        let (seq_rc, _) = Sequence::normalize("test", b"ACGT").unwrap();

        let index = SequenceIndex {
            name: "seq1".to_string(),
            forward_sa: SuffixArray::from(vec![3, 0, 1, 2]),
            reverse_sa: SuffixArray::from(vec![0, 3, 2, 1]),
            sequence: seq.clone(),
            sequence_rc: seq_rc,
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
        let index_file = TargetRegistry::new(Vec::new());
        assert!(index_file.entries().is_empty());
    }

    #[test]
    fn test_index_file_multiple_sequences() {
        let (seq1, _) = Sequence::normalize("test", b"A").unwrap();
        let (seq1_rc, _) = Sequence::normalize("test", b"U").unwrap();
        let (seq2, _) = Sequence::normalize("test", b"AC").unwrap();
        let (seq2_rc, _) = Sequence::normalize("test", b"GU").unwrap();

        let index_file = TargetRegistry::new(vec![
            SequenceIndex {
                name: "seq1".to_string(),
                forward_sa: SuffixArray::from(vec![0]),
                reverse_sa: SuffixArray::from(vec![0]),
                sequence: seq1,
                sequence_rc: seq1_rc,
            },
            SequenceIndex {
                name: "seq2".to_string(),
                forward_sa: SuffixArray::from(vec![0, 1]),
                reverse_sa: SuffixArray::from(vec![1, 0]),
                sequence: seq2,
                sequence_rc: seq2_rc,
            },
        ]);
        assert_eq!(index_file.entries().len(), 2);
        assert_eq!(index_file.entries()[0].name, "seq1");
        assert_eq!(index_file.entries()[1].name, "seq2");
    }

    #[test]
    fn test_index_file_serialization() {
        let (seq, _) = Sequence::normalize("test", b"AU").unwrap();
        let (seq_rc, _) = Sequence::normalize("test", b"AU").unwrap();

        let index_file = TargetRegistry::new(vec![SequenceIndex {
            name: "test".to_string(),
            forward_sa: SuffixArray::from(vec![0, 1]),
            reverse_sa: SuffixArray::from(vec![1, 0]),
            sequence: seq,
            sequence_rc: seq_rc,
        }]);

        let encoded = bincode::serialize(&index_file).expect("Serialization failed");
        let decoded: TargetRegistry =
            bincode::deserialize(&encoded).expect("Deserialization failed");

        assert_eq!(decoded.entries().len(), 1);
        assert_eq!(decoded.entries()[0].name, "test");
    }

    // ==================== process_sequences Tests ====================

    #[test]
    fn test_process_sequences_single_sequence() {
        let fasta_content = ">seq1\nACGT\n";
        let temp_file = create_temp_fasta(fasta_content);

        let result = process_sequences(temp_file.path()).expect("Processing failed");

        assert_eq!(result.entries().len(), 1);
        assert_eq!(result.entries()[0].name, "seq1");
        // Verify sequence length
        assert_eq!(result.entries()[0].sequence.len(), 4);
        // Suffix arrays should have same length as sequence
        assert_eq!(result.entries()[0].forward_sa.len(), 4);
        assert_eq!(result.entries()[0].reverse_sa.len(), 4);
    }

    #[test]
    fn test_process_sequences_multiple_sequences() {
        let fasta_content = ">seq1\nACGT\n>seq2\nTGCA\n>seq3\nAAAA\n";
        let temp_file = create_temp_fasta(fasta_content);

        let result = process_sequences(temp_file.path()).expect("Processing failed");

        assert_eq!(result.entries().len(), 3);
        // Note: parallel processing may change order, so check by finding sequences
        let names: Vec<&str> = result.entries().iter().map(|s| s.name.as_str()).collect();
        assert!(names.contains(&"seq1"));
        assert!(names.contains(&"seq2"));
        assert!(names.contains(&"seq3"));
    }

    #[test]
    fn test_process_sequences_lowercase_conversion() {
        let fasta_content = ">test\nAcGtAcGt\n";
        let temp_file = create_temp_fasta(fasta_content);

        let result = process_sequences(temp_file.path()).expect("Processing failed");

        // Verify sequence length (normalized)
        assert_eq!(result.entries()[0].sequence.len(), 8);
    }

    #[test]
    fn test_process_sequences_already_lowercase() {
        let fasta_content = ">test\nacgt\n";
        let temp_file = create_temp_fasta(fasta_content);

        let result = process_sequences(temp_file.path()).expect("Processing failed");

        // Verify sequence length
        assert_eq!(result.entries()[0].sequence.len(), 4);
    }

    #[test]
    fn test_process_sequences_multiline_sequence() {
        let fasta_content = ">seq1\nACGT\nTGCA\nAAAA\n";
        let temp_file = create_temp_fasta(fasta_content);

        let result = process_sequences(temp_file.path()).expect("Processing failed");

        assert_eq!(result.entries().len(), 1);
        // Verify sequence length (4+4+4=12)
        assert_eq!(result.entries()[0].sequence.len(), 12);
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
        let seq_idx = &result.entries()[0];

        // Check forward SA contains all indices
        let mut forward_sorted: Vec<u32> = seq_idx.forward_sa.iter().copied().collect();
        forward_sorted.sort();
        assert_eq!(forward_sorted, vec![0u32, 1, 2, 3]);

        // Check reverse SA contains all indices
        let mut reverse_sorted: Vec<u32> = seq_idx.reverse_sa.iter().copied().collect();
        reverse_sorted.sort();
        assert_eq!(reverse_sorted, vec![0u32, 1, 2, 3]);
    }

    // ==================== write_index_file Tests ====================

    #[test]
    fn test_write_index_file_success() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let output_path = temp_dir.path().join("test.idx");

        let (seq, _) = Sequence::normalize("test", b"ACG").unwrap();
        let (seq_rc, _) = Sequence::normalize("test", b"CGU").unwrap();

        let index = TargetRegistry::new(vec![SequenceIndex {
            name: "test".to_string(),
            forward_sa: SuffixArray::from(vec![0, 1, 2]),
            reverse_sa: SuffixArray::from(vec![2, 1, 0]),
            sequence: seq,
            sequence_rc: seq_rc,
        }]);

        let result = write_index_file(&index, &output_path);
        assert!(result.is_ok());
        assert!(output_path.exists());
    }

    #[test]
    fn test_write_index_file_invalid_path() {
        let index = TargetRegistry::new(Vec::new());

        let result = write_index_file(&index, "/nonexistent/directory/file.idx");
        assert!(result.is_err());
    }

    #[test]
    fn test_write_index_file_empty_index() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let output_path = temp_dir.path().join("empty.idx");

        let index = TargetRegistry::new(Vec::new());

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

        let (seq, _) = Sequence::normalize("test", b"ACGT").unwrap();
        let (seq_rc, _) = Sequence::normalize("test", b"ACGT").unwrap();

        let original = TargetRegistry::new(vec![SequenceIndex {
            name: "loaded_seq".to_string(),
            forward_sa: SuffixArray::from(vec![3, 0, 1, 2]),
            reverse_sa: SuffixArray::from(vec![0, 3, 2, 1]),
            sequence: seq.clone(),
            sequence_rc: seq_rc,
        }]);

        write_index_file(&original, &file_path).expect("Write failed");

        let loaded = load_index_file(&file_path).expect("Load failed");

        assert_eq!(loaded.entries().len(), 1);
        assert_eq!(loaded.entries()[0].name, "loaded_seq");
        assert_eq!(&loaded.entries()[0].forward_sa[..], &[3, 0, 1, 2]);
        assert_eq!(&loaded.entries()[0].reverse_sa[..], &[0, 3, 2, 1]);
        assert_eq!(loaded.entries()[0].sequence, seq);
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
        assert_eq!(loaded.entries().len(), 1);
        assert_eq!(loaded.entries()[0].name, "seq1");
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
        assert_eq!(loaded.entries().len(), 2);
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
        assert_eq!(loaded.entries().len(), 1);
        assert_eq!(loaded.entries()[0].name, "myseq description here");
        // Verify sequence length (24 bases)
        assert_eq!(loaded.entries()[0].sequence.len(), 24);
        assert_eq!(loaded.entries()[0].forward_sa.len(), 24);
        assert_eq!(loaded.entries()[0].reverse_sa.len(), 24);
    }

    #[test]
    fn test_write_load_roundtrip_preserves_data() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("roundtrip.idx");

        let (seq1, _) = Sequence::normalize("test", b"AAAAAA").unwrap();
        let (seq1_rc, _) = Sequence::normalize("test", b"UUUUUU").unwrap();
        let (seq2, _) = Sequence::normalize("test", b"C").unwrap();
        let (seq2_rc, _) = Sequence::normalize("test", b"G").unwrap();

        let original = TargetRegistry::new(vec![
            SequenceIndex {
                name: "first".to_string(),
                forward_sa: SuffixArray::from(vec![5, 4, 3, 2, 1, 0]),
                reverse_sa: SuffixArray::from(vec![0, 1, 2, 3, 4, 5]),
                sequence: seq1.clone(),
                sequence_rc: seq1_rc.clone(),
            },
            SequenceIndex {
                name: "second".to_string(),
                forward_sa: SuffixArray::from(vec![0]),
                reverse_sa: SuffixArray::from(vec![0]),
                sequence: seq2.clone(),
                sequence_rc: seq2_rc.clone(),
            },
        ]);

        write_index_file(&original, &file_path).expect("Write failed");
        let loaded = load_index_file(&file_path).expect("Load failed");

        assert_eq!(loaded.entries().len(), original.entries().len());
        for (orig, load) in original.entries().iter().zip(loaded.entries().iter()) {
            assert_eq!(orig.name, load.name);
            assert_eq!(orig.forward_sa, load.forward_sa);
            assert_eq!(orig.reverse_sa, load.reverse_sa);
            assert_eq!(orig.sequence, load.sequence);
        }
    }
}
