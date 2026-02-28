use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use memmap2::Mmap;
use needletail::parse_fastx_file;

use crate::index::io::{validate_output_path, validate_readable_file};
use crate::index::sa::SuffixArray;
use crate::seq::Sequence;
use crate::types::Base;

// =============================================================================
// FORMAT CONSTANTS
// =============================================================================

const FILE_MAGIC: [u8; 8] = *b"RSIDX5\0\0";

/// Fixed file header: magic[8] + target_count[4] + reserved[4]
const FILE_HEADER_BYTES: usize = 16;

/// Per-target metadata entry: name_len[4] + seq_len[4] + offset[8]
const META_ENTRY_FIXED_BYTES: usize = 16;

const DATA_ALIGN: usize = 8;

/// Number of zero-valued u64 entries appended after the real SA data.
/// These act as sentinels so `sa[suffix_pos + offset]` never goes out of
/// bounds during binary search character lookups, eliminating a branch
/// in the innermost hot loop.
pub const SA_CHAR_PADDING: usize = 256;

// =============================================================================
// PUBLIC TYPES
// =============================================================================

#[derive(Clone, Debug)]
pub struct TargetMeta {
    pub name: String,
    pub sequence_len: usize,
    /// Start offset of this target's block within the global combined_seq.
    /// Each block is: fwd_comp[seq_len] + Gap + rc_comp[seq_len] + Gap
    pub global_offset: u64,
}

pub struct TargetStore {
    path: PathBuf,
    mmap: Mmap,
    targets: Vec<TargetMeta>,
    /// Precomputed per-target offsets into global combined_seq
    offsets: Vec<u64>,
    /// Precomputed per-target forward sequence lengths
    seq_lens: Vec<u32>,
    /// Byte offset into mmap where the global combined_seq begins
    seq_data_offset: usize,
    /// Total number of Base entries in global combined_seq (including padding)
    seq_byte_count: usize,
    /// Byte offset into mmap where the global combined_sa begins
    sa_data_offset: usize,
    /// Total number of u64 entries in global combined_sa (including padding)
    sa_entry_count: usize,
    /// Number of real SA entries (excluding SA_CHAR_PADDING sentinels)
    sa_real_len: usize,
}

/// Global suffix array view for seed search — zero-copy from mmap.
pub struct GlobalView<'a> {
    pub combined_seq: &'a [Base],
    pub combined_sa: &'a [u64],
    pub sa_real_len: usize,
    /// Start offset of each target block in the global combined_seq.
    /// Length = target_count.
    pub offsets: &'a [u64],
    /// Original forward sequence length of each target.
    /// Length = target_count.
    pub seq_lens: &'a [u32],
}

/// Per-target sequence slices for output and extension.
pub struct TargetSeqs<'a> {
    pub name: &'a str,
    pub fwd_transformed: &'a [Base],
    pub rc_transformed: &'a [Base],
    pub seq_len: usize,
}

// =============================================================================
// BUILD PIPELINE
// =============================================================================

struct PendingTargetRecord {
    id: String,
    seq: Vec<u8>,
}

struct NormalizedTarget {
    id: String,
    seq_len: u32,
    fwd_comp: Vec<Base>,
    rc_comp: Vec<Base>,
}

impl TargetStore {
    /// Build a new RSIDX5 index from a FASTA file.
    ///
    /// Pipeline:
    /// 1. Parse all FASTA records, normalize, complement
    /// 2. Concatenate into global combined_seq with Gap separators
    /// 3. Build single global SA
    /// 4. Write flat binary: header → metadata → seq → SA
    pub fn build_from_fasta(input: &Path, output: &Path) -> Result<()> {
        validate_readable_file(input)?;
        validate_output_path(output)?;

        // Phase 1: Parse and normalize all records
        let mut seen = HashSet::new();
        let mut pending = Vec::new();

        let mut reader = parse_fastx_file(input)
            .with_context(|| format!("Failed to open FASTA/FASTQ file: {}", input.display()))?;

        while let Some(record) = reader.next() {
            let rec = record.with_context(|| {
                format!(
                    "Failed to parse FASTA/FASTQ record from {}",
                    input.display()
                )
            })?;

            let id = String::from_utf8_lossy(rec.id()).into_owned();
            if id.trim().is_empty() {
                bail!("Encountered empty FASTA record id in {}", input.display());
            }
            if !seen.insert(id.clone()) {
                bail!("Duplicate FASTA record id '{}' in {}", id, input.display());
            }

            pending.push(PendingTargetRecord {
                id,
                seq: rec.seq().as_ref().to_vec(),
            });
        }

        // Normalize all records (parallel)
        use rayon::prelude::*;
        let normalized: Vec<Result<Option<NormalizedTarget>>> = pending
            .into_par_iter()
            .map(normalize_target_record)
            .collect();

        let mut targets: Vec<NormalizedTarget> = Vec::new();
        for result in normalized {
            if let Some(t) = result? {
                targets.push(t);
            }
        }

        if targets.is_empty() {
            bail!(
                "All sequences were empty after normalization in {}",
                input.display()
            );
        }

        // Phase 2: Build global combined_seq and track offsets
        let total_bases: usize = targets.iter().map(|t| 2 * t.seq_len as usize + 2).sum();
        let mut combined_bases: Vec<Base> = Vec::with_capacity(total_bases);
        let mut offsets: Vec<u64> = Vec::with_capacity(targets.len());
        let mut seq_lens: Vec<u32> = Vec::with_capacity(targets.len());

        for t in &targets {
            offsets.push(combined_bases.len() as u64);
            seq_lens.push(t.seq_len);
            combined_bases.extend_from_slice(&t.fwd_comp);
            combined_bases.push(Base::Gap);
            combined_bases.extend_from_slice(&t.rc_comp);
            combined_bases.push(Base::Gap);
        }

        // Phase 3: Build single global SA directly from &[Base]
        // (avoids cloning combined_bases and the extra to_bytes() allocation)
        let combined_sa = SuffixArray::build_from_bases(&combined_bases)
            .context("Failed to build global suffix array")?;

        // Append SA_CHAR_PADDING sentinels to seq
        combined_bases.resize(combined_bases.len() + SA_CHAR_PADDING, Base::Gap);

        // Pad the SA with zero entries
        let mut padded_sa = combined_sa.into_inner();
        padded_sa.resize(padded_sa.len() + SA_CHAR_PADDING, 0u64);

        // Phase 4: Write flat binary
        let output_name = output
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "target.idx".to_string());
        let tmp_path = output.with_file_name(format!("{output_name}.tmp"));

        let file = File::create(&tmp_path)
            .with_context(|| format!("Failed to create temp index file: {}", tmp_path.display()))?;
        // Use 8 MiB buffer for bulk sequential writes (default 8 KiB is too small)
        let mut writer = BufWriter::with_capacity(8 << 20, file);

        let target_count = u32::try_from(targets.len()).context("Target count exceeds u32::MAX")?;

        // Header: magic[8] + target_count[4] + reserved[4]
        writer
            .write_all(&FILE_MAGIC)
            .context("Failed to write magic")?;
        writer
            .write_all(&target_count.to_le_bytes())
            .context("Failed to write target count")?;
        writer
            .write_all(&0u32.to_le_bytes())
            .context("Failed to write reserved field")?;

        // Per-target metadata
        for (i, t) in targets.iter().enumerate() {
            let name_bytes = t.id.as_bytes();
            let name_len =
                u32::try_from(name_bytes.len()).context("Target name length exceeds u32::MAX")?;
            writer.write_all(&name_len.to_le_bytes())?;
            writer.write_all(&t.seq_len.to_le_bytes())?;
            writer.write_all(&offsets[i].to_le_bytes())?;
            writer.write_all(name_bytes)?;
        }

        // Pad to 8-byte alignment
        let mut cursor = FILE_HEADER_BYTES;
        for t in &targets {
            cursor += META_ENTRY_FIXED_BYTES + t.id.len();
        }
        let aligned_cursor = align_up(cursor, DATA_ALIGN);
        let pad_len = aligned_cursor - cursor;
        if pad_len > 0 {
            writer.write_all(&vec![0u8; pad_len])?;
        }

        // Global combined_seq: count[8] + data
        let seq_byte_count = combined_bases.len() as u64;
        writer.write_all(&seq_byte_count.to_le_bytes())?;
        // SAFETY: Base is #[repr(u8)], so &[Base] is layout-compatible with &[u8]
        let seq_bytes = unsafe {
            std::slice::from_raw_parts(combined_bases.as_ptr().cast::<u8>(), combined_bases.len())
        };
        writer.write_all(seq_bytes)?;

        // Pad to 8-byte alignment before SA
        let after_seq = aligned_cursor + 8 + combined_bases.len();
        let sa_start = align_up(after_seq, DATA_ALIGN);
        let sa_pad = sa_start - after_seq;
        if sa_pad > 0 {
            writer.write_all(&vec![0u8; sa_pad])?;
        }

        // Global combined_sa: count[8] + data
        let sa_entry_count = padded_sa.len() as u64;
        writer.write_all(&sa_entry_count.to_le_bytes())?;
        let sa_bytes = unsafe {
            std::slice::from_raw_parts(padded_sa.as_ptr().cast::<u8>(), padded_sa.len() * 8)
        };
        writer.write_all(sa_bytes)?;

        writer.flush().context("Failed to flush index file")?;
        drop(writer);

        std::fs::rename(&tmp_path, output).with_context(|| {
            format!(
                "Failed to finalize index file: {} -> {}",
                tmp_path.display(),
                output.display()
            )
        })?;

        Ok(())
    }

    pub fn open(path: &Path) -> Result<Self> {
        validate_readable_file(path)?;
        let file = File::open(path)
            .with_context(|| format!("Failed to open index file: {}", path.display()))?;
        let mmap = unsafe { Mmap::map(&file) }
            .with_context(|| format!("Failed to memory-map index file: {}", path.display()))?;

        let bytes = mmap.as_ref();
        if bytes.len() < FILE_HEADER_BYTES {
            bail!(
                "Index file too small ({} bytes): {}",
                bytes.len(),
                path.display()
            );
        }

        if bytes[..FILE_MAGIC.len()] != FILE_MAGIC {
            bail!("Invalid index magic in {}", path.display());
        }

        let target_count =
            u32::from_le_bytes(bytes[8..12].try_into().expect("header slice")) as usize;

        // Parse per-target metadata
        let mut cursor = FILE_HEADER_BYTES;
        let mut targets = Vec::with_capacity(target_count);

        for _ in 0..target_count {
            if cursor + META_ENTRY_FIXED_BYTES > bytes.len() {
                bail!("Truncated metadata entry in {}", path.display());
            }

            let name_len =
                u32::from_le_bytes(bytes[cursor..cursor + 4].try_into().expect("name_len"))
                    as usize;
            let seq_len =
                u32::from_le_bytes(bytes[cursor + 4..cursor + 8].try_into().expect("seq_len"))
                    as usize;
            let global_offset =
                u64::from_le_bytes(bytes[cursor + 8..cursor + 16].try_into().expect("offset"));
            cursor += META_ENTRY_FIXED_BYTES;

            let name_end = cursor
                .checked_add(name_len)
                .ok_or_else(|| anyhow!("Name overflow in {}", path.display()))?;
            if name_end > bytes.len() {
                bail!("Truncated target name in {}", path.display());
            }
            let name = std::str::from_utf8(&bytes[cursor..name_end])
                .context("Invalid UTF-8 in target name")?
                .to_owned();
            cursor = name_end;

            targets.push(TargetMeta {
                name,
                sequence_len: seq_len,
                global_offset,
            });
        }

        // Align to data section
        cursor = align_up(cursor, DATA_ALIGN);

        // Read seq_byte_count
        if cursor + 8 > bytes.len() {
            bail!("Truncated seq header in {}", path.display());
        }
        let seq_byte_count =
            u64::from_le_bytes(bytes[cursor..cursor + 8].try_into().expect("seq count")) as usize;
        let seq_data_offset = cursor + 8;
        let seq_data_end = seq_data_offset
            .checked_add(seq_byte_count)
            .ok_or_else(|| anyhow!("Seq data overflow in {}", path.display()))?;
        if seq_data_end > bytes.len() {
            bail!("Truncated seq data in {}", path.display());
        }

        // Align to SA section
        let sa_section_start = align_up(seq_data_end, DATA_ALIGN);
        if sa_section_start + 8 > bytes.len() {
            bail!("Truncated SA header in {}", path.display());
        }
        let sa_entry_count = u64::from_le_bytes(
            bytes[sa_section_start..sa_section_start + 8]
                .try_into()
                .expect("sa count"),
        ) as usize;
        let sa_data_offset = sa_section_start + 8;
        let sa_data_end = sa_data_offset
            .checked_add(sa_entry_count * 8)
            .ok_or_else(|| anyhow!("SA data overflow in {}", path.display()))?;
        if sa_data_end > bytes.len() {
            bail!("Truncated SA data in {}", path.display());
        }

        // Compute real SA length: total entries - padding
        let sa_real_len = sa_entry_count.saturating_sub(SA_CHAR_PADDING);

        // Precompute offset/seq_len arrays for global_view
        let offsets: Vec<u64> = targets.iter().map(|t| t.global_offset).collect();
        let seq_lens: Vec<u32> = targets.iter().map(|t| t.sequence_len as u32).collect();

        Ok(Self {
            path: path.to_path_buf(),
            mmap,
            targets,
            offsets,
            seq_lens,
            seq_data_offset,
            seq_byte_count,
            sa_data_offset,
            sa_entry_count,
            sa_real_len,
        })
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.targets.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.targets.is_empty()
    }

    #[inline]
    pub fn get_name(&self, idx: u32) -> &str {
        &self.targets[idx as usize].name
    }

    #[inline]
    pub fn get_sequence_len(&self, idx: usize) -> usize {
        self.targets[idx].sequence_len
    }

    pub fn index_of(&self, name: &str) -> Option<u32> {
        self.targets
            .iter()
            .position(|t| t.name == name)
            .map(|i| i as u32)
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn iter_meta(&self) -> impl Iterator<Item = (u32, &TargetMeta)> {
        self.targets.iter().enumerate().map(|(i, m)| (i as u32, m))
    }

    /// Get a global view of the combined SA and sequence for seed search.
    pub fn global_view(&self) -> GlobalView<'_> {
        let bytes = self.mmap.as_ref();

        // SAFETY: Base is #[repr(u8)], byte slice from mmap is valid
        let combined_seq = unsafe {
            std::slice::from_raw_parts(
                bytes[self.seq_data_offset..].as_ptr().cast::<Base>(),
                self.seq_byte_count,
            )
        };

        // SAFETY: u64 alignment guaranteed by DATA_ALIGN=8 padding in writer
        let combined_sa = unsafe {
            std::slice::from_raw_parts(
                bytes[self.sa_data_offset..].as_ptr().cast::<u64>(),
                self.sa_entry_count,
            )
        };

        GlobalView {
            combined_seq,
            combined_sa,
            sa_real_len: self.sa_real_len,
            offsets: &self.offsets,
            seq_lens: &self.seq_lens,
        }
    }

    /// Get per-target sequence slices for output/extension — zero-copy from mmap.
    pub fn target_seqs(&self, idx: usize) -> Result<TargetSeqs<'_>> {
        let meta = self.targets.get(idx).ok_or_else(|| {
            anyhow!(
                "Target index out of bounds: {} >= {}",
                idx,
                self.targets.len()
            )
        })?;

        let bytes = self.mmap.as_ref();
        let combined_seq = unsafe {
            std::slice::from_raw_parts(
                bytes[self.seq_data_offset..].as_ptr().cast::<Base>(),
                self.seq_byte_count,
            )
        };

        let offset = meta.global_offset as usize;
        let seq_len = meta.sequence_len;

        // Layout within global: fwd_comp[seq_len] + Gap + rc_comp[seq_len] + Gap
        let fwd_end = offset + seq_len;
        let rc_start = offset + seq_len + 1;
        let rc_end = rc_start + seq_len;

        if rc_end > combined_seq.len() {
            bail!(
                "Target '{}' sequence data out of bounds (offset={}, seq_len={}, total={})",
                meta.name,
                offset,
                seq_len,
                combined_seq.len()
            );
        }

        Ok(TargetSeqs {
            name: &meta.name,
            fwd_transformed: &combined_seq[offset..fwd_end],
            rc_transformed: &combined_seq[rc_start..rc_end],
            seq_len,
        })
    }
}

fn normalize_target_record(record: PendingTargetRecord) -> Result<Option<NormalizedTarget>> {
    let PendingTargetRecord { id, seq } = record;
    let (sequence, stats) = Sequence::normalize(&id, &seq)
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

    let sequence_rc = sequence.reverse_complement();

    // Complement sequences for canonical/wobble matching in SA
    let fwd_comp: Vec<Base> = sequence.iter().map(|b| b.complement()).collect();
    let rc_comp: Vec<Base> = sequence_rc.iter().map(|b| b.complement()).collect();

    let seq_len =
        u32::try_from(sequence.len()).context("Target sequence length too large for u32")?;

    Ok(Some(NormalizedTarget {
        id,
        seq_len,
        fwd_comp,
        rc_comp,
    }))
}

#[inline]
fn align_up(value: usize, align: usize) -> usize {
    debug_assert!(align.is_power_of_two());
    (value + (align - 1)) & !(align - 1)
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use tempfile::tempdir;

    use super::TargetStore;
    use crate::registry::TargetRegistry;

    #[test]
    fn roundtrip_build_open_global_view() {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("targets.fa");
        let index_path = dir.path().join("targets.idx");

        let mut fasta = std::fs::File::create(&fasta_path).unwrap();
        writeln!(fasta, ">chrA\nACGUGA").unwrap();
        writeln!(fasta, ">chrB\nUUUGCA").unwrap();
        drop(fasta);

        TargetStore::build_from_fasta(&fasta_path, &index_path).unwrap();
        let store = TargetStore::open(&index_path).unwrap();
        let expected = TargetRegistry::from_fasta(&fasta_path).unwrap();

        assert_eq!(store.len(), expected.len());

        // Verify global view is accessible
        let global = store.global_view();
        assert_eq!(global.offsets.len(), store.len());
        assert_eq!(global.seq_lens.len(), store.len());
        assert!(global.sa_real_len > 0);

        // Verify target_seqs for each target
        for i in 0..store.len() {
            let seqs = store.target_seqs(i).unwrap();
            let exp = expected.get(i as u32);
            assert_eq!(seqs.name, exp.name.as_str());
            assert_eq!(seqs.seq_len, exp.seq_len);
            assert_eq!(seqs.fwd_transformed.len(), exp.seq_len);
            assert_eq!(seqs.rc_transformed.len(), exp.seq_len);
        }
    }

    #[test]
    fn open_rejects_invalid_magic() {
        let dir = tempdir().unwrap();
        let bad_path = dir.path().join("bad.idx");
        // Must be at least FILE_HEADER_BYTES (16) to reach magic check
        std::fs::write(&bad_path, b"not-a-valid-idx!").unwrap();
        match TargetStore::open(&bad_path) {
            Ok(_) => panic!("Expected invalid magic error"),
            Err(err) => assert!(err.to_string().contains("Invalid index magic")),
        }
    }

    #[test]
    fn global_offsets_are_contiguous() {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("targets.fa");
        let index_path = dir.path().join("targets.idx");

        let mut fasta = std::fs::File::create(&fasta_path).unwrap();
        writeln!(fasta, ">t1\nACGU").unwrap();
        writeln!(fasta, ">t2\nGGCC").unwrap();
        writeln!(fasta, ">t3\nAA").unwrap();
        drop(fasta);

        TargetStore::build_from_fasta(&fasta_path, &index_path).unwrap();
        let store = TargetStore::open(&index_path).unwrap();
        let global = store.global_view();

        // Verify offsets are contiguous: offset[i+1] = offset[i] + 2*seq_len[i] + 2
        for i in 0..store.len() - 1 {
            let expected_next = global.offsets[i] + 2 * global.seq_lens[i] as u64 + 2;
            assert_eq!(
                global.offsets[i + 1],
                expected_next,
                "Offset mismatch at target {}",
                i
            );
        }
    }
}
