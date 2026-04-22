use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use memmap2::Mmap;
use needletail::parse_fastx_file;

use crate::index::io::validate_output_path;
use crate::index::sa::SuffixArray;
use crate::seed::SeedView;
use crate::seq::Sequence;
use crate::types::{Base, Strand};

// =============================================================================
// FORMAT CONSTANTS
// =============================================================================

const FILE_MAGIC: [u8; 8] = *b"RSIDX6\0\0";

/// Fixed file header: magic[8] + target_count[4] + reserved[4]
const FILE_HEADER_BYTES: usize = 16;

/// Per-target metadata entry: name_len[4] + seq_len[4] + offset[8]
const META_ENTRY_FIXED_BYTES: usize = 16;

const DATA_ALIGN: usize = 8;
const ZERO_PAD: [u8; DATA_ALIGN] = [0u8; DATA_ALIGN];

/// Number of zero-valued u64 entries appended after the real SA data.
/// These act as sentinels so `sa[suffix_pos + offset]` never goes out of
/// bounds during binary search character lookups, eliminating a branch
/// in the innermost hot loop.
pub const SA_CHAR_PADDING: usize = 256;

// =============================================================================
// PUBLIC TYPES
// =============================================================================

pub struct TargetStore {
    mmap: Mmap,
    names: Vec<String>,
    /// Start offset of each target block in global combined_seq.
    offsets: Vec<usize>,
    /// Original forward sequence length of each target.
    seq_lens: Vec<usize>,
    /// Byte offset into mmap where the global combined_seq begins.
    seq_data_offset: usize,
    /// Total number of Base entries in global combined_seq (including padding).
    seq_byte_count: usize,
    /// Byte offset into mmap where the global combined_sa begins.
    sa_data_offset: usize,
    /// Total number of u64 entries in global combined_sa (including padding).
    sa_entry_count: usize,
    /// Number of real SA entries (excluding SA_CHAR_PADDING sentinels).
    sa_real_len: usize,
}


impl TargetStore {
    /// Build a new RSIDX6 index from a FASTA file.
    ///
    /// Pipeline:
    /// 1. Parse all FASTA records, normalize, complement
    /// 2. Concatenate into global combined_seq with Gap separators
    /// 3. Build single global SA
    /// 4. Write flat binary: header -> metadata -> seq -> SA
    pub fn build(input: &Path, output: &Path) -> Result<()> {
        validate_output_path(output)?;

        let mut seen = HashSet::new();
        let mut targets: Vec<(String, u32, u64)> = Vec::new(); // (id, seq_len, global_offset)
        let mut combined_bases: Vec<Base> = Vec::new();
        let mut metadata_bytes = FILE_HEADER_BYTES;

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

            let (sequence, stats) = Sequence::normalize(&id, rec.seq().as_ref())
                .with_context(|| format!("Failed to normalize sequence '{}'", id))?;

            if sequence.is_empty() {
                log::warn!(
                    "Skipping empty sequence after normalization: '{}' (removed_gaps={}, converted_to_n={})",
                    id,
                    stats.removed_gaps,
                    stats.converted_to_n
                );
                continue;
            }

            if stats.removed_gaps > 0 || stats.converted_to_n > 0 {
                log::debug!(
                    "Normalized sequence '{}': removed_gaps={}, converted_to_n={}",
                    id,
                    stats.removed_gaps,
                    stats.converted_to_n
                );
            }

            let seq_len = u32::try_from(sequence.len())
                .context("Target sequence length too large for u32")?;
            let offset = combined_bases.len() as u64;
            metadata_bytes = metadata_bytes
                .checked_add(META_ENTRY_FIXED_BYTES + id.len())
                .ok_or_else(|| anyhow!("Metadata size overflow while building index"))?;
            targets.push((id, seq_len, offset));

            // Store transformed sequence layout:
            // fwd_comp[seq_len] + Gap + rc_comp[seq_len] + Gap
            // where *_comp is complemented for canonical/wobble seed matching.
            combined_bases.reserve(2 * seq_len as usize + 2);
            combined_bases.extend(sequence.iter().copied().map(Base::complement));
            combined_bases.push(Base::Gap);
            combined_bases.extend(sequence[..].iter().rev().copied());
            combined_bases.push(Base::Gap);
        }

        if targets.is_empty() {
            bail!(
                "All sequences were empty after normalization in {}",
                input.display()
            );
        }

        // Build single global SA directly from &[Base]
        // (avoids cloning combined_bases and extra conversion allocations)
        let combined_sa = SuffixArray::try_from(combined_bases.as_slice())
            .context("Failed to build global suffix array")?;

        // Append sentinels to seq/SA for branchless SA character lookup.
        combined_bases.resize(combined_bases.len() + SA_CHAR_PADDING, Base::Gap);
        let mut padded_sa = combined_sa.into_inner();
        padded_sa.resize(padded_sa.len() + SA_CHAR_PADDING, 0u64);

        write_index_file(
            output,
            &targets,
            metadata_bytes,
            &combined_bases,
            &padded_sa,
        )
    }

    pub fn open(path: &Path) -> Result<Self> {
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
        let mut names = Vec::with_capacity(target_count);
        let mut seq_lens = Vec::with_capacity(target_count);
        let mut offsets = Vec::with_capacity(target_count);

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
                u64::from_le_bytes(bytes[cursor + 8..cursor + 16].try_into().expect("offset"))
                    as usize;
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

            names.push(name);
            seq_lens.push(seq_len);
            offsets.push(global_offset);
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
        let sa_byte_count = sa_entry_count
            .checked_mul(8)
            .ok_or_else(|| anyhow!("SA data overflow in {}", path.display()))?;
        let sa_data_end = sa_data_offset
            .checked_add(sa_byte_count)
            .ok_or_else(|| anyhow!("SA data overflow in {}", path.display()))?;
        if sa_data_end > bytes.len() {
            bail!("Truncated SA data in {}", path.display());
        }

        // Compute real SA length: total entries - padding
        let sa_real_len = sa_entry_count.saturating_sub(SA_CHAR_PADDING);

        Ok(Self {
            mmap,
            names,
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
        self.names.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    #[inline]
    pub fn get_name(&self, idx: usize) -> &str {
        &self.names[idx]
    }

    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }

    /// Get a seed-search view of the combined SA and sequence.
    pub fn view(&self) -> SeedView<'_> {
        SeedView {
            combined_seq: self.combined_seq(),
            combined_sa: self.combined_sa(),
            sa_real_len: self.sa_real_len,
            offsets: &self.offsets,
            seq_lens: &self.seq_lens,
        }
    }

    /// Return the forward and reverse-complement slices for a target entry.
    ///
    /// Block layout: `fwd_comp[seq_len] + Gap + rc_comp[seq_len] + Gap`.
    #[inline]
    pub fn target_slices(&self, target_idx: usize) -> (&[Base], &[Base], usize) {
        let seq_len = self.seq_lens[target_idx];
        let offset = self.offsets[target_idx];
        let combined_seq = self.combined_seq();
        let t_fwd = &combined_seq[offset..offset + seq_len];
        let t_rc = &combined_seq[offset + seq_len + 1..offset + 2 * seq_len + 1];
        (t_fwd, t_rc, seq_len)
    }

    /// Map a local position within a target block to strand and view-normalized start.
    ///
    /// Block layout: `fwd_comp[seq_len] + Gap + rc[seq_len] + Gap`.
    /// Returns `None` if the position falls on a Gap or the seed overflows the block.
    #[inline]
    pub fn map_target_pos(
        local_pos: usize,
        seq_len: usize,
        seed_len: usize,
    ) -> Option<(Strand, usize)> {
        if local_pos < seq_len {
            if local_pos + seed_len > seq_len {
                return None;
            }
            return Some((Strand::Reverse, seq_len - (local_pos + seed_len)));
        }

        let rc_start = seq_len + 1;
        let rc_end = rc_start + seq_len;
        if local_pos >= rc_start && local_pos < rc_end {
            let rc_pos = local_pos - rc_start;
            if rc_pos + seed_len > seq_len {
                return None;
            }
            return Some((Strand::Forward, seq_len - (rc_pos + seed_len)));
        }

        None
    }

    /// Get per-target sequence slices for output/extension.
    ///
    /// Returns `(name, fwd_transformed, rc_transformed, seq_len)` where each
    /// transformed slice is from the mmap-backed combined sequence.
    pub fn target_seqs(&self, idx: usize) -> Result<(&str, &[Base], &[Base], usize)> {
        let target_count = self.names.len();
        if idx >= target_count {
            bail!("Target index out of bounds: {} >= {}", idx, target_count);
        }

        let name = self.names[idx].as_str();
        let seq_len = self.seq_lens[idx];
        let offset = self.offsets[idx];
        let combined_seq = self.combined_seq();

        // Layout within global: fwd_comp[seq_len] + Gap + rc_comp[seq_len] + Gap
        let fwd_end = offset + seq_len;
        let rc_start = offset + seq_len + 1;
        let rc_end = rc_start + seq_len;

        if rc_end > combined_seq.len() {
            bail!(
                "Target '{}' sequence data out of bounds (offset={}, seq_len={}, total={})",
                name,
                offset,
                seq_len,
                combined_seq.len()
            );
        }

        Ok((
            name,
            &combined_seq[offset..fwd_end],
            &combined_seq[rc_start..rc_end],
            seq_len,
        ))
    }

    #[inline]
    fn combined_seq(&self) -> &[Base] {
        let bytes = self.mmap.as_ref();
        // SAFETY: Base is #[repr(u8)] and seq_data_offset/seq_byte_count are validated on open.
        unsafe {
            std::slice::from_raw_parts(
                bytes[self.seq_data_offset..].as_ptr().cast::<Base>(),
                self.seq_byte_count,
            )
        }
    }

    #[inline]
    fn combined_sa(&self) -> &[u64] {
        let bytes = self.mmap.as_ref();
        // SAFETY: u64 alignment guaranteed by DATA_ALIGN=8 padding in writer.
        unsafe {
            std::slice::from_raw_parts(
                bytes[self.sa_data_offset..].as_ptr().cast::<u64>(),
                self.sa_entry_count,
            )
        }
    }
}

fn write_index_file(
    output: &Path,
    targets: &[(String, u32, u64)],
    metadata_bytes: usize,
    seq: &[Base],
    sa: &[u64],
) -> Result<()> {
    let output_name = output
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "target.idx".to_string());
    let tmp_path = output.with_file_name(format!("{output_name}.tmp"));

    let file = File::create(&tmp_path)
        .with_context(|| format!("Failed to create temp index file: {}", tmp_path.display()))?;
    // Use 8 MiB buffer for bulk sequential writes (default 8 KiB is too small).
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
    for (id, seq_len, offset) in targets.iter() {
        let name_bytes = id.as_bytes();
        let name_len =
            u32::try_from(name_bytes.len()).context("Target name length exceeds u32::MAX")?;
        writer.write_all(&name_len.to_le_bytes())?;
        writer.write_all(&seq_len.to_le_bytes())?;
        writer.write_all(&offset.to_le_bytes())?;
        writer.write_all(name_bytes)?;
    }

    // Pad to 8-byte alignment
    let aligned_cursor = align_up(metadata_bytes, DATA_ALIGN);
    let pad_len = aligned_cursor - metadata_bytes;
    if pad_len > 0 {
        writer.write_all(&ZERO_PAD[..pad_len])?;
    }

    // Global combined_seq: count[8] + data
    let seq_byte_count = seq.len() as u64;
    writer.write_all(&seq_byte_count.to_le_bytes())?;
    // SAFETY: Base is #[repr(u8)], so &[Base] is layout-compatible with &[u8].
    let seq_bytes = unsafe { std::slice::from_raw_parts(seq.as_ptr().cast::<u8>(), seq.len()) };
    writer.write_all(seq_bytes)?;

    // Pad to 8-byte alignment before SA
    let after_seq = aligned_cursor + 8 + seq.len();
    let sa_start = align_up(after_seq, DATA_ALIGN);
    let sa_pad = sa_start - after_seq;
    if sa_pad > 0 {
        writer.write_all(&ZERO_PAD[..sa_pad])?;
    }

    // Global combined_sa: count[8] + data
    let sa_entry_count = sa.len() as u64;
    writer.write_all(&sa_entry_count.to_le_bytes())?;
    let sa_bytes = unsafe { std::slice::from_raw_parts(sa.as_ptr().cast::<u8>(), sa.len() * 8) };
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

    #[test]
    fn roundtrip_build_open_target_view() {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("targets.fa");
        let index_path = dir.path().join("targets.idx");

        let mut fasta = std::fs::File::create(&fasta_path).unwrap();
        writeln!(fasta, ">chrA\nACGUGA").unwrap();
        writeln!(fasta, ">chrB\nUUUGCA").unwrap();
        drop(fasta);

        TargetStore::build(&fasta_path, &index_path).unwrap();
        let store = TargetStore::open(&index_path).unwrap();
        let expected = [("chrA", 6usize), ("chrB", 6usize)];
        assert_eq!(store.len(), expected.len());

        // Verify target view is accessible
        let target = store.view();
        assert_eq!(target.offsets.len(), store.len());
        assert_eq!(target.seq_lens.len(), store.len());
        assert!(target.sa_real_len > 0);

        // Verify target_seqs for each target
        for (i, (expected_name, expected_seq_len)) in expected.iter().enumerate().take(store.len())
        {
            let (name, fwd, rc, seq_len) = store.target_seqs(i).unwrap();
            assert_eq!(name, *expected_name);
            assert_eq!(seq_len, *expected_seq_len);
            assert_eq!(fwd.len(), *expected_seq_len);
            assert_eq!(rc.len(), *expected_seq_len);
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
    fn map_target_pos_normalizes_into_strand_view() {
        use super::TargetStore;
        use crate::types::Strand;

        let seq_len = 5;
        let seed_len = 2;

        assert_eq!(
            TargetStore::map_target_pos(1, seq_len, seed_len),
            Some((Strand::Reverse, 2))
        );
        assert_eq!(
            TargetStore::map_target_pos(seq_len + 1 + 2, seq_len, seed_len),
            Some((Strand::Forward, 1))
        );
        assert_eq!(TargetStore::map_target_pos(seq_len, seq_len, seed_len), None);
        assert_eq!(
            TargetStore::map_target_pos(2 * seq_len + 1, seq_len, seed_len),
            None
        );
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

        TargetStore::build(&fasta_path, &index_path).unwrap();
        let store = TargetStore::open(&index_path).unwrap();
        let target = store.view();

        // Verify offsets are contiguous: offset[i+1] = offset[i] + 2*seq_len[i] + 2
        for i in 0..store.len() - 1 {
            let expected_next = target.offsets[i] + 2 * target.seq_lens[i] + 2;
            assert_eq!(
                target.offsets[i + 1],
                expected_next,
                "Offset mismatch at target {}",
                i
            );
        }
    }
}
