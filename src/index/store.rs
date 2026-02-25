use std::collections::HashSet;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::mem::{align_of, size_of};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use memmap2::Mmap;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use rkyv::rancor::Error as RkyvError;
use rkyv::util::AlignedVec;

use crate::index::io::{validate_output_path, validate_readable_file};
use crate::index::sa::SuffixArray;
use crate::seq::Sequence;
use crate::types::Base;

const FILE_MAGIC: [u8; 8] = *b"RSIDX4\0\0";
const FILE_HEADER_BYTES: usize = 12; // magic[8] + target_count[4]
const ENTRY_HEADER_BYTES: usize = 16; // name_len[u32] + seq_len[u32] + payload_len[u64]
const PAYLOAD_ALIGN: usize = 8;
const TARGET_COUNT_OFFSET: u64 = 8;

/// Number of zero-valued u64 entries appended after the real SA data.
/// These act as sentinels so `sa[suffix_pos + offset]` never goes out of
/// bounds during binary search character lookups, eliminating a branch
/// in the innermost hot loop.
pub const SA_CHAR_PADDING: usize = 256;

#[derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)]
struct TargetChunk {
    /// Forward sequence ++ [Gap] ++ reverse-complement sequence
    combined_seq: Vec<Base>,
    /// Suffix array built on combined_seq (bit-packed u64 containing character data)
    combined_sa: Vec<u64>,
    /// Length of the original (forward) sequence
    seq_len: u32,
}

#[derive(Clone, Debug)]
pub struct TargetMeta {
    pub name: String,
    pub sequence_len: usize,
    payload_offset: usize,
    payload_len: usize,
}

pub struct TargetStore {
    path: PathBuf,
    mmap: Mmap,
    targets: Vec<TargetMeta>,
}

pub struct TargetView<'a> {
    pub name: &'a str,
    /// Forward sequence ++ [Gap] ++ reverse-complement sequence (+ padding)
    pub combined_seq: &'a [Base],
    /// Suffix array built on combined_seq (packed u64 containing character data, + padding)
    pub combined_sa: &'a [u64],
    /// Number of real SA entries (excluding padding sentinels)
    pub sa_real_len: usize,
    /// Length of the original (forward) sequence
    pub seq_len: usize,
}

struct PendingTargetRecord {
    id: String,
    seq: Vec<u8>,
}

struct BuiltTargetRecord {
    id: String,
    seq_len: u32,
    payload: AlignedVec,
}

impl TargetStore {
    pub fn build_from_fasta(input: &Path, output: &Path) -> Result<()> {
        validate_readable_file(input)?;
        validate_output_path(output)?;

        let output_name = output
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "target.idx".to_string());
        let tmp_path = output.with_file_name(format!("{output_name}.tmp"));

        let file = File::create(&tmp_path)
            .with_context(|| format!("Failed to create temp index file: {}", tmp_path.display()))?;
        let mut writer = BufWriter::new(file);
        writer
            .write_all(&FILE_MAGIC)
            .context("Failed to write index magic")?;
        writer
            .write_all(&0u32.to_le_bytes())
            .context("Failed to write target count placeholder")?;

        let mut target_count = 0u32;
        let mut cursor = FILE_HEADER_BYTES;
        let mut seen = HashSet::new();
        let batch_size = rayon::current_num_threads().max(1) * 4;
        let mut pending = Vec::with_capacity(batch_size);

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

            if pending.len() >= batch_size {
                flush_pending_batch(&mut pending, &mut writer, &mut cursor, &mut target_count)?;
            }
        }

        flush_pending_batch(&mut pending, &mut writer, &mut cursor, &mut target_count)?;

        if target_count == 0 {
            bail!(
                "All sequences were empty after normalization in {}",
                input.display()
            );
        }

        writer.flush().context("Failed to flush index file")?;
        drop(writer);

        {
            let mut patch_file = OpenOptions::new()
                .write(true)
                .open(&tmp_path)
                .with_context(|| {
                    format!("Failed to reopen temp index file: {}", tmp_path.display())
                })?;
            patch_file
                .seek(SeekFrom::Start(TARGET_COUNT_OFFSET))
                .context("Failed to seek to target count field")?;
            patch_file
                .write_all(&target_count.to_le_bytes())
                .context("Failed to patch target count field")?;
            patch_file
                .flush()
                .context("Failed to flush patched header")?;
        }

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
        // SAFETY: The file descriptor stays alive for the mmap creation call and
        // the resulting read-only mapping is stored immutably in TargetStore.
        // We do not mutate the mapped file through this process.
        let mmap = unsafe { Mmap::map(&file) }
            .with_context(|| format!("Failed to memory-map index file: {}", path.display()))?;

        let targets = parse_index_metadata(mmap.as_ref(), path)?;
        Ok(Self {
            path: path.to_path_buf(),
            mmap,
            targets,
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

    pub fn target_view(&self, idx: usize) -> Result<TargetView<'_>> {
        let meta = self.targets.get(idx).ok_or_else(|| {
            anyhow!(
                "Target index out of bounds: {} >= {}",
                idx,
                self.targets.len()
            )
        })?;

        let payload_end = meta
            .payload_offset
            .checked_add(meta.payload_len)
            .ok_or_else(|| anyhow!("Payload span overflow for '{}'", meta.name))?;
        let payload = self
            .mmap
            .get(meta.payload_offset..payload_end)
            .ok_or_else(|| anyhow!("Payload out of bounds for '{}'", meta.name))?;

        let chunk = access_archived_chunk(payload)
            .with_context(|| format!("Failed to access target payload '{}'", meta.name))?;

        let combined_seq = archived_base_slice_as_native(chunk.combined_seq.as_slice());
        let combined_sa = archived_u64_slice_as_native(chunk.combined_sa.as_slice());
        let seq_len = chunk.seq_len.to_native() as usize;
        let min_seq_len = 2 * seq_len + 2; // fwd + Gap + rc + sentinel
        let min_sa_len = 2 * seq_len + 1; // SA built before sentinel

        if seq_len != meta.sequence_len
            || combined_seq.len() < min_seq_len
            || combined_sa.len() < min_sa_len
        {
            bail!(
                "Corrupt payload for '{}': expected seq_len={}, combined_seq>={}, combined_sa>={}, got seq_len={}, combined_seq={}, combined_sa={}",
                meta.name,
                meta.sequence_len,
                min_seq_len,
                min_sa_len,
                seq_len,
                combined_seq.len(),
                combined_sa.len()
            );
        }

        Ok(TargetView {
            name: &meta.name,
            combined_seq,
            combined_sa,
            sa_real_len: min_sa_len,
            seq_len,
        })
    }

    pub fn iter_meta(&self) -> impl Iterator<Item = (u32, &TargetMeta)> {
        self.targets.iter().enumerate().map(|(i, m)| (i as u32, m))
    }
}

fn flush_pending_batch(
    pending: &mut Vec<PendingTargetRecord>,
    writer: &mut BufWriter<File>,
    cursor: &mut usize,
    target_count: &mut u32,
) -> Result<()> {
    if pending.is_empty() {
        return Ok(());
    }

    let batch = std::mem::replace(pending, Vec::with_capacity(pending.capacity().max(1)));
    let batch_results: Vec<Result<Option<BuiltTargetRecord>>> =
        batch.into_par_iter().map(build_target_record).collect();

    for record in batch_results {
        if let Some(record) = record? {
            write_target_entry(writer, cursor, record)?;
            *target_count = target_count
                .checked_add(1)
                .ok_or_else(|| anyhow!("Target count exceeded u32::MAX"))?;
        }
    }

    Ok(())
}

fn build_target_record(record: PendingTargetRecord) -> Result<Option<BuiltTargetRecord>> {
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

    // Complement the sequences before building the combined SA.
    // This allows direct canonical/wobble matching in the transformed alphabet.
    let sequence_comp: Vec<Base> = sequence.iter().map(|b| b.complement()).collect();
    let sequence_rc_comp: Vec<Base> = sequence_rc.iter().map(|b| b.complement()).collect();

    // Build combined sequence: fwd_comp ++ [Gap] ++ rc_comp
    let seq_len = sequence.len();
    let seq_len_u32 =
        u32::try_from(seq_len).context("Target sequence length is too large for index header")?;
    let mut combined_bases: Vec<Base> = Vec::with_capacity(2 * seq_len + 2);
    combined_bases.extend_from_slice(&sequence_comp);
    combined_bases.push(Base::Gap);
    combined_bases.extend_from_slice(&sequence_rc_comp);
    let combined_seq = Sequence::from(combined_bases);

    let combined_sa = SuffixArray::try_from(&combined_seq)
        .with_context(|| format!("Failed to build combined SA for '{}'", id))?;

    // Append sentinel Gap byte after SA construction, plus SA_CHAR_PADDING
    // extra Gap bytes so that `sa[suffix_pos + offset]` never goes out of
    // bounds during binary search character lookups.
    let mut combined_seq_vec: Vec<Base> = combined_seq.iter().copied().collect();
    combined_seq_vec.resize(combined_seq_vec.len() + 1 + SA_CHAR_PADDING, Base::Gap);

    // Pad the SA with zero entries (pos=0, char=Gap) so unchecked access is
    // safe for any offset up to SA_CHAR_PADDING.
    let mut padded_sa = combined_sa.into_inner();
    padded_sa.resize(padded_sa.len() + SA_CHAR_PADDING, 0u64);

    let chunk = TargetChunk {
        combined_seq: combined_seq_vec,
        combined_sa: padded_sa,
        seq_len: seq_len_u32,
    };

    let payload = rkyv::to_bytes::<RkyvError>(&chunk)
        .with_context(|| format!("Failed to archive target chunk '{}'", id))?;

    Ok(Some(BuiltTargetRecord {
        id,
        seq_len: seq_len_u32,
        payload,
    }))
}

fn write_target_entry(
    writer: &mut BufWriter<File>,
    cursor: &mut usize,
    record: BuiltTargetRecord,
) -> Result<()> {
    let BuiltTargetRecord {
        id,
        seq_len,
        payload,
    } = record;
    let name_len = u32::try_from(id.len()).context("Target id is too long for index header")?;
    let payload_len =
        u64::try_from(payload.len()).context("Archived payload is too large for index header")?;

    writer
        .write_all(&name_len.to_le_bytes())
        .context("Failed to write target name length")?;
    writer
        .write_all(&seq_len.to_le_bytes())
        .context("Failed to write target sequence length")?;
    writer
        .write_all(&payload_len.to_le_bytes())
        .context("Failed to write target payload length")?;
    writer
        .write_all(id.as_bytes())
        .context("Failed to write target name bytes")?;

    *cursor = cursor
        .checked_add(ENTRY_HEADER_BYTES)
        .and_then(|v| v.checked_add(id.len()))
        .ok_or_else(|| anyhow!("Index file cursor overflow while writing headers"))?;

    let aligned_payload_offset = align_up(*cursor, PAYLOAD_ALIGN);
    let padding_len = aligned_payload_offset - *cursor;
    if padding_len > 0 {
        writer
            .write_all(&[0u8; PAYLOAD_ALIGN][..padding_len])
            .context("Failed to write payload alignment padding")?;
    }
    writer
        .write_all(payload.as_slice())
        .context("Failed to write target payload bytes")?;
    *cursor = aligned_payload_offset
        .checked_add(payload.len())
        .ok_or_else(|| anyhow!("Index file cursor overflow while writing payload"))?;

    Ok(())
}

fn parse_index_metadata(bytes: &[u8], path: &Path) -> Result<Vec<TargetMeta>> {
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

    let target_count = u32::from_le_bytes(
        bytes[FILE_MAGIC.len()..FILE_HEADER_BYTES]
            .try_into()
            .expect("header slice has fixed width"),
    ) as usize;
    let max_possible_entries = (bytes.len() - FILE_HEADER_BYTES) / ENTRY_HEADER_BYTES;
    if target_count > max_possible_entries {
        bail!(
            "Corrupt index header in {}: target count {} exceeds possible entries {}",
            path.display(),
            target_count,
            max_possible_entries
        );
    }

    let mut cursor = FILE_HEADER_BYTES;
    let mut targets = Vec::with_capacity(target_count);

    for _ in 0..target_count {
        if cursor + ENTRY_HEADER_BYTES > bytes.len() {
            bail!("Truncated entry header in {}", path.display());
        }

        let name_len = u32::from_le_bytes(
            bytes[cursor..cursor + 4]
                .try_into()
                .expect("name_len slice has fixed width"),
        ) as usize;
        let seq_len = u32::from_le_bytes(
            bytes[cursor + 4..cursor + 8]
                .try_into()
                .expect("seq_len slice has fixed width"),
        ) as usize;
        let payload_len = u64::from_le_bytes(
            bytes[cursor + 8..cursor + 16]
                .try_into()
                .expect("payload_len slice has fixed width"),
        ) as usize;
        cursor += ENTRY_HEADER_BYTES;

        let name_end = cursor
            .checked_add(name_len)
            .ok_or_else(|| anyhow!("Entry name overflow while parsing {}", path.display()))?;
        if name_end > bytes.len() {
            bail!("Truncated target name in {}", path.display());
        }
        let name = std::str::from_utf8(&bytes[cursor..name_end])
            .context("Invalid UTF-8 in target name")?
            .to_owned();
        cursor = name_end;

        let payload_offset = align_up(cursor, PAYLOAD_ALIGN);
        if payload_offset > bytes.len() {
            bail!("Truncated payload alignment padding in {}", path.display());
        }
        let payload_end = cursor
            .checked_add(payload_offset - cursor)
            .and_then(|v| v.checked_add(payload_len))
            .ok_or_else(|| anyhow!("Entry payload overflow while parsing {}", path.display()))?;
        if payload_end > bytes.len() {
            bail!("Truncated target payload in {}", path.display());
        }

        targets.push(TargetMeta {
            name,
            sequence_len: seq_len,
            payload_offset,
            payload_len,
        });
        cursor = payload_end;
    }

    if cursor != bytes.len() {
        bail!(
            "Unexpected trailing bytes in index file {} (parsed {} of {} bytes)",
            path.display(),
            cursor,
            bytes.len()
        );
    }

    Ok(targets)
}

#[inline]
fn align_up(value: usize, align: usize) -> usize {
    debug_assert!(align.is_power_of_two());
    (value + (align - 1)) & !(align - 1)
}

#[inline]
fn access_archived_chunk(payload: &[u8]) -> Result<&ArchivedTargetChunk> {
    #[cfg(debug_assertions)]
    {
        rkyv::access::<ArchivedTargetChunk, RkyvError>(payload).context("rkyv validation failed")
    }

    #[cfg(not(debug_assertions))]
    {
        // SAFETY: Payload comes from a bounded chunk span parsed from this
        // index format. Index files are produced by `build_from_fasta` in this
        // crate, and release builds prioritize load throughput over byte-level
        // validation on every open.
        Ok(unsafe { rkyv::access_unchecked::<ArchivedTargetChunk>(payload) })
    }
}

#[inline]
fn archived_base_slice_as_native(slice: &[crate::types::ArchivedBase]) -> &[Base] {
    debug_assert_eq!(size_of::<crate::types::ArchivedBase>(), size_of::<Base>());
    debug_assert_eq!(align_of::<crate::types::ArchivedBase>(), align_of::<Base>());
    // SAFETY: `Base` is `#[repr(u8)]` and `ArchivedBase` is generated by rkyv
    // for this exact enum, with validated payload bytes (`rkyv::access`).
    // The two types are layout-compatible for read-only access.
    unsafe { std::slice::from_raw_parts(slice.as_ptr().cast::<Base>(), slice.len()) }
}

#[inline]
fn archived_u64_slice_as_native(slice: &[rkyv::primitive::ArchivedU64]) -> &[u64] {
    #[cfg(not(target_endian = "little"))]
    compile_error!("TargetStore zero-copy u64 view currently requires little-endian targets");
    debug_assert_eq!(size_of::<rkyv::primitive::ArchivedU64>(), size_of::<u64>());
    debug_assert_eq!(
        align_of::<rkyv::primitive::ArchivedU64>(),
        align_of::<u64>()
    );
    // SAFETY: `ArchivedU64` is a transparent little-endian wrapper over `u64`
    // in the archived payload. This project targets little-endian machines.
    unsafe { std::slice::from_raw_parts(slice.as_ptr().cast::<u64>(), slice.len()) }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use tempfile::tempdir;

    use super::TargetStore;
    use crate::registry::TargetRegistry;

    #[test]
    fn roundtrip_build_open_load_matches_registry() {
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
        for i in 0..expected.len() {
            let got = store.target_view(i).unwrap();
            let exp = expected.get(i as u32);
            assert_eq!(got.name, exp.name.as_str());
            assert_eq!(got.seq_len, exp.seq_len);
            // Suffix Arrays will differ because TargetStore builds on COMPLEMENT
            // whereas Registry builds on original for backwards-compat in tests.
            // We just verify the store loaded correctly and real length matches.
            assert_eq!(got.sa_real_len, exp.combined_sa.len());
            // TargetStore adds SA_CHAR_PADDING sentinels
            assert_eq!(
                got.combined_sa.len(),
                exp.combined_sa.len() + crate::index::store::SA_CHAR_PADDING
            );
        }
    }

    #[test]
    fn open_rejects_invalid_magic() {
        let dir = tempdir().unwrap();
        let bad_path = dir.path().join("bad.idx");
        std::fs::write(&bad_path, b"not-an-index").unwrap();
        match TargetStore::open(&bad_path) {
            Ok(_) => panic!("Expected invalid magic error"),
            Err(err) => assert!(err.to_string().contains("Invalid index magic")),
        }
    }
}
