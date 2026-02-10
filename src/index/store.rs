use std::collections::HashSet;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::mem::{align_of, size_of};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use memmap2::Mmap;
use needletail::parse_fastx_file;
use rkyv::rancor::Error as RkyvError;

use crate::index::io::{validate_output_path, validate_readable_file};
use crate::index::sa::SuffixArray;
use crate::seq::Sequence;
use crate::types::{ArchivedBase, Base};

const FILE_MAGIC: [u8; 8] = *b"RSIDX2\0\0";
const FILE_HEADER_BYTES: usize = 12; // magic[8] + target_count[4]
const ENTRY_HEADER_BYTES: usize = 16; // name_len[u32] + seq_len[u32] + payload_len[u64]
const PAYLOAD_ALIGN: usize = 8;
const TARGET_COUNT_OFFSET: u64 = 8;

#[derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)]
struct TargetChunk {
    sequence: Vec<Base>,
    sequence_rc: Vec<Base>,
    forward_sa: Vec<u32>,
    reverse_sa: Vec<u32>,
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
    pub sequence: &'a [Base],
    pub sequence_rc: &'a [Base],
    pub forward_sa: &'a [u32],
    pub reverse_sa: &'a [u32],
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

            let (sequence, stats) = Sequence::normalize(&id, &rec.seq())
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

            if sequence.len() > u32::MAX as usize {
                bail!(
                    "Sequence '{}' too long for u32 suffix array ({} bases > {} max)",
                    id,
                    sequence.len(),
                    u32::MAX
                );
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
            let forward_sa = SuffixArray::try_from(&sequence)
                .with_context(|| format!("Failed to build forward SA for '{}'", id))?;
            let reverse_sa = SuffixArray::try_from(&sequence_rc)
                .with_context(|| format!("Failed to build reverse SA for '{}'", id))?;

            let chunk = TargetChunk {
                sequence: sequence.iter().copied().collect(),
                sequence_rc: sequence_rc.iter().copied().collect(),
                forward_sa: forward_sa.into_inner(),
                reverse_sa: reverse_sa.into_inner(),
            };

            let payload = rkyv::to_bytes::<RkyvError>(&chunk)
                .with_context(|| format!("Failed to archive target chunk '{}'", id))?;

            let name_len =
                u32::try_from(id.len()).context("Target id is too long for index header")?;
            let seq_len = u32::try_from(sequence.len())
                .context("Target sequence length is too large for index header")?;
            let payload_len = u64::try_from(payload.len())
                .context("Archived payload is too large for index header")?;

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

            cursor = cursor
                .checked_add(ENTRY_HEADER_BYTES)
                .and_then(|v| v.checked_add(id.len()))
                .ok_or_else(|| anyhow!("Index file cursor overflow while writing headers"))?;

            let aligned_payload_offset = align_up(cursor, PAYLOAD_ALIGN);
            let padding_len = aligned_payload_offset - cursor;
            if padding_len > 0 {
                writer
                    .write_all(&[0u8; PAYLOAD_ALIGN][..padding_len])
                    .context("Failed to write payload alignment padding")?;
            }
            writer
                .write_all(payload.as_slice())
                .context("Failed to write target payload bytes")?;
            cursor = aligned_payload_offset
                .checked_add(payload.len())
                .ok_or_else(|| anyhow!("Index file cursor overflow while writing payload"))?;

            target_count = target_count
                .checked_add(1)
                .ok_or_else(|| anyhow!("Target count exceeded u32::MAX"))?;
        }

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

        let sequence = chunk.sequence.as_slice();
        let sequence_rc = chunk.sequence_rc.as_slice();
        let forward_sa = chunk.forward_sa.as_slice();
        let reverse_sa = chunk.reverse_sa.as_slice();

        let sequence = archived_base_slice_as_native(sequence);
        let sequence_rc = archived_base_slice_as_native(sequence_rc);
        let forward_sa = archived_u32_slice_as_native(forward_sa);
        let reverse_sa = archived_u32_slice_as_native(reverse_sa);

        if sequence.len() != meta.sequence_len
            || sequence_rc.len() != meta.sequence_len
            || forward_sa.len() != meta.sequence_len
            || reverse_sa.len() != meta.sequence_len
        {
            bail!(
                "Corrupt payload for '{}': expected length {}, got seq={}, seq_rc={}, sa_fwd={}, sa_rev={}",
                meta.name,
                meta.sequence_len,
                sequence.len(),
                sequence_rc.len(),
                forward_sa.len(),
                reverse_sa.len()
            );
        }

        Ok(TargetView {
            name: &meta.name,
            sequence,
            sequence_rc,
            forward_sa,
            reverse_sa,
        })
    }

    pub fn iter_meta(&self) -> impl Iterator<Item = (u32, &TargetMeta)> {
        self.targets.iter().enumerate().map(|(i, m)| (i as u32, m))
    }
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
fn archived_base_slice_as_native(slice: &[ArchivedBase]) -> &[Base] {
    debug_assert_eq!(size_of::<ArchivedBase>(), size_of::<Base>());
    debug_assert_eq!(align_of::<ArchivedBase>(), align_of::<Base>());
    // SAFETY: `Base` is `#[repr(u8)]` and `ArchivedBase` is generated by rkyv
    // for this exact enum, with validated payload bytes (`rkyv::access`).
    // The two types are layout-compatible for read-only access.
    unsafe { std::slice::from_raw_parts(slice.as_ptr().cast::<Base>(), slice.len()) }
}

#[inline]
fn archived_u32_slice_as_native(slice: &[rkyv::primitive::ArchivedU32]) -> &[u32] {
    #[cfg(not(target_endian = "little"))]
    compile_error!("TargetStore zero-copy u32 view currently requires little-endian targets");
    debug_assert_eq!(size_of::<rkyv::primitive::ArchivedU32>(), size_of::<u32>());
    debug_assert_eq!(
        align_of::<rkyv::primitive::ArchivedU32>(),
        align_of::<u32>()
    );
    // SAFETY: `ArchivedU32` is a transparent little-endian wrapper over `u32`
    // in the archived payload. This project targets little-endian machines.
    unsafe { std::slice::from_raw_parts(slice.as_ptr().cast::<u32>(), slice.len()) }
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
            assert_eq!(got.sequence, &exp.sequence[..]);
            assert_eq!(got.sequence_rc, &exp.sequence_rc[..]);
            assert_eq!(got.forward_sa, &exp.forward_sa[..]);
            assert_eq!(got.reverse_sa, &exp.reverse_sa[..]);
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
