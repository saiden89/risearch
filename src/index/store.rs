use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use fs_err::File;
use memmap2::Mmap;

use crate::fastx::{normalize_record, read_and_validate_fasta};
use crate::index::io::validate_output_path;
use crate::index::sa::SuffixArray;
use crate::index::view::RegistryView;
use crate::types::{Base, Strand};

const TARGET_REGISTRY_VERSION: u32 = 1;

/// Number of zero-valued entries appended after the real sequence and SA data.
/// These sentinels keep `sa[suffix_pos + offset]` branchless in the seed hot path.
pub const SA_CHAR_PADDING: usize = 256;

/// Runtime handle for a target index.
///
/// The mapped file contains a private rkyv `TargetStore`; this type owns the
/// mmap lifetime and carries the small native sidecars needed by `RegistryView`.
pub struct TargetRegistry {
    mmap: Mmap,
    offsets: Vec<usize>,
    seq_lens: Vec<usize>,
    sa_real_len: usize,
}

#[derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)]
struct TargetStore {
    version: u32,
    targets: Vec<TargetRecord>,
    combined_seq: Vec<u8>,
    combined_sa: Vec<u64>,
}

#[derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)]
struct TargetRecord {
    name: String,
    offset: u64,
    seq_len: u64,
}

impl TargetRegistry {
    /// Build a target index from a FASTA file and write it as an rkyv archive.
    pub fn build(input: &Path, output: &Path) -> Result<()> {
        validate_output_path(output)?;

        let records = read_and_validate_fasta(input)?;
        let mut targets = Vec::new();
        let mut combined_bases = Vec::new();

        for (id, raw_seq) in records {
            let Some(sequence) = normalize_record(&id, &raw_seq)? else {
                continue;
            };

            let seq_len = sequence.len() as u64;
            let offset = combined_bases.len() as u64;

            targets.push(TargetRecord {
                name: id,
                offset,
                seq_len,
            });

            combined_bases.reserve(2 * sequence.len() + 2);
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

        let combined_sa = SuffixArray::try_from(combined_bases.as_slice())
            .context("Failed to build global suffix array")?;

        combined_bases.resize(combined_bases.len() + SA_CHAR_PADDING, Base::Gap);
        let mut combined_sa = combined_sa.into_inner();
        combined_sa.resize(combined_sa.len() + SA_CHAR_PADDING, 0u64);

        let store = TargetStore {
            version: TARGET_REGISTRY_VERSION,
            targets,
            combined_seq: combined_bases.into_iter().map(Base::as_u8).collect(),
            combined_sa,
        };

        write_target_registry(output, &store)
    }

    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open index file: {}", path.display()))?;
        let mmap = unsafe { Mmap::map(&file) }
            .with_context(|| format!("Failed to memory-map index file: {}", path.display()))?;

        let root = rkyv::access::<ArchivedTargetStore, rkyv::rancor::Error>(mmap.as_ref())
            .with_context(|| format!("Invalid target index archive: {}", path.display()))?;
        let (offsets, seq_lens, sa_real_len) = validate_store(root, path)?;

        Ok(Self {
            mmap,
            offsets,
            seq_lens,
            sa_real_len,
        })
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.root().targets.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.root().targets.is_empty()
    }

    #[inline]
    pub fn get_name(&self, idx: usize) -> &str {
        self.root().targets[idx].name.as_str()
    }

    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.root()
            .targets
            .iter()
            .position(|target| target.name.as_str() == name)
    }

    /// Get a seed-search view of the combined SA and sequence.
    pub fn view(&self) -> RegistryView<'_> {
        RegistryView {
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
    ///
    /// # Panics
    /// Panics if `idx >= self.len()`.
    pub fn target_seqs(&self, idx: usize) -> (&str, &[Base], &[Base], usize) {
        let name = self.get_name(idx);
        let seq_len = self.seq_lens[idx];
        let offset = self.offsets[idx];
        let combined_seq = self.combined_seq();
        let fwd_end = offset + seq_len;
        let rc_start = fwd_end + 1;
        let rc_end = rc_start + seq_len;
        (
            name,
            &combined_seq[offset..fwd_end],
            &combined_seq[rc_start..rc_end],
            seq_len,
        )
    }

    #[inline]
    fn root(&self) -> &ArchivedTargetStore {
        // SAFETY: `open` validates the archive before constructing `TargetRegistry`.
        unsafe { rkyv::access_unchecked::<ArchivedTargetStore>(self.mmap.as_ref()) }
    }

    #[inline]
    fn combined_seq(&self) -> &[Base] {
        let bytes = self.root().combined_seq.as_slice();
        // SAFETY: build path only writes Base::as_u8() values; discriminants are always valid.
        unsafe { std::slice::from_raw_parts(bytes.as_ptr().cast::<Base>(), bytes.len()) }
    }

    #[inline]
    fn combined_sa(&self) -> &[u64] {
        archived_u64_as_native(self.root().combined_sa.as_slice())
    }
}

fn write_target_registry(output: &Path, store: &TargetStore) -> Result<()> {
    let output_name = output
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "target.idx".to_string());
    let tmp_path = output.with_file_name(format!("{output_name}.tmp"));

    let bytes =
        rkyv::to_bytes::<rkyv::rancor::Error>(store).context("Failed to serialize target index")?;
    fs_err::write(&tmp_path, bytes.as_slice())
        .with_context(|| format!("Failed to write target index: {}", tmp_path.display()))?;
    fs_err::rename(&tmp_path, output).with_context(|| {
        format!(
            "Failed to finalize index file: {} -> {}",
            tmp_path.display(),
            output.display()
        )
    })?;
    Ok(())
}

fn validate_store(
    root: &ArchivedTargetStore,
    path: &Path,
) -> Result<(Vec<usize>, Vec<usize>, usize)> {
    if root.version.to_native() != TARGET_REGISTRY_VERSION {
        bail!(
            "Unsupported target index version {} in {}",
            root.version.to_native(),
            path.display()
        );
    }

    let seq = root.combined_seq.as_slice();
    let sa = root.combined_sa.as_slice();
    if seq.len() < SA_CHAR_PADDING || sa.len() < SA_CHAR_PADDING {
        bail!("Target index missing sentinel padding: {}", path.display());
    }
    if seq.len() != sa.len() {
        bail!(
            "Target index sequence/SA length mismatch (seq={}, sa={}): {}",
            seq.len(),
            sa.len(),
            path.display()
        );
    }
    let real_len = sa.len() - SA_CHAR_PADDING;
    let mut offsets = Vec::with_capacity(root.targets.len());
    let mut seq_lens = Vec::with_capacity(root.targets.len());
    let mut expected_offset = 0usize;

    for target in root.targets.iter() {
        let name = target.name.as_str();
        let offset = usize::try_from(target.offset.to_native())
            .context("Target offset does not fit in usize")?;
        let seq_len = target.seq_len.to_native() as usize;

        if offset != expected_offset {
            bail!(
                "Non-contiguous target offset for '{}' (got {}, expected {}) in {}",
                name,
                offset,
                expected_offset,
                path.display()
            );
        }

        let block_len = seq_len
            .checked_mul(2)
            .and_then(|len| len.checked_add(2))
            .ok_or_else(|| anyhow!("Target '{}' block length overflow", name))?;
        let block_end = offset
            .checked_add(block_len)
            .ok_or_else(|| anyhow!("Target '{}' block offset overflow", name))?;
        if block_end > real_len {
            bail!(
                "Target '{}' block exceeds sequence data (end={}, real_len={}) in {}",
                name,
                block_end,
                real_len,
                path.display()
            );
        }

        offsets.push(offset);
        seq_lens.push(seq_len);
        expected_offset = block_end;
    }

    if expected_offset != real_len {
        bail!(
            "Target blocks cover {} bases but index has {} real bases in {}",
            expected_offset,
            real_len,
            path.display()
        );
    }

    Ok((offsets, seq_lens, real_len))
}

#[cfg(target_endian = "little")]
#[inline]
fn archived_u64_as_native(values: &[rkyv::primitive::ArchivedU64]) -> &[u64] {
    // SAFETY: rkyv's aligned little-endian archived u64 is a transparent-sized,
    // 8-aligned wrapper over native u64 bytes on little-endian targets.
    unsafe { std::slice::from_raw_parts(values.as_ptr().cast::<u64>(), values.len()) }
}

#[cfg(not(target_endian = "little"))]
#[inline]
fn archived_u64_as_native(_values: &[rkyv::primitive::ArchivedU64]) -> &[u64] {
    panic!("mmap-backed target indexes currently require a little-endian target")
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use tempfile::tempdir;

    use super::{
        write_target_registry, TargetRecord, TargetRegistry, TargetStore, TARGET_REGISTRY_VERSION,
    };

    #[test]
    fn roundtrip_build_open_target_view() {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("targets.fa");
        let index_path = dir.path().join("targets.idx");

        let mut fasta = fs_err::File::create(&fasta_path).unwrap();
        writeln!(fasta, ">chrA\nACGUGA").unwrap();
        writeln!(fasta, ">chrB\nUUUGCA").unwrap();
        drop(fasta);

        TargetRegistry::build(&fasta_path, &index_path).unwrap();
        let store = TargetRegistry::open(&index_path).unwrap();
        let expected = [("chrA", 6usize), ("chrB", 6usize)];
        assert_eq!(store.len(), expected.len());

        let target = store.view();
        assert_eq!(target.offsets.len(), store.len());
        assert_eq!(target.seq_lens.len(), store.len());
        assert!(target.sa_real_len > 0);

        for (i, (expected_name, expected_seq_len)) in expected.iter().enumerate().take(store.len())
        {
            let (name, fwd, rc, seq_len) = store.target_seqs(i);
            assert_eq!(name, *expected_name);
            assert_eq!(seq_len, *expected_seq_len);
            assert_eq!(fwd.len(), *expected_seq_len);
            assert_eq!(rc.len(), *expected_seq_len);
        }
    }

    #[test]
    fn open_rejects_invalid_archive() {
        let dir = tempdir().unwrap();
        let bad_path = dir.path().join("bad.idx");
        fs_err::write(&bad_path, b"not-a-valid-rkyv-archive").unwrap();
        match TargetRegistry::open(&bad_path) {
            Ok(_) => panic!("Expected invalid archive error"),
            Err(err) => assert!(err.to_string().contains("Invalid target index archive")),
        }
    }

    #[test]
    fn open_rejects_unsupported_version() {
        let dir = tempdir().unwrap();
        let index_path = dir.path().join("targets.idx");
        let store = tiny_store(999);

        write_target_registry(&index_path, &store).unwrap();

        match TargetRegistry::open(&index_path) {
            Ok(_) => panic!("Expected unsupported version error"),
            Err(err) => assert!(err.to_string().contains("Unsupported target index version")),
        }
    }

    #[test]
    fn open_rejects_missing_sa_padding() {
        let dir = tempdir().unwrap();
        let index_path = dir.path().join("targets.idx");
        let mut store = tiny_store(TARGET_REGISTRY_VERSION);
        store.combined_sa.pop();

        write_target_registry(&index_path, &store).unwrap();

        match TargetRegistry::open(&index_path) {
            Ok(_) => panic!("Expected padding/length error"),
            Err(err) => assert!(err.to_string().contains("sequence/SA length mismatch")),
        }
    }

    #[test]
    fn map_target_pos_normalizes_into_strand_view() {
        use crate::types::Strand;

        let seq_len = 5;
        let seed_len = 2;

        assert_eq!(
            TargetRegistry::map_target_pos(1, seq_len, seed_len),
            Some((Strand::Reverse, 2))
        );
        assert_eq!(
            TargetRegistry::map_target_pos(seq_len + 1 + 2, seq_len, seed_len),
            Some((Strand::Forward, 1))
        );
        assert_eq!(
            TargetRegistry::map_target_pos(seq_len, seq_len, seed_len),
            None
        );
        assert_eq!(
            TargetRegistry::map_target_pos(2 * seq_len + 1, seq_len, seed_len),
            None
        );
    }

    #[test]
    fn global_offsets_are_contiguous() {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("targets.fa");
        let index_path = dir.path().join("targets.idx");

        let mut fasta = fs_err::File::create(&fasta_path).unwrap();
        writeln!(fasta, ">t1\nACGU").unwrap();
        writeln!(fasta, ">t2\nGGCC").unwrap();
        writeln!(fasta, ">t3\nAA").unwrap();
        drop(fasta);

        TargetRegistry::build(&fasta_path, &index_path).unwrap();
        let store = TargetRegistry::open(&index_path).unwrap();
        let target = store.view();

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

    fn tiny_store(version: u32) -> TargetStore {
        let mut combined_seq = vec![1, 0, 1, 0];
        combined_seq.resize(combined_seq.len() + super::SA_CHAR_PADDING, 0);
        let mut combined_sa = vec![0, 1, 2, 3];
        combined_sa.resize(combined_sa.len() + super::SA_CHAR_PADDING, 0);

        TargetStore {
            version,
            targets: vec![TargetRecord {
                name: "t1".to_string(),
                offset: 0,
                seq_len: 1,
            }],
            combined_seq,
            combined_sa,
        }
    }
}
