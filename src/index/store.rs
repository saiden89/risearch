use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use fs_err::File;
use memmap2::Mmap;
use rayon::prelude::*;

use crate::fastx::{normalize_record, read_and_validate_fasta};
use crate::index::sa::{SuffixIndex, SuffixIndexView};
use crate::index::view::TargetView;
use crate::types::{Base, Strand};

const TARGET_REGISTRY_VERSION: u32 = 3;

/// Runtime handle for a target index.
///
/// The mapped file contains a private rkyv `TargetStore`; this type owns the
/// mmap lifetime and carries the native offset directory needed by `TargetView`.
pub struct TargetRegistry {
    mmap: Mmap,
    offsets: Vec<usize>,
}

#[derive(rkyv::Archive, rkyv::Serialize)]
struct TargetStore {
    version: u32,
    targets: Vec<TargetRecord>,
    suffix_index: SuffixIndex,
}

#[derive(rkyv::Archive, rkyv::Serialize)]
struct TargetRecord {
    name: String,
    offset: u64,
}

impl TargetRegistry {
    /// Build a target index from a FASTA file and write it as an rkyv archive.
    ///
    /// `threads` controls suffix-index construction parallelism.
    /// `None` (or `Some(0)`) means auto.
    pub fn build(input: &Path, output: &Path, threads: Option<usize>) -> Result<()> {
        validate_output_path(output)?;

        let records = read_and_validate_fasta(input)?;
        let mut targets = Vec::new();
        let mut combined_bases = Vec::new();

        for (id, raw_seq) in records {
            let Some(sequence) = normalize_record(&id, &raw_seq)? else {
                continue;
            };

            let offset = combined_bases.len() as u64;

            targets.push(TargetRecord { name: id, offset });

            combined_bases.reserve(2 * sequence.len() + 2);
            combined_bases.extend(sequence[..].iter().rev().copied());
            combined_bases.push(Base::Gap);
            combined_bases.extend(sequence.iter().copied().map(Base::complement));
            combined_bases.push(Base::Gap);
        }

        if targets.is_empty() {
            bail!(
                "All sequences were empty after normalization in {}",
                input.display()
            );
        }

        let suffix_index = SuffixIndex::build(combined_bases, threads)
            .context("Failed to build global suffix array")?;

        let store = TargetStore {
            version: TARGET_REGISTRY_VERSION,
            targets,
            suffix_index,
        };

        write_target_registry(output, &store)
    }

    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open index file: {}", path.display()))?;
        let mmap = unsafe { Mmap::map(&file) }
            .with_context(|| format!("Failed to memory-map index file: {}", path.display()))?;

        // An index written by an older version fails bytecheck before the version
        // check below, so the rebuild hint has to be on this error too.
        let root = rkyv::access::<ArchivedTargetStore, rkyv::rancor::Error>(mmap.as_ref())
            .with_context(|| {
                format!(
                    "Invalid target index archive: {}; rebuild it with `risearch index`",
                    path.display()
                )
            })?;
        let offsets = validate_store(root, path)?;

        Ok(Self { mmap, offsets })
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
    pub fn view(&self) -> TargetView<'_> {
        TargetView::new(self.suffixes(), &self.offsets)
    }

    /// Return the selected physical target strand in duplex-column order.
    ///
    /// For input `T` written 5' to 3', Forward is `R(T)` and Reverse is `C(T)`.
    #[inline]
    pub fn target(&self, target_idx: usize, strand: Strand) -> &[Base] {
        self.view().target(target_idx, strand)
    }

    #[inline]
    fn root(&self) -> &ArchivedTargetStore {
        // SAFETY: `open` validates the archive before constructing `TargetRegistry`.
        unsafe { rkyv::access_unchecked::<ArchivedTargetStore>(self.mmap.as_ref()) }
    }

    #[inline]
    fn suffixes(&self) -> SuffixIndexView<'_> {
        let stored = &self.root().suffix_index;
        let bytes = stored.sequence.as_slice();
        let suffix_array = archived_u64_as_native(stored.suffix_array.as_slice());
        // SAFETY: `open` validates every stored byte as a Base discriminant.
        unsafe { SuffixIndexView::from_bytes_unchecked(bytes, suffix_array) }
    }
}

fn validate_output_path(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            let md = fs_err::metadata(parent).with_context(|| {
                format!("Output directory does not exist: {}", parent.display())
            })?;
            if !md.is_dir() {
                bail!("Output parent is not a directory: {}", parent.display());
            }
        }
    }
    Ok(())
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

fn validate_store(root: &ArchivedTargetStore, path: &Path) -> Result<Vec<usize>> {
    if root.version.to_native() != TARGET_REGISTRY_VERSION {
        bail!(
            "Unsupported target index version {} in {}; rebuild it with `risearch index`",
            root.version.to_native(),
            path.display()
        );
    }

    let seq = root.suffix_index.sequence.as_slice();
    let sa = root.suffix_index.suffix_array.as_slice();
    if seq.len() != sa.len() {
        bail!(
            "Target index sequence/SA length mismatch (seq={}, sa={}): {}",
            seq.len(),
            sa.len(),
            path.display()
        );
    }
    let data_len = sa.len();
    // Only the sequence is scanned: reinterpreting these bytes as `Base` needs
    // every discriminant to be valid, while an out-of-range suffix position is
    // already harmless because `base_unchecked` clamps past-the-end reads.
    // Reduced as a per-chunk max rather than a short-circuiting `any`, which
    // cannot vectorize; at index scale that is the difference between ~0.1 and
    // ~0.8 G instructions per open.
    let max_rank = seq
        .par_chunks(1 << 20)
        .map(|chunk| chunk.iter().copied().fold(0u8, u8::max))
        .max()
        .unwrap_or(0);
    if max_rank > Base::U.as_u8() {
        bail!(
            "Target index contains an invalid base rank: {}",
            path.display()
        );
    }
    let mut offsets = Vec::with_capacity(root.targets.len());
    for target in root.targets.iter() {
        let offset = usize::try_from(target.offset.to_native())
            .context("Target offset does not fit in usize")?;
        offsets.push(offset);
    }

    if offsets.first() != Some(&0) {
        bail!(
            "First target block does not start at offset 0 in {}",
            path.display()
        );
    }

    for (target_idx, target) in root.targets.iter().enumerate() {
        let name = target.name.as_str();
        let offset = offsets[target_idx];
        let block_end = offsets.get(target_idx + 1).copied().unwrap_or(data_len);
        if block_end > data_len {
            bail!(
                "Target '{}' block exceeds sequence data (end={}, data_len={}) in {}",
                name,
                block_end,
                data_len,
                path.display()
            );
        }
        let block_len = block_end.checked_sub(offset).ok_or_else(|| {
            anyhow!(
                "Target '{}' offset {} is not before block end {} in {}",
                name,
                offset,
                block_end,
                path.display()
            )
        })?;
        if block_len < 2 || (block_len - 2) % 2 != 0 {
            bail!(
                "Target '{}' has invalid block length {} in {}",
                name,
                block_len,
                path.display()
            );
        }
        let seq_len = (block_len - 2) / 2;
        if seq[offset + seq_len] != Base::Gap.as_u8() || seq[block_end - 1] != Base::Gap.as_u8() {
            bail!(
                "Target '{}' block is missing a separator Gap in {}",
                name,
                path.display()
            );
        }
    }

    Ok(offsets)
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

    use crate::index::sa::SuffixIndex;
    use crate::types::{Base, Strand};

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

        TargetRegistry::build(&fasta_path, &index_path, None).unwrap();
        let store = TargetRegistry::open(&index_path).unwrap();
        let expected = [("chrA", 6usize), ("chrB", 6usize)];
        assert_eq!(store.len(), expected.len());

        let target = store.view();
        assert!(!target.suffixes().is_empty());
        assert_eq!(store.offsets.len(), store.len());
        let forward = target.map_target_range(0, Strand::Forward, 1..4);
        assert_eq!(forward, 2..5);
        assert_eq!(target.map_target_range(0, Strand::Forward, forward), 1..4);
        assert_eq!(target.map_target_range(0, Strand::Reverse, 1..4), 1..4);

        for (i, (expected_name, expected_seq_len)) in expected.iter().enumerate().take(store.len())
        {
            assert_eq!(store.get_name(i), *expected_name);
            for strand in [Strand::Forward, Strand::Reverse] {
                assert_eq!(store.target(i, strand).len(), *expected_seq_len);
            }
        }

        assert_eq!(
            store.target(0, Strand::Forward),
            &[Base::A, Base::G, Base::U, Base::G, Base::C, Base::A]
        );
        assert_eq!(
            store.target(0, Strand::Reverse),
            &[Base::U, Base::G, Base::C, Base::A, Base::C, Base::U]
        );
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
    fn open_rejects_v1_index() {
        let dir = tempdir().unwrap();
        let index_path = dir.path().join("targets.idx");
        let store = tiny_store(1);

        write_target_registry(&index_path, &store).unwrap();

        let Err(err) = TargetRegistry::open(&index_path) else {
            panic!("Expected v1 index rejection");
        };
        let message = err.to_string();
        assert!(message.contains("Unsupported target index version 1"));
        assert!(message.contains("rebuild it with `risearch index"));
    }

    #[test]
    fn open_rejects_v2_index() {
        let dir = tempdir().unwrap();
        let index_path = dir.path().join("targets.idx");
        let store = tiny_store(2);

        write_target_registry(&index_path, &store).unwrap();

        let Err(err) = TargetRegistry::open(&index_path) else {
            panic!("Expected v2 index rejection");
        };
        let message = err.to_string();
        assert!(message.contains("Unsupported target index version 2"));
        assert!(message.contains("rebuild it with `risearch index"));
    }

    #[test]
    fn open_rejects_sequence_sa_length_mismatch() {
        let dir = tempdir().unwrap();
        let index_path = dir.path().join("targets.idx");
        let mut store = tiny_store(TARGET_REGISTRY_VERSION);
        store.suffix_index.suffix_array.pop();

        write_target_registry(&index_path, &store).unwrap();

        match TargetRegistry::open(&index_path) {
            Ok(_) => panic!("Expected padding/length error"),
            Err(err) => assert!(err.to_string().contains("sequence/SA length mismatch")),
        }
    }

    #[test]
    fn open_rejects_invalid_base_rank() {
        let dir = tempdir().unwrap();
        let index_path = dir.path().join("targets.idx");
        let mut store = tiny_store(TARGET_REGISTRY_VERSION);
        store.suffix_index.sequence[0] = Base::U.as_u8() + 1;

        write_target_registry(&index_path, &store).unwrap();

        match TargetRegistry::open(&index_path) {
            Ok(_) => panic!("Expected invalid base rejection"),
            Err(err) => assert!(err.to_string().contains("invalid base rank")),
        }
    }

    #[test]
    fn open_tolerates_out_of_range_suffix_position() {
        let dir = tempdir().unwrap();
        let index_path = dir.path().join("targets.idx");
        let mut store = tiny_store(TARGET_REGISTRY_VERSION);
        let past_end = store.suffix_index.sequence.len() as u64;
        store.suffix_index.suffix_array[0] = past_end;

        write_target_registry(&index_path, &store).unwrap();

        // Bounds are enforced at read time by the clamp in `base_unchecked`, so
        // opening must not pay an O(n) scan of the suffix array to prove it.
        let store = TargetRegistry::open(&index_path).expect("open tolerates a stale position");
        let suffixes = store.view().suffixes();
        assert_eq!(suffixes.suffix_positions(0..1), &[past_end]);
        // SAFETY: index 0 is inside the suffix array.
        assert_eq!(unsafe { suffixes.base_unchecked(0, 0) }, Base::Gap);
    }

    #[test]
    fn open_rejects_missing_separator_gap() {
        // Both halves separately: the mid separator, then the trailing one that
        // makes the last stored byte a Gap.
        for corrupt_idx in [1, 3] {
            let dir = tempdir().unwrap();
            let index_path = dir.path().join("targets.idx");
            let mut store = tiny_store(TARGET_REGISTRY_VERSION);
            store.suffix_index.sequence[corrupt_idx] = Base::A.as_u8();

            write_target_registry(&index_path, &store).unwrap();

            match TargetRegistry::open(&index_path) {
                Ok(_) => panic!("Expected separator Gap rejection for byte {corrupt_idx}"),
                Err(err) => assert!(err.to_string().contains("missing a separator Gap")),
            }
        }
    }

    #[test]
    fn open_rejects_v2_layout_archive() {
        #[derive(rkyv::Archive, rkyv::Serialize)]
        struct TargetRecordV2 {
            name: String,
            offset: u64,
            seq_len: u64,
        }

        #[derive(rkyv::Archive, rkyv::Serialize)]
        struct TargetStoreV2 {
            version: u32,
            targets: Vec<TargetRecordV2>,
            combined_seq: Vec<u8>,
            combined_sa: Vec<u64>,
        }

        let dir = tempdir().unwrap();
        let index_path = dir.path().join("targets.idx");
        // A name over rkyv's inline-string threshold; short names still validate
        // and reach the version check instead.
        let store = TargetStoreV2 {
            version: 2,
            targets: vec![TargetRecordV2 {
                name: "ENSG00000139618.15_transcript".to_string(),
                offset: 0,
                seq_len: 1,
            }],
            combined_seq: vec![1, 0, 1, 0],
            combined_sa: vec![0, 1, 2, 3],
        };
        let bytes = rkyv::to_bytes::<rkyv::rancor::Error>(&store).unwrap();
        fs_err::write(&index_path, bytes.as_slice()).unwrap();

        let Err(err) = TargetRegistry::open(&index_path) else {
            panic!("Expected v2 layout rejection");
        };
        assert!(format!("{err:#}").contains("rebuild it with `risearch index`"));
    }

    #[test]
    fn map_seed_pos_resolves_the_target_block() {
        let sequence = [Base::Gap; 24];
        let offsets = [0, 12];
        let suffix_index = SuffixIndex::build(sequence.to_vec(), Some(1)).unwrap();
        let target = crate::index::view::TargetView::new(suffix_index.view(), &offsets);
        let seed_len = 2;

        for (block_idx, block_start) in offsets.iter().copied().enumerate() {
            assert_eq!(
                target.map_seed_pos(block_start, seed_len),
                Some((block_idx, Strand::Forward, 0))
            );
            assert_eq!(
                target.map_seed_pos(block_start + 6, seed_len),
                Some((block_idx, Strand::Reverse, 0))
            );
            assert_eq!(
                target.map_seed_pos(block_start + 6 + 2, seed_len),
                Some((block_idx, Strand::Reverse, 2))
            );
            assert_eq!(target.map_seed_pos(block_start + 5, seed_len), None);
            assert_eq!(target.map_seed_pos(block_start + 11, seed_len), None);
        }
    }

    #[test]
    fn target_view_maps_positions_in_both_strands() {
        let sequence = [Base::Gap; 12];
        let offsets = [0];
        let suffix_index = SuffixIndex::build(sequence.to_vec(), Some(1)).unwrap();
        let target = crate::index::view::TargetView::new(suffix_index.view(), &offsets);
        let seq_len = 5;
        let seed_len = 2;

        assert_eq!(
            target.map_seed_pos(1, seed_len),
            Some((0, Strand::Forward, 1))
        );
        assert_eq!(
            target.map_seed_pos(seq_len + 1 + 2, seed_len),
            Some((0, Strand::Reverse, 2))
        );
        assert_eq!(target.map_seed_pos(seq_len - 1, seed_len), None);
        assert_eq!(target.map_seed_pos(2 * seq_len, seed_len), None);
        assert_eq!(target.map_seed_pos(seq_len, seed_len), None);
        assert_eq!(target.map_seed_pos(2 * seq_len + 1, seed_len), None);
    }

    #[test]
    fn global_offsets_delimit_target_blocks() {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("targets.fa");
        let index_path = dir.path().join("targets.idx");

        let mut fasta = fs_err::File::create(&fasta_path).unwrap();
        writeln!(fasta, ">t1\nACGU").unwrap();
        writeln!(fasta, ">t2\nGGCC").unwrap();
        writeln!(fasta, ">t3\nAA").unwrap();
        drop(fasta);

        TargetRegistry::build(&fasta_path, &index_path, None).unwrap();
        let store = TargetRegistry::open(&index_path).unwrap();

        assert_eq!(store.offsets, [0, 10, 20]);
    }

    fn tiny_store(version: u32) -> TargetStore {
        TargetStore {
            version,
            targets: vec![TargetRecord {
                name: "t1".to_string(),
                offset: 0,
            }],
            suffix_index: SuffixIndex {
                sequence: vec![1, 0, 1, 0],
                suffix_array: vec![0, 1, 2, 3],
            },
        }
    }
}
