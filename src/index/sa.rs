//! Suffix array backend utilities.

use std::ops::Deref;

use anyhow::{anyhow, Error};
use libsais::SuffixArrayConstruction;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::types::{Base, PackedSaEntry};
use crate::Sequence;

/// Map Base → SA sort byte for lexicographic ordering.
/// Matches the order used by RIsearch2: `0, a, c, g, n, u`.
#[inline(always)]
fn sa_sort_byte(b: Base) -> u8 {
    match b {
        Base::Gap => 0,
        Base::A => b'a',
        Base::C => b'c',
        Base::G => b'g',
        Base::N => b'n',
        Base::U => b'u',
    }
}

/// A suffix array stored as bit-packed Vec<u64>.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct SuffixArray(Vec<u64>);

impl SuffixArray {
    /// Consume and return inner Vec.
    pub fn into_inner(self) -> Vec<u64> {
        self.0
    }

    /// Access as slice of packed entries.
    pub fn as_packed(&self) -> &[PackedSaEntry] {
        // SAFETY: PackedSaEntry is #[repr(transparent)] around u64
        unsafe { std::slice::from_raw_parts(self.0.as_ptr() as *const PackedSaEntry, self.0.len()) }
    }

    /// Build a suffix array directly from a `&[Base]` slice.
    ///
    /// This avoids extra allocations: the sort-byte conversion is done
    /// in-place into a single temporary buffer, and packing is parallelized
    /// via rayon.
    pub fn build_from_bases(bases: &[Base]) -> Result<Self, Error> {
        // Convert Base → sort-byte in a single allocation
        let mut sort_bytes: Vec<u8> = Vec::with_capacity(bases.len());
        // SAFETY: we immediately write all `len` bytes via ptr::write
        unsafe { sort_bytes.set_len(bases.len()) };
        // Parallelize the byte conversion for large sequences
        sort_bytes
            .par_chunks_mut(1 << 20) // 1 MiB chunks
            .enumerate()
            .for_each(|(chunk_idx, chunk)| {
                let base_offset = chunk_idx * (1 << 20);
                for (j, dst) in chunk.iter_mut().enumerate() {
                    *dst = sa_sort_byte(bases[base_offset + j]);
                }
            });

        // Build the suffix array (single-threaded libsais)
        let sa_raw = SuffixArrayConstruction::for_text(sort_bytes.as_slice())
            .in_owned_buffer()
            .single_threaded()
            .run()
            .map_err(|e| anyhow!("{e:?}"))?;

        // Extract the SA vec first (releases borrow on sort_bytes), then free
        let sa_vec: Vec<i64> = sa_raw.into_vec();
        drop(sort_bytes);
        let packed: Vec<u64> = sa_vec
            .into_par_iter()
            .enumerate()
            .map(|(i, pos_i64)| {
                let pos = pos_i64 as usize;
                let char_at_i = bases.get(i).copied().unwrap_or(Base::Gap);
                PackedSaEntry::new(pos, char_at_i).raw()
            })
            .collect();

        Ok(Self(packed))
    }
}

impl Deref for SuffixArray {
    type Target = [u64];
    fn deref(&self) -> &[u64] {
        &self.0
    }
}

impl From<Vec<u64>> for SuffixArray {
    fn from(sa: Vec<u64>) -> Self {
        Self(sa)
    }
}

impl TryFrom<&Sequence> for SuffixArray {
    type Error = Error;

    fn try_from(seq: &Sequence) -> Result<Self, Self::Error> {
        Self::build_from_bases(seq)
    }
}
