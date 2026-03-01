//! Suffix array backend utilities.

use std::ops::Deref;

use anyhow::{anyhow, Error};
use libsais::SuffixArrayConstruction;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::types::Base;

/// Map Base → SA sort byte for lexicographic ordering.
/// Matches the order used by RIsearch2: `0, a, c, g, n, u`.
// TODO: Unify Base discriminant order with SA lexicographic order so this
// conversion layer is unnecessary and ordering bugs can't drift between modules.
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

/// A suffix array stored as `u64` suffix positions.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct SuffixArray(Vec<u64>);

impl SuffixArray {
    /// Consume and return inner Vec.
    pub fn into_inner(self) -> Vec<u64> {
        self.0
    }
}

impl Deref for SuffixArray {
    type Target = [u64];
    fn deref(&self) -> &[u64] {
        &self.0
    }
}

impl TryFrom<&[Base]> for SuffixArray {
    type Error = Error;

    fn try_from(bases: &[Base]) -> Result<Self, Self::Error> {
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
        let positions: Vec<u64> = sa_vec
            .into_par_iter()
            .map(|pos_i64| {
                debug_assert!(pos_i64 >= 0, "libsais returned negative suffix position");
                pos_i64 as u64
            })
            .collect();

        Ok(Self(positions))
    }
}
