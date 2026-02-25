//! Suffix array backend utilities.

use std::ops::Deref;

use anyhow::{anyhow, Error};
use libsais::SuffixArrayConstruction;
use serde::{Deserialize, Serialize};

use crate::types::PackedSaEntry;
use crate::Sequence;

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
        let bytes = seq.to_bytes();
        let sa_raw = SuffixArrayConstruction::for_text(bytes.as_slice())
            .in_owned_buffer()
            .single_threaded()
            .run()
            .map_err(|e| anyhow!("{e:?}"))?;

        let sa_vec: Vec<i64> = sa_raw.into_vec();
        let packed = sa_vec
            .into_iter()
            .enumerate()
            .map(|(i, pos_i64)| {
                let pos = pos_i64 as usize;
                // Pack the character at string position `i` into the entry at index `i`
                let char_at_i = *seq.get(i).unwrap_or(&crate::types::Base::Gap);
                PackedSaEntry::new(pos, char_at_i).raw()
            })
            .collect();

        Ok(Self(packed))
    }
}
