//! Suffix array newtype and utilities.

use std::ops::Deref;

use libsais::SuffixArrayConstruction;
use serde::{Deserialize, Serialize};

pub use crate::index::io::{load_index_file, write_index_file};
pub use crate::index::sa::*;
pub use crate::registry::TargetRegistry;
use crate::Sequence;

/// A suffix array stored as Vec<u32>.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct SuffixArray(Vec<u32>);

impl SuffixArray {
    /// Build from a sequence.
    pub fn build(seq: &Sequence) -> Self {
        let bytes = seq.to_bytes();
        Self::try_build(&bytes).expect("SA construction failed")
    }

    /// Build from bytes with error handling.
    pub fn try_build(bytes: &[u8]) -> Result<Self, Box<dyn std::error::Error + Send + Sync>> {
        SuffixArrayConstruction::for_text(bytes)
            .in_owned_buffer()
            .single_threaded()
            .run()
            .map(|sa| Self(sa.into_vec().into_iter().map(|x: i64| x as u32).collect()))
            .map_err(|e| format!("{e:?}").into())
    }

    /// Consume and return inner Vec.
    pub fn into_inner(self) -> Vec<u32> {
        self.0
    }
}

impl Deref for SuffixArray {
    type Target = [u32];
    fn deref(&self) -> &[u32] {
        &self.0
    }
}

impl From<Vec<u32>> for SuffixArray {
    fn from(sa: Vec<u32>) -> Self {
        Self(sa)
    }
}
