//! Suffix array backend utilities.

use std::ops::Deref;

use anyhow::{anyhow, Error};
use libsais::SuffixArrayConstruction;
use serde::{Deserialize, Serialize};

use crate::Sequence;

/// A suffix array stored as Vec<u32>.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct SuffixArray(Vec<u32>);

impl SuffixArray {
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

impl TryFrom<&Sequence> for SuffixArray {
    type Error = Error;

    fn try_from(seq: &Sequence) -> Result<Self, Self::Error> {
        let bytes = seq.to_bytes();
        SuffixArrayConstruction::for_text(bytes.as_slice())
            .in_owned_buffer()
            .single_threaded()
            .run()
            .map(|sa| Self(sa.into_vec().into_iter().map(|x: i64| x as u32).collect()))
            .map_err(|e| anyhow!("{e:?}"))
    }
}
