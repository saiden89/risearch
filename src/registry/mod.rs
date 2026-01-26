use serde::{Deserialize, Serialize};

use crate::index::sa::SequenceIndex;
use crate::sa::SuffixArray;
use crate::seq::Sequence;
use crate::types::Base;

pub trait RegistryEntry {
    fn name(&self) -> &str;
}

#[derive(Serialize, Deserialize)]
pub struct Registry<T> {
    entries: Vec<T>,
}

impl<T> Registry<T> {
    pub fn new(entries: Vec<T>) -> Self {
        Self { entries }
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn get(&self, idx: u32) -> &T {
        &self.entries[idx as usize]
    }

    pub fn entries(&self) -> &[T] {
        &self.entries
    }

    pub fn into_entries(self) -> Vec<T> {
        self.entries
    }
}

impl<T: RegistryEntry> Registry<T> {
    pub fn get_name(&self, idx: u32) -> &str {
        self.entries[idx as usize].name()
    }

    pub fn index_of(&self, name: &str) -> Option<u32> {
        self.entries
            .iter()
            .position(|e| e.name() == name)
            .map(|i| i as u32)
    }

    pub fn iter(&self) -> impl Iterator<Item = (u32, &T)> {
        self.entries
            .iter()
            .enumerate()
            .map(|(i, e)| (i as u32, e))
    }
}

pub struct QueryEntry {
    pub index: SequenceIndex,
    n_prefix: Vec<u32>,
    has_n_any: bool,
}

impl QueryEntry {
    pub fn new(index: SequenceIndex) -> Self {
        let q_len = index.sequence.len();
        let mut n_prefix = Vec::with_capacity(q_len + 1);
        n_prefix.push(0);
        let mut n_total = 0;
        for &base in index.sequence.iter() {
            if base == Base::N {
                n_total += 1;
            }
            n_prefix.push(n_total);
        }
        let has_n_any = n_total != 0;
        Self {
            index,
            n_prefix,
            has_n_any,
        }
    }

    pub fn name(&self) -> &str {
        &self.index.name
    }

    pub fn sequence(&self) -> &Sequence {
        &self.index.sequence
    }

    pub fn sequence_rc(&self) -> &Sequence {
        &self.index.sequence_rc
    }

    pub fn reverse_sa(&self) -> &SuffixArray {
        &self.index.reverse_sa
    }

    pub fn n_prefix(&self) -> &[u32] {
        &self.n_prefix
    }

    pub fn has_n_any(&self) -> bool {
        self.has_n_any
    }
}

impl RegistryEntry for QueryEntry {
    fn name(&self) -> &str {
        self.name()
    }
}

pub type QueryRegistry = Registry<QueryEntry>;
pub type TargetEntry = SequenceIndex;
pub type TargetRegistry = Registry<TargetEntry>;

impl QueryRegistry {
    pub fn from_indices(indices: Vec<SequenceIndex>) -> Self {
        Self::new(indices.into_iter().map(QueryEntry::new).collect())
    }

    pub fn from_names(names: Vec<String>) -> Self {
        let entries = names
            .into_iter()
            .map(|name| QueryEntry::new(SequenceIndex {
                name,
                forward_sa: SuffixArray::from(Vec::new()),
                reverse_sa: SuffixArray::from(Vec::new()),
                sequence: Sequence::from(Vec::new()),
                sequence_rc: Sequence::from(Vec::new()),
            }))
            .collect();
        Self { entries }
    }
}

impl RegistryEntry for SequenceIndex {
    fn name(&self) -> &str {
        &self.name
    }
}

impl TargetRegistry {
    pub fn get_sequence(&self, seq_idx: usize) -> &Sequence {
        &self.entries[seq_idx].sequence
    }

    pub fn get_sequence_rc(&self, seq_idx: usize) -> &Sequence {
        &self.entries[seq_idx].sequence_rc
    }

    pub fn get_sequence_len(&self, seq_idx: usize) -> usize {
        self.entries[seq_idx].sequence.len()
    }

    pub fn from_names(names: Vec<String>) -> Self {
        let entries = names
            .into_iter()
            .map(|name| SequenceIndex {
                name,
                forward_sa: SuffixArray::from(Vec::new()),
                reverse_sa: SuffixArray::from(Vec::new()),
                sequence: Sequence::from(Vec::new()),
                sequence_rc: Sequence::from(Vec::new()),
            })
            .collect();
        Self { entries }
    }
}
