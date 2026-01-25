pub mod alignment;
pub mod sa;
pub(crate) mod search;
pub mod searcher;
pub mod spec;

pub use alignment::build_seed_alignment;
pub use sa::{build_suffix_array, complement_sequence};
pub(crate) use search::find_seeds;
pub use searcher::{SeedMatch, SeedSearcher};
pub use spec::{MismatchSpec, SeedSpec};

#[cfg(test)]
mod spec_tests;

use crate::types::Strand;

/// A candidate seed match found during suffix array search.
#[derive(Debug, Clone)]
pub struct SeedCandidate {
    /// Position in query sequence (0-based)
    pub query_pos: usize,
    /// Index of target sequence in the index
    pub target_idx: usize,
    /// Start position in target sequence (0-based)
    pub target_start: usize,
    /// Length of the seed match
    pub len: usize,
    /// Strand of the match (forward or reverse)
    pub strand: Strand,
}
