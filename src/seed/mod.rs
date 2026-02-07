pub(crate) mod search;
pub(crate) mod searcher;

pub(crate) use search::find_seeds;

use crate::types::{SeedLen, Strand, TargetId};

/// A candidate seed match found during suffix array search.
#[derive(Debug, Clone)]
pub struct SeedHit {
    /// Position in query sequence (0-based)
    pub query_pos: usize,
    /// Index of target sequence in the index
    pub target_id: TargetId,
    /// Start position in target sequence (0-based)
    pub target_start: usize,
    /// Length of the seed match
    pub seed_len: SeedLen,
    /// Strand of the match (forward or reverse)
    pub strand: Strand,
}
