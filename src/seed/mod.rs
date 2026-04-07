mod searcher;

pub use searcher::collect;

use crate::types::{SeedLen, Strand, TargetId};

/// A candidate seed match found during suffix array search.
#[derive(Debug, Clone)]
pub struct SeedHit {
    /// Index of query in the registry
    pub query_idx: u32,
    /// Start position in query sequence (0-based)
    pub query_start: usize,
    /// Index of target sequence in the index
    pub target_id: TargetId,
    /// Start position in the strand-selected target view used by seeding/DP (0-based)
    pub target_start: usize,
    /// Length of the seed match
    pub len: SeedLen,
    /// Strand of the match (forward or reverse)
    pub strand: Strand,
}
