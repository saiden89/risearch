mod engine;
mod parallel_sa;

pub use engine::SeedingEngine;

use crate::types::Strand;

/// A candidate seed match found during suffix array search.
#[derive(Debug, Clone)]
pub struct SeedHit {
    /// Index of query in the registry
    pub query_idx: usize,
    /// Start position in query sequence (0-based)
    pub query_start: usize,
    /// Index of target sequence in the index
    pub target_idx: usize,
    /// Start position in the strand-selected target view used by seeding/DP (0-based)
    pub target_start: usize,
    /// Length of the seed match
    pub len: usize,
    /// Strand of the match (forward or reverse)
    pub strand: Strand,
}
