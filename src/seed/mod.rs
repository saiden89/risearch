pub(crate) mod search;
pub mod searcher;

pub(crate) use search::{for_each_seed_one_target, TargetSeedView};

use crate::types::{SeedLen, Strand, TargetId};

/// A candidate seed match found during suffix array search.
#[derive(Debug, Clone)]
pub struct SeedHit {
    /// Interval-group identifier from seed SA traversal.
    /// Hits emitted from the same group share identical seed sequence content.
    pub group_id: u32,
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
