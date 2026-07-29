mod engine;
mod parallel_sa;

pub use engine::SeedingEngine;

use std::ops::Range;

use crate::types::{Base, Strand};

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

impl SeedHit {
    /// Whether the match cannot be grown by one more pairing column on either
    /// side. A non-maximal seed is a shorter copy of a longer one, so extending
    /// it only rediscovers the same duplex.
    ///
    /// `query` and `target` are the same views seeding used; `seed_interval`
    /// bounds the query region seeds were drawn from.
    pub fn is_maximal(
        &self,
        query: &[Base],
        target: &[Base],
        seed_interval: &Range<usize>,
        seed_wobble: bool,
    ) -> bool {
        let pairs = |qi: usize, ti: usize| query[qi].pair_type(target[ti]).is_match(seed_wobble);
        let (q_start, t_start, len) = (self.query_start, self.target_start, self.len);

        if q_start > seed_interval.start && t_start > 0 && pairs(q_start - 1, t_start - 1) {
            return false;
        }

        if q_start + len < seed_interval.end
            && t_start + len < target.len()
            && pairs(q_start + len, t_start + len)
        {
            return false;
        }

        true
    }
}
