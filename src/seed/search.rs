use crate::config::SeedConfig;
use crate::sa::{SaIndexFile, SequenceIndex};
use crate::seq::reverse_complement_dna;
use crate::types::Strand;

use super::{build_suffix_array, SeedCandidate, SeedMatch, SeedSearcher};

/// Pre-computed query data to avoid rebuilding per-target.
pub struct QueryPrep {
    /// Normalized query (lowercase DNA, U->T)
    pub q_norm: Vec<u8>,
    /// Reverse complement of normalized query
    pub q_rc: Vec<u8>,
    /// Suffix array of reverse complement
    pub q_rc_sa: Vec<u32>,
    /// Seed interval start (0-based)
    pub start0: usize,
    /// Seed interval end (1-based, exclusive)
    pub end1: usize,
    /// Minimum seed length
    pub mi_len: usize,
    /// Bitmap: has_n[i] = true if position i contains 'n'
    /// Prefix sum of N positions for O(1) N-checking
    n_prefix: Vec<u32>,
    /// Fast path when query has no Ns
    has_n_any: bool,
}

impl QueryPrep {
    /// Build query preprocessing data (called once per query)
    pub fn new(query: &[u8], config: &SeedConfig) -> Option<Self> {
        let q_len = query.len();

        // Normalize query to lowercase DNA (t not u) and build N prefix sums.
        let mut q_norm = Vec::with_capacity(q_len);
        let mut n_prefix = Vec::with_capacity(q_len + 1);
        n_prefix.push(0);
        for &b in query {
            let lower = b.to_ascii_lowercase();
            let norm = if lower == b'u' { b't' } else { lower };
            q_norm.push(norm);
            let last = *n_prefix.last().unwrap();
            n_prefix.push(last + u32::from(norm == b'n'));
        }
        let has_n_any = *n_prefix.last().unwrap() != 0;

        // Get seed interval bounds
        let (start1, end1, mi_len) = config.seed.normalize(q_len).ok()?;

        // Build query RC and its SA once
        let q_rc = reverse_complement_dna(&q_norm);
        let q_rc_sa = build_suffix_array(&q_rc);

        Some(Self {
            q_norm,
            q_rc,
            q_rc_sa,
            start0: start1 - 1,
            end1,
            mi_len,
            n_prefix,
            has_n_any,
        })
    }

    /// Check if any position in range [start, start+len) contains 'n'
    #[inline]
    pub fn contains_n(&self, start: usize, len: usize) -> bool {
        if !self.has_n_any {
            return false;
        }
        let end = start + len;
        self.n_prefix[end] != self.n_prefix[start]
    }
}

/// Find all seed matches between a query and a single target sequence.
///
/// This is the low-level, parallelizable unit of work. Callers can use Rayon
/// to parallelize over queries or targets as needed.
///
/// Appends candidates to the provided Vec (avoids allocation per target).
pub fn find_seeds_in_target_into(
    prep: &QueryPrep,
    target: &SequenceIndex,
    target_idx: usize,
    config: &SeedConfig,
    candidates: &mut Vec<SeedCandidate>,
    matches: &mut Vec<SeedMatch>,
) {
    let q_len = prep.q_norm.len();
    let start0 = prep.start0;
    let end1 = prep.end1;
    let mi_len = prep.mi_len;

    // Forward strand: use pre-built forward_sa
    let t_sa = &target.forward_sa;
    matches.clear(); // Reuse allocation
    let searcher = SeedSearcher::new(&prep.q_rc_sa, &prep.q_rc, t_sa, &target.sequence, config);
    searcher.find_seeds_range_into(mi_len, q_len, matches);

    for m in matches.iter() {
        let seed_len = m.depth;
        for &q_rc_pos_i32 in &prep.q_rc_sa[m.query_interval.start..m.query_interval.end] {
            let q_rc_pos = q_rc_pos_i32 as usize;
            if q_rc_pos + seed_len > q_len {
                continue;
            }
            let q_pos = q_len - q_rc_pos - seed_len;
            // Check: seed must start >= start0 AND end <= end0
            // seed_end = q_pos + seed_len - 1, so check q_pos + seed_len <= end0 + 1 = end1
            if q_pos < start0 || q_pos + seed_len > end1 {
                continue;
            }
            if prep.contains_n(q_pos, seed_len) {
                continue;
            }

            for &t_pos_i32 in &t_sa[m.target_interval.start..m.target_interval.end] {
                let t_pos = t_pos_i32 as usize;
                if t_pos + seed_len > target.sequence.len() {
                    continue;
                }
                candidates.push(SeedCandidate {
                    query_pos: q_pos,
                    target_idx,
                    target_start: t_pos,
                    len: seed_len,
                    strand: Strand::Forward,
                });
            }
        }
    }

    // Reverse strand: use pre-built reverse_sa and sequence_rc
    let t_rc_sa = &target.reverse_sa;
    let t_rc = &target.sequence_rc;
    matches.clear(); // Reuse allocation
    let searcher = SeedSearcher::new(&prep.q_rc_sa, &prep.q_rc, t_rc_sa, t_rc, config);
    searcher.find_seeds_range_into(mi_len, q_len, matches);

    for m in matches.iter() {
        let seed_len = m.depth;
        for &q_rc_pos_i32 in &prep.q_rc_sa[m.query_interval.start..m.query_interval.end] {
            let q_rc_pos = q_rc_pos_i32 as usize;
            if q_rc_pos + seed_len > q_len {
                continue;
            }
            let q_pos = q_len - q_rc_pos - seed_len;
            // Check: seed must start >= start0 AND end <= end0
            if q_pos < start0 || q_pos + seed_len > end1 {
                continue;
            }
            if prep.contains_n(q_pos, seed_len) {
                continue;
            }

            for &t_pos_i32 in &t_rc_sa[m.target_interval.start..m.target_interval.end] {
                let t_pos = t_pos_i32 as usize;
                if t_pos + seed_len > t_rc.len() {
                    continue;
                }
                candidates.push(SeedCandidate {
                    query_pos: q_pos,
                    target_idx,
                    target_start: t_pos,
                    len: seed_len,
                    strand: Strand::Reverse,
                });
            }
        }
    }
}

/// Find all seed matches between a query and a single target sequence.
/// Convenience wrapper that returns a new Vec.
pub fn find_seeds_in_target(
    prep: &QueryPrep,
    target: &SequenceIndex,
    target_idx: usize,
    config: &SeedConfig,
) -> Vec<SeedCandidate> {
    let mut candidates = Vec::new();
    let mut matches = Vec::new();
    find_seeds_in_target_into(prep, target, target_idx, config, &mut candidates, &mut matches);
    candidates
}

/// Find all seed matches between a query and all targets in an index.
///
/// Convenience wrapper that iterates over all targets. For parallelization,
/// use `find_seeds_in_target` directly with Rayon.
pub fn find_seeds(query: &[u8], index: &SaIndexFile, config: &SeedConfig) -> Vec<SeedCandidate> {
    let mut candidates = Vec::new();
    let mut matches = Vec::new();
    find_seeds_into(query, index, config, &mut candidates, &mut matches);
    candidates
}

/// Find all seed matches, reusing provided Vecs to avoid allocation.
///
/// Clears `candidates` and `matches` before filling.
pub fn find_seeds_into(
    query: &[u8],
    index: &SaIndexFile,
    config: &SeedConfig,
    candidates: &mut Vec<SeedCandidate>,
    matches: &mut Vec<SeedMatch>,
) {
    candidates.clear();
    matches.clear();

    // Pre-compute query data once (SA, RC, etc.)
    let Some(prep) = QueryPrep::new(query, config) else {
        return;
    };

    for (idx, target) in index.sequences.iter().enumerate() {
        find_seeds_in_target_into(&prep, target, idx, config, candidates, matches);
    }
}
