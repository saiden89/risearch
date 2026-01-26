use crate::config::SeedConfig;
use crate::sa::{SaIndexFile, SuffixArray};
use crate::seq::Sequence;
use crate::types::{Base, Strand};

use super::{SeedCandidate, SeedMatch, SeedSearcher};

/// Pre-computed query data to avoid rebuilding per-target.
struct QueryCache {
    /// Normalized query sequence
    pub q_norm: Sequence,
    /// Reverse complement of normalized query
    pub q_rc: Sequence,
    /// Suffix array of reverse complement
    pub q_rc_sa: SuffixArray,
    /// Seed interval start (0-based)
    pub start0: usize,
    /// Seed interval end (1-based, exclusive)
    pub end1: usize,
    /// Minimum seed length
    pub mi_len: usize,
    /// Prefix sum of N positions for O(1) N-checking
    n_prefix: Vec<u32>,
    /// Fast path when query has no Ns
    has_n_any: bool,
}

impl QueryCache {
    /// Build query preprocessing data (called once per query)
    fn new(query: &Sequence, config: &SeedConfig) -> Option<Self> {
        let q_norm = query.clone();
        let q_len = q_norm.len();

        let mut n_prefix = Vec::with_capacity(q_len + 1);
        n_prefix.push(0);
        let mut n_total = 0;
        for &base in q_norm.iter() {
            if base == Base::N {
                n_total += 1;
            }
            n_prefix.push(n_total);
        }
        let has_n_any = n_total != 0;

        // Get seed interval bounds
        let (start1, end1, mi_len) = config.seed.normalize(q_len).ok()?;

        // Build query RC and its SA once
        let q_rc = q_norm.reverse_complement();
        let q_rc_sa = SuffixArray::build(&q_rc);

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

    /// Check if any position in range [start, start+len) contains 'N'
    #[inline]
    fn has_n_in_range(&self, start: usize, len: usize) -> bool {
        if !self.has_n_any {
            return false;
        }
        let end = start + len;
        self.n_prefix[end] != self.n_prefix[start]
    }
}

/// Find all seed matches, reusing provided Vecs to avoid allocation.
///
/// Clears `candidates` and `matches` before filling.
pub(crate) fn find_seeds(
    query: &Sequence,
    index: &SaIndexFile,
    config: &SeedConfig,
    candidates: &mut Vec<SeedCandidate>,
    matches: &mut Vec<SeedMatch>,
) {
    candidates.clear();
    matches.clear();

    // Pre-compute query data once (SA, RC, etc.)
    let Some(prep) = QueryCache::new(query, config) else {
        return;
    };

    for (idx, target) in index.sequences.iter().enumerate() {
        collect_target_seeds(
            &prep,
            idx,
            config,
            candidates,
            matches,
            Strand::Forward,
            &target.forward_sa,
            &target.sequence,
        );
        collect_target_seeds(
            &prep,
            idx,
            config,
            candidates,
            matches,
            Strand::Reverse,
            &target.reverse_sa,
            &target.sequence_rc,
        );
    }
}

/// Collect seeds from a single target strand, reusing the provided scratch buffers.
fn collect_target_seeds(
    prep: &QueryCache,
    target_idx: usize,
    config: &SeedConfig,
    candidates: &mut Vec<SeedCandidate>,
    matches: &mut Vec<SeedMatch>,
    strand: Strand,
    t_sa: &SuffixArray,
    t_seq: &Sequence,
) {
    let q_len = prep.q_norm.len();
    let start0 = prep.start0;
    let end1 = prep.end1;
    let mi_len = prep.mi_len;

    matches.clear();
    let searcher = SeedSearcher::new(&prep.q_rc_sa, &prep.q_rc, t_sa, t_seq, config);
    searcher.search_length_range(mi_len, q_len, matches);

    for m in matches.iter() {
        let seed_len = m.depth;
        for &q_rc_pos_i32 in &prep.q_rc_sa[m.query_interval.start..m.query_interval.end] {
            let q_rc_pos = q_rc_pos_i32 as usize;
            if q_rc_pos + seed_len > q_len {
                continue;
            }
            let q_pos = q_len - q_rc_pos - seed_len;
            if q_pos < start0 || q_pos + seed_len > end1 {
                continue;
            }
            if prep.has_n_in_range(q_pos, seed_len) {
                continue;
            }

            for &t_pos_i32 in &t_sa[m.target_interval.start..m.target_interval.end] {
                let t_pos = t_pos_i32 as usize;
                if t_pos + seed_len > t_seq.len() {
                    continue;
                }
                candidates.push(SeedCandidate {
                    query_pos: q_pos,
                    target_idx,
                    target_start: t_pos,
                    len: seed_len,
                    strand,
                });
            }
        }
    }
}
