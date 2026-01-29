use crate::config::SeedConfig;
use crate::registry::QueryEntry;
use crate::sa::{SuffixArray, TargetRegistry};
use crate::seq::Sequence;
use crate::types::Strand;

use super::{SeedHit, SeedMatch, SeedSearcher};

/// Pre-computed query data to avoid rebuilding per-target.
struct QueryCache<'a> {
    /// Normalized query sequence
    pub q_norm: &'a Sequence,
    /// Reverse complement of normalized query
    pub q_rc: &'a Sequence,
    /// Suffix array of reverse complement
    pub q_rc_sa: &'a SuffixArray,
    /// Seed interval start (0-based)
    pub start0: usize,
    /// Seed interval end (1-based, exclusive)
    pub end1: usize,
    /// Minimum seed length
    pub mi_len: usize,
    /// Prefix sum of N positions for O(1) N-checking
    n_prefix: &'a [u32],
    /// Fast path when query has no Ns
    has_n_any: bool,
}

impl<'a> QueryCache<'a> {
    /// Build query preprocessing data (called once per query)
    fn new(query: &'a QueryEntry, config: &SeedConfig) -> Option<Self> {
        let q_norm = query.sequence();
        let q_len = q_norm.len();

        let n_prefix = query.n_prefix();
        let has_n_any = query.has_n_any();

        // Get seed interval bounds
        let (start1, end1, mi_len) = config.seed.normalize(q_len).ok()?;

        // Build query RC and its SA once
        let q_rc = query.sequence_rc();
        let q_rc_sa = query.reverse_sa();

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
    query: &QueryEntry,
    index: &TargetRegistry,
    config: &SeedConfig,
    candidates: &mut Vec<SeedHit>,
    matches: &mut Vec<SeedMatch>,
) {
    candidates.clear();
    matches.clear();

    // Pre-compute query data once (SA, RC, etc.)
    let Some(prep) = QueryCache::new(query, config) else {
        return;
    };

    for (idx, target) in index.entries().iter().enumerate() {
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
#[allow(clippy::too_many_arguments)]
fn collect_target_seeds(
    query: &QueryCache,
    target_idx: usize,
    config: &SeedConfig,
    candidates: &mut Vec<SeedHit>,
    matches: &mut Vec<SeedMatch>,
    strand: Strand,
    t_sa: &SuffixArray,
    t_seq: &Sequence,
) {
    let q_len = query.q_norm.len();
    let start0 = query.start0;
    let end1 = query.end1;
    let mi_len = query.mi_len;

    matches.clear();
    let searcher = SeedSearcher::new(query.q_rc_sa, query.q_rc, t_sa, t_seq, config);
    searcher.search_length_range(mi_len, q_len, matches);

    for m in matches.iter() {
        let seed_len = m.depth;
        for &q_rc_pos_i32 in &query.q_rc_sa[m.query_interval.start..m.query_interval.end] {
            let q_rc_pos = q_rc_pos_i32 as usize;
            if q_rc_pos + seed_len > q_len {
                continue;
            }
            let q_pos = q_len - q_rc_pos - seed_len;
            if q_pos < start0 || q_pos + seed_len > end1 {
                continue;
            }
            if query.has_n_in_range(q_pos, seed_len) {
                continue;
            }

            for &t_pos_i32 in &t_sa[m.target_interval.start..m.target_interval.end] {
                let t_pos = t_pos_i32 as usize;
                if t_pos + seed_len > t_seq.len() {
                    continue;
                }
                candidates.push(SeedHit {
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
