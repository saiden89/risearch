use crate::config::SeedConfig;
use crate::registry::QueryData;
use crate::sa::{SuffixArray, TargetRegistry};
use crate::seq::Sequence;
use crate::types::Strand;

use super::{SeedHit, SeedMatch, SeedSearcher};

/// Config-dependent view of query data for seed search.
///
/// Borrows immutable `QueryData` and adds seed interval bounds
/// computed from `SeedConfig`.
struct QueryView<'a> {
    /// Reference to the query's immutable data
    data: &'a QueryData,
    /// Seed interval start (0-based)
    start0: usize,
    /// Seed interval end (1-based, exclusive)
    end1: usize,
    /// Minimum seed length
    mi_len: usize,
}

impl<'a> QueryView<'a> {
    /// Build query view with config-dependent interval bounds.
    fn new(data: &'a QueryData, config: &SeedConfig) -> Option<Self> {
        let q_len = data.sequence().len();

        // Get seed interval bounds
        let (start1, end1, mi_len) = config.seed.normalize(q_len).ok()?;

        Some(Self {
            data,
            start0: start1 - 1,
            end1,
            mi_len,
        })
    }

    /// Check if any position in range [start, start+len) contains 'N'
    #[inline]
    fn has_n_in_range(&self, start: usize, len: usize) -> bool {
        if !self.data.has_n_any() {
            return false;
        }
        let end = start + len;
        self.data.n_prefix()[end] != self.data.n_prefix()[start]
    }
}

/// Find all seed matches, reusing provided Vecs to avoid allocation.
///
/// Clears `candidates` and `matches` before filling.
pub(crate) fn find_seeds(
    query: &QueryData,
    index: &TargetRegistry,
    config: &SeedConfig,
    candidates: &mut Vec<SeedHit>,
    matches: &mut Vec<SeedMatch>,
) {
    candidates.clear();
    matches.clear();

    // Pre-compute query data once (SA, RC, etc.)
    let Some(view) = QueryView::new(query, config) else {
        return;
    };

    for (idx, target) in index.entries().iter().enumerate() {
        collect_target_seeds(
            &view,
            idx,
            config,
            candidates,
            matches,
            Strand::Forward,
            &target.forward_sa,
            &target.sequence,
        );
        collect_target_seeds(
            &view,
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
    view: &QueryView,
    target_idx: usize,
    config: &SeedConfig,
    candidates: &mut Vec<SeedHit>,
    matches: &mut Vec<SeedMatch>,
    strand: Strand,
    t_sa: &SuffixArray,
    t_seq: &Sequence,
) {
    let q_len = view.data.sequence().len();
    let start0 = view.start0;
    let end1 = view.end1;
    let mi_len = view.mi_len;

    matches.clear();
    let searcher = SeedSearcher::new(
        view.data.reverse_sa(),
        view.data.sequence_rc(),
        t_sa,
        t_seq,
        config,
    );
    searcher.search_length_range(mi_len, q_len, matches);

    for m in matches.iter() {
        let seed_len = m.depth;
        for &q_rc_pos_i32 in &view.data.reverse_sa()[m.query_interval.start..m.query_interval.end] {
            let q_rc_pos = q_rc_pos_i32 as usize;
            if q_rc_pos + seed_len > q_len {
                continue;
            }
            let q_pos = q_len - q_rc_pos - seed_len;
            if q_pos < start0 || q_pos + seed_len > end1 {
                continue;
            }
            if view.has_n_in_range(q_pos, seed_len) {
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
