use crate::config::SeedConfig;
use crate::index::sa::SuffixArray;
use crate::registry::{QueryData, TargetRegistry};
use crate::seq::Sequence;
use crate::types::{SeedLen, Strand, TargetId};

use super::searcher::{SeedMatch, SeedSearcher};
use super::SeedHit;

/// Check if any position in range [start, start+len) contains 'N'
#[inline]
fn has_n_in_range(query: &QueryData, start: usize, len: usize) -> bool {
    if !query.has_n_any() {
        return false;
    }
    let end = start + len;
    query.n_prefix()[end] != query.n_prefix()[start]
}

/// Find all seed matches, reusing provided Vec to avoid allocation.
///
/// Clears `candidates` before filling.
pub(crate) fn find_seeds(
    query: &QueryData,
    index: &TargetRegistry,
    config: &SeedConfig,
    candidates: &mut Vec<SeedHit>,
) {
    candidates.clear();

    // Scratch buffer for SeedMatch results - reused across loop iterations
    let mut matches = Vec::with_capacity(1024);

    for (idx, target) in index.entries().iter().enumerate() {
        collect_target_seeds(
            query,
            idx,
            config,
            candidates,
            &mut matches,
            Strand::Forward,
            &target.forward_sa,
            &target.sequence,
        );
        collect_target_seeds(
            query,
            idx,
            config,
            candidates,
            &mut matches,
            Strand::Reverse,
            &target.reverse_sa,
            &target.sequence_rc,
        );
    }
}

/// Collect seeds from a single target strand, reusing the provided scratch buffers.
#[allow(clippy::too_many_arguments)]
fn collect_target_seeds(
    query: &QueryData,
    target_idx: usize,
    config: &SeedConfig,
    candidates: &mut Vec<SeedHit>,
    matches: &mut Vec<SeedMatch>,
    strand: Strand,
    t_sa: &SuffixArray,
    t_seq: &Sequence,
) {
    let q_len = query.sequence().len();
    let interval = query.seed_interval();
    let start = interval.start;
    let end = interval.end;
    let min_len = query.min_seed_len();

    matches.clear();
    let searcher = SeedSearcher::new(
        query.reverse_sa(),
        query.sequence_rc(),
        t_sa,
        t_seq,
        config,
    );
    searcher.search_length_range(min_len, q_len, matches);

    for m in matches.iter() {
        let seed_len = m.seed_len;
        for &q_rc_pos_i32 in &query.reverse_sa()[m.query_interval.start..m.query_interval.end] {
            let q_rc_pos = q_rc_pos_i32 as usize;
            if q_rc_pos + seed_len > q_len {
                continue;
            }
            let q_pos = q_len - q_rc_pos - seed_len;
            if q_pos < start || q_pos + seed_len > end {
                continue;
            }
            if has_n_in_range(query, q_pos, seed_len) {
                continue;
            }

            for &t_pos_i32 in &t_sa[m.target_interval.start..m.target_interval.end] {
                let t_pos = t_pos_i32 as usize;
                if t_pos + seed_len > t_seq.len() {
                    continue;
                }
                candidates.push(SeedHit {
                    query_pos: q_pos,
                    target_id: TargetId(target_idx as u32),
                    target_start: t_pos,
                    seed_len: SeedLen::new(seed_len)
                        .expect("seed length from search must be positive and fit in u16"),
                    strand,
                });
            }
        }
    }
}
