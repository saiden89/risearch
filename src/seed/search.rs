use crate::config::SeedConfig;
use crate::registry::{QueryData, TargetRegistry};
use crate::types::{Base, SeedLen, Strand, TargetId};

use super::searcher::SeedSearcher;
use super::SeedHit;

pub(crate) struct TargetSeedView<'a> {
    pub combined_seq: &'a [Base],
    pub combined_sa: &'a [u32],
    pub seq_len: usize,
}

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

    for (idx, target) in index.entries().iter().enumerate() {
        let target_view = TargetSeedView {
            combined_seq: &target.combined_seq,
            combined_sa: &target.combined_sa,
            seq_len: target.seq_len,
        };
        for_each_seed_one_target(query, idx as u32, &target_view, config, |seed| {
            candidates.push(seed);
        });
    }
}

pub(crate) fn for_each_seed_one_target<F: FnMut(SeedHit)>(
    query: &QueryData,
    target_idx: u32,
    target: &TargetSeedView<'_>,
    config: &SeedConfig,
    mut on_seed: F,
) {
    let q_len = query.sequence().len();
    let interval = query.seed_interval();
    let q_start = interval.start;
    let q_end = interval.end;
    let min_len = query.min_seed_len();
    let seq_len = target.seq_len;

    let mut matches = Vec::with_capacity(1024);

    let searcher = SeedSearcher::new(
        query.reverse_sa(),
        query.sequence_rc(),
        target.combined_sa,
        target.combined_seq,
        config,
    );
    searcher.search_length_range(min_len, q_len, &mut matches);

    for m in matches.iter() {
        let seed_len = m.seed_len;
        for &q_rc_pos_i32 in &query.reverse_sa()[m.query_interval.start..m.query_interval.end] {
            let q_rc_pos = q_rc_pos_i32 as usize;
            if q_rc_pos + seed_len > q_len {
                continue;
            }
            let q_pos = q_len - q_rc_pos - seed_len;
            if q_pos < q_start || q_pos + seed_len > q_end {
                continue;
            }
            if has_n_in_range(query, q_pos, seed_len) {
                continue;
            }

            for &t_pos_i32 in &target.combined_sa[m.target_interval.start..m.target_interval.end] {
                let t_pos = t_pos_i32 as usize;

                // Determine strand from position in combined sequence.
                // t_pos == seq_len is the Gap separator — skip it.
                let (strand, target_start) = if t_pos < seq_len {
                    if t_pos + seed_len > seq_len {
                        continue;
                    }
                    (Strand::Forward, t_pos)
                } else if t_pos > seq_len {
                    let rc_pos = t_pos - seq_len - 1;
                    if rc_pos + seed_len > seq_len {
                        continue;
                    }
                    (Strand::Reverse, rc_pos)
                } else {
                    continue; // Gap separator position
                };

                on_seed(SeedHit {
                    query_pos: q_pos,
                    target_id: TargetId(target_idx),
                    target_start,
                    seed_len: SeedLen::new(seed_len)
                        .expect("seed length from search must be positive and fit in u16"),
                    strand,
                });
            }
        }
    }
}
