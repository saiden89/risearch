use crate::config::SeedConfig;
use crate::registry::QueryData;
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

pub(crate) fn for_each_seed_one_target<F: FnMut(SeedHit)>(
    query: &QueryData,
    target_idx: u32,
    target: &TargetSeedView<'_>,
    config: &SeedConfig,
    mut on_seed: F,
) {
    let q_len = query.sequence().len();
    let q_sa = query.reverse_sa();
    let target_id = TargetId(target_idx);
    let interval = query.seed_interval();
    let q_start = interval.start;
    let q_end = interval.end;
    let min_len = query.min_seed_len();
    let max_len = q_end.saturating_sub(q_start);
    let seq_len = target.seq_len;

    let searcher = SeedSearcher::new(
        q_sa,
        query.sequence_rc(),
        target.combined_sa,
        target.combined_seq,
        config,
    );
    searcher.for_each_length_range(min_len, max_len, |m| {
        let seed_len = m.seed_len;
        let seed_len_typed = SeedLen::new(seed_len)
            .expect("seed length from search must be positive and fit in u16");
        for &q_rc_pos_i32 in &q_sa[m.query_interval.start..m.query_interval.end] {
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
                    target_id,
                    target_start,
                    seed_len: seed_len_typed,
                    strand,
                });
            }
        }
    });
}
