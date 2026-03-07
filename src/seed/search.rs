use crate::config::SeedConfig;
use crate::index::store::{GlobalView, SA_CHAR_PADDING};
use crate::registry::Query;
use crate::types::{SeedLen, Strand, TargetId};

use super::searcher::SeedSearcher;
use super::SeedHit;

/// Check if any position in range [start, start+len) contains 'N'
#[inline]
fn has_n_in_range(query: &Query, start: usize, len: usize) -> bool {
    if !query.has_n_any() {
        return false;
    }
    let end = start + len;
    query.n_prefix()[end] != query.n_prefix()[start]
}

/// Map a target position in the global combined target buffer to the
/// strand-selected target view consumed by seeding/DP.
#[inline]
fn map_target_view_start(
    local_pos: usize,
    seq_len: usize,
    seed_len: usize,
) -> Option<(Strand, usize)> {
    if local_pos < seq_len {
        if local_pos + seed_len > seq_len {
            return None;
        }
        return Some((Strand::Reverse, seq_len - (local_pos + seed_len)));
    }

    let reverse_block_start = seq_len + 1;
    let reverse_block_end = reverse_block_start + seq_len;
    if (reverse_block_start..reverse_block_end).contains(&local_pos) {
        let reverse_pos = local_pos - reverse_block_start;
        if reverse_pos + seed_len > seq_len {
            return None;
        }
        return Some((Strand::Forward, seq_len - (reverse_pos + seed_len)));
    }

    None
}

/// Enumerate seeds for a single query against the global index.
///
/// The global combined_seq contains all targets concatenated with Gap separators.
/// After the SeedSearcher finds matches in the global SA, we remap each target
/// position back to a specific target using binary search on the offset table.
pub(crate) fn for_each_seed<F: FnMut(SeedHit)>(
    query: &Query,
    global: &GlobalView<'_>,
    config: &SeedConfig,
    mut on_seed: F,
) {
    let interval = query.seed_interval();
    let q_offset = interval.start;
    let q_end = interval.end;
    let q_seed_len = query.seed_sequence().len();
    let min_len = config
        .seed
        .normalize(query.sequence().len())
        .expect("prepared query must be used with a compatible seed config")
        .2;
    let max_len = q_end.saturating_sub(interval.start);

    let q_sa_real_len = query.sa().len();
    if q_sa_real_len == 0 {
        return;
    }
    let mut padded_q_sa = Vec::with_capacity(q_sa_real_len + SA_CHAR_PADDING);
    padded_q_sa.extend_from_slice(query.sa());
    padded_q_sa.resize(q_sa_real_len + SA_CHAR_PADDING, 0u64);
    let mut padded_q_seq = Vec::with_capacity(q_seed_len + SA_CHAR_PADDING);
    padded_q_seq.extend_from_slice(query.seed_sequence().as_slice());
    padded_q_seq.resize(q_seed_len + SA_CHAR_PADDING, crate::types::Base::Gap);

    let searcher = SeedSearcher::new(
        &padded_q_sa,
        &padded_q_seq,
        0,
        q_sa_real_len,
        global.combined_sa,
        global.combined_seq,
        global.sa_real_len,
        config,
    );

    let offsets = global.offsets;
    let seq_lens = global.seq_lens;
    searcher.for_each_length_range(min_len, max_len, |m| {
        let seed_len = m.seed_len;
        let Some(seed_len_typed) = SeedLen::new(seed_len) else {
            return;
        };

        for &q_sa_pos in &padded_q_sa[m.query_interval.start..m.query_interval.end] {
            let local_q_pos = q_sa_pos as usize;
            if local_q_pos + seed_len > q_seed_len {
                continue;
            }
            let q_pos = q_offset + local_q_pos;
            if has_n_in_range(query, q_pos, seed_len) {
                continue;
            }

            for &t_sa_pos in &global.combined_sa[m.target_interval.start..m.target_interval.end] {
                let t_pos = t_sa_pos as usize;

                // Remap global position to target index via binary search on offsets
                let target_idx = match offsets.partition_point(|&o| o <= t_pos as u64) {
                    0 => continue, // before first target
                    i => i - 1,
                };
                let local_pos = t_pos - offsets[target_idx] as usize;
                let seq_len = seq_lens[target_idx] as usize;

                // Block layout: fwd_comp[seq_len] + Gap + reverse[seq_len] + Gap.
                // Normalize into the strand-selected target view so `target_start`
                // has one stable meaning from `SeedHit` onward.
                let Some((strand, target_start)) =
                    map_target_view_start(local_pos, seq_len, seed_len)
                else {
                    continue;
                };

                on_seed(SeedHit {
                    query_start: q_pos,
                    target_id: TargetId(target_idx as u32),
                    target_start,
                    len: seed_len_typed,
                    strand,
                });
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::map_target_view_start;
    use crate::types::Strand;

    #[test]
    fn target_view_start_is_normalized_into_selected_view() {
        let seq_len = 5;
        let seed_len = 2;

        assert_eq!(
            map_target_view_start(1, seq_len, seed_len),
            Some((Strand::Reverse, 2))
        );
        assert_eq!(
            map_target_view_start(seq_len + 1 + 2, seq_len, seed_len),
            Some((Strand::Forward, 1))
        );
        assert_eq!(map_target_view_start(seq_len, seq_len, seed_len), None);
        assert_eq!(
            map_target_view_start(2 * seq_len + 1, seq_len, seed_len),
            None
        );
    }
}
