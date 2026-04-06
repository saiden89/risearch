use crate::config::SeedConfig;
use crate::index::store::{TargetStore, TargetView};
use crate::registry::QueryRegistry;
use crate::types::{SeedLen, TargetId};

use super::searcher::SeedSearcher;
use super::SeedHit;

/// Check if any position in range [start, start+len) contains 'N'
#[inline]
fn has_n_in_range(n_prefix: &[u32], start: usize, len: usize) -> bool {
    n_prefix[start + len] != n_prefix[start]
}

/// Remap a global SA position to an index via binary search on an offset table.
/// Returns `None` if the position falls before the first entry.
#[inline]
fn remap(offsets: &[u64], global_pos: usize) -> Option<(usize, usize)> {
    let idx = offsets
        .partition_point(|&o| o <= global_pos as u64)
        .checked_sub(1)?;
    Some((idx, global_pos - offsets[idx] as usize))
}

/// Find all seeds across all queries × all targets in a single SA traversal.
///
/// Returns non-empty per-query seed buckets ready for parallel extension.
pub(crate) fn collect_seeds(
    queries: &QueryRegistry,
    targets: &TargetStore,
    config: &SeedConfig,
) -> Vec<(u32, Vec<SeedHit>)> {
    let qview = queries.view();
    let target = targets.target_view();
    if qview.len == 0 {
        return Vec::new();
    }

    let mut global_min = usize::MAX;
    let mut global_max = 0usize;
    for q in queries.entries() {
        global_min = global_min.min(q.min_seed_len);
        global_max = global_max.max(q.max_seed_len);
    }

    if global_min > global_max {
        return Vec::new();
    }

    let searcher = SeedSearcher::new(
        (qview.combined_sa, qview.combined_seed_seq, qview.len),
        0,
        (target.combined_sa, target.combined_seq, target.sa_real_len),
        config,
    );

    let mut seeds_by_query: Vec<Vec<SeedHit>> = (0..queries.len()).map(|_| Vec::new()).collect();
    searcher.for_each_length_range(global_min, global_max, |m| {
        let seed_len = m.seed_len;
        let Some(seed_len_typed) = SeedLen::new(seed_len) else {
            return;
        };

        for &q_sa_pos in &qview.combined_sa[m.query_interval.start..m.query_interval.end] {
            let Some((qi, q_local_pos)) = remap(qview.offsets, q_sa_pos as usize) else {
                continue;
            };
            if q_local_pos + seed_len > qview.seed_seq_lens[qi] as usize {
                continue;
            }

            let query = queries.get(qi as u32);
            if seed_len < query.min_seed_len || seed_len > query.max_seed_len {
                continue;
            }
            let q_pos = query.seed_interval.start + q_local_pos;
            if query.has_n_any() && has_n_in_range(query.n_prefix(), q_pos, seed_len) {
                continue;
            }

            for &t_sa_pos in &target.combined_sa[m.target_interval.start..m.target_interval.end] {
                let Some((ti, t_local_pos)) = remap(target.offsets, t_sa_pos as usize) else {
                    continue;
                };
                let Some((strand, target_start)) =
                    TargetView::map_target_pos(t_local_pos, target.seq_lens[ti] as usize, seed_len)
                else {
                    continue;
                };

                seeds_by_query[qi].push(SeedHit {
                    query_idx: qi as u32,
                    query_start: q_pos,
                    target_id: TargetId(ti as u32),
                    target_start,
                    len: seed_len_typed,
                    strand,
                });
            }
        }
    });

    seeds_by_query
        .into_iter()
        .enumerate()
        .filter_map(|(qi, seeds)| (!seeds.is_empty()).then_some((qi as u32, seeds)))
        .collect()
}
