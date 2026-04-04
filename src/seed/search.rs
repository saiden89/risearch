use crate::config::SeedConfig;
use crate::index::store::TargetView;
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
    let idx = offsets.partition_point(|&o| o <= global_pos as u64).checked_sub(1)?;
    Some((idx, global_pos - offsets[idx] as usize))
}


/// Pre-computed per-query seed length bounds, avoiding repeated `normalize()` calls.
struct QuerySeedBounds {
    min_len: usize,
    max_len: usize,
}

/// Find all seeds across all queries × all targets in a single SA traversal.
///
/// Returns non-empty per-query seed buckets ready for parallel extension.
pub(crate) fn collect_seeds(
    queries: &QueryRegistry,
    target: &TargetView<'_>,
    config: &SeedConfig,
) -> Vec<(u32, Vec<SeedHit>)> {
    let qv = queries.query_view();
    if qv.sa_real_len == 0 {
        return Vec::new();
    }

    // Pre-compute per-query seed bounds once.
    let mut global_min = usize::MAX;
    let mut global_max = 0usize;
    let mut bounds = Vec::with_capacity(queries.len());
    for q in queries.entries() {
        let min_len = config.seed.normalize(q.sequence().len()).expect("prepared query").2;
        let max_len = q.seed_interval().end.saturating_sub(q.seed_interval().start);
        global_min = global_min.min(min_len);
        global_max = global_max.max(max_len);
        bounds.push(QuerySeedBounds { min_len, max_len });
    }

    if global_min > global_max {
        return Vec::new();
    }

    let searcher = SeedSearcher::new(
        (qv.combined_sa, qv.combined_seed_seq, qv.sa_real_len),
        0,
        (target.combined_sa, target.combined_seq, target.sa_real_len),
        config,
    );

    let mut seeds_by_query: Vec<Vec<SeedHit>> =
        (0..queries.len()).map(|_| Vec::new()).collect();
    searcher.for_each_length_range(global_min, global_max, |m| {
        let seed_len = m.seed_len;
        let Some(seed_len_typed) = SeedLen::new(seed_len) else {
            return;
        };

        for &q_sa_pos in &qv.combined_sa[m.query_interval.start..m.query_interval.end] {
            let Some((qi, q_local_pos)) = remap(qv.offsets, q_sa_pos as usize) else {
                continue;
            };
            if q_local_pos + seed_len > qv.seed_seq_lens[qi] as usize {
                continue;
            }
            let qb = &bounds[qi];
            if seed_len < qb.min_len || seed_len > qb.max_len {
                continue;
            }

            let query = queries.get(qi as u32);
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

