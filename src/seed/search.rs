use crate::config::SeedConfig;
use crate::index::store::{GlobalView, SA_CHAR_PADDING};
use crate::registry::QueryData;
use crate::types::{Base, SeedLen, Strand, TargetId};

use super::searcher::SeedSearcher;
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

/// Build a C-style partial query SA for seed traversal.
///
/// Suffixes that cannot satisfy the minimum seed length within the configured
/// seed interval are rewritten to `pos = q_len` (sentinel / invalid). We then
/// stable-sort by original SA rank marker (`idx`) so invalid suffixes are
/// grouped first, mirroring C's `sa_create_partial_reverse` pre-pruning.
///
fn build_partial_query_sa(
    query_sa: &[u64],
    seed_start: usize,
    seed_end: usize,
    min_len: usize,
) -> Vec<u64> {
    let q_len = query_sa.len();

    let max_valid_start = seed_end.saturating_sub(min_len);
    let invalid_pos = q_len;

    let mut keyed: Vec<(usize, usize)> = Vec::with_capacity(q_len);
    for (i, &sa_pos) in query_sa.iter().enumerate() {
        let pos = sa_pos as usize;
        let valid = pos >= seed_start && pos <= max_valid_start;
        let idx_key = if valid { i + 1 } else { 0 };
        let out_pos = if valid { pos } else { invalid_pos };
        keyed.push((idx_key, out_pos));
    }

    keyed.sort_unstable_by_key(|(idx_key, _)| *idx_key);

    let mut out = Vec::with_capacity(q_len);
    for (_, pos) in keyed {
        out.push(pos as u64);
    }
    out
}

/// Enumerate seeds for a single query against the global index.
///
/// The global combined_seq contains all targets concatenated with Gap separators.
/// After the SeedSearcher finds matches in the global SA, we remap each target
/// position back to a specific target using binary search on the offset table.
pub(crate) fn for_each_seed<F: FnMut(SeedHit)>(
    query: &QueryData,
    global: &GlobalView<'_>,
    config: &SeedConfig,
    mut on_seed: F,
) {
    let q_len = query.sequence().len();
    let q_sa = query.sa();
    let interval = query.seed_interval();
    let q_start = interval.start;
    let q_end = interval.end;
    let min_len = query.min_seed_len();
    let max_len = q_end.saturating_sub(q_start);

    // Build C-style pre-pruned query SA, then pad for unchecked lookup.
    let partial_q_sa = build_partial_query_sa(q_sa, q_start, q_end, min_len);
    let q_sa_real_len = partial_q_sa.len();
    let q_sa_start = partial_q_sa
        .iter()
        .position(|&sa_pos| sa_pos as usize != q_len)
        .unwrap_or(q_sa_real_len);
    if q_sa_start == q_sa_real_len {
        return;
    }
    let mut padded_q_sa = Vec::with_capacity(q_sa_real_len + SA_CHAR_PADDING);
    padded_q_sa.extend_from_slice(&partial_q_sa);
    padded_q_sa.resize(q_sa_real_len + SA_CHAR_PADDING, 0u64);
    let mut padded_q_seq = Vec::with_capacity(q_len + SA_CHAR_PADDING);
    padded_q_seq.extend_from_slice(query.sequence().as_slice());
    padded_q_seq.resize(q_len + SA_CHAR_PADDING, Base::Gap);

    let searcher = SeedSearcher::new(
        &padded_q_sa,
        &padded_q_seq,
        q_sa_start,
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
            let q_pos = q_sa_pos as usize;
            if q_pos + seed_len > q_len {
                continue;
            }
            if q_pos < q_start || q_pos + seed_len > q_end {
                continue;
            }
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

                // Determine strand from local position within target block.
                // Block layout: fwd_comp[seq_len] + Gap + rc_comp[seq_len] + Gap
                let (strand, t_start) = if local_pos < seq_len {
                    if local_pos + seed_len > seq_len {
                        continue;
                    }
                    (Strand::Reverse, local_pos)
                } else if local_pos > seq_len && local_pos < 2 * seq_len + 1 {
                    let rc_pos = local_pos - seq_len - 1;
                    if rc_pos + seed_len > seq_len {
                        continue;
                    }
                    (Strand::Forward, rc_pos)
                } else {
                    continue; // Gap separator position
                };

                on_seed(SeedHit {
                    query_start: q_pos,
                    target_id: TargetId(target_idx as u32),
                    target_start: t_start,
                    len: seed_len_typed,
                    strand,
                });
            }
        }
    });
}
