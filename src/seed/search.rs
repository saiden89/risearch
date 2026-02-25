use crate::config::SeedConfig;
use crate::index::store::SA_CHAR_PADDING;
use crate::registry::QueryData;
use crate::types::{Base, PackedSaEntry, SeedLen, Strand, TargetId};

use super::searcher::SeedSearcher;
use super::SeedHit;

pub(crate) struct TargetSeedView<'a> {
    pub combined_seq: &'a [Base],
    pub combined_sa: &'a [u64],
    pub sa_real_len: usize,
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

/// Build a C-style partial query SA for seed traversal.
///
/// Suffixes that cannot satisfy the minimum seed length within the configured
/// seed interval are rewritten to `pos = q_len` (sentinel / invalid). We then
/// stable-sort by original SA rank marker (`idx`) so invalid suffixes are
/// grouped first, mirroring C's `sa_create_partial_reverse` pre-pruning.
///
/// The packed base nibble is rewritten after sorting so `sa[pos+offset]`
/// character lookups remain valid for binary partitioning.
fn build_partial_query_sa(
    query_seq: &[Base],
    query_sa: &[u64],
    seed_start: usize,
    seed_end: usize,
    min_len: usize,
) -> Vec<u64> {
    let q_len = query_seq.len();
    debug_assert_eq!(query_sa.len(), q_len);

    let max_valid_start = seed_end.saturating_sub(min_len);
    let invalid_pos = q_len;

    let mut keyed: Vec<(usize, usize)> = Vec::with_capacity(q_len);
    for (i, &packed) in query_sa.iter().enumerate() {
        let pos = PackedSaEntry::from(packed).pos();
        let valid = pos >= seed_start && pos <= max_valid_start;
        let idx_key = if valid { i + 1 } else { 0 };
        let out_pos = if valid { pos } else { invalid_pos };
        keyed.push((idx_key, out_pos));
    }

    keyed.sort_unstable_by_key(|(idx_key, _)| *idx_key);

    let mut out = Vec::with_capacity(q_len);
    for (i, (_, pos)) in keyed.into_iter().enumerate() {
        out.push(PackedSaEntry::new(pos, query_seq[i]).raw());
    }
    out
}

pub(crate) fn for_each_seed_one_target<F: FnMut(SeedHit)>(
    query: &QueryData,
    target_idx: u32,
    target: &TargetSeedView<'_>,
    config: &SeedConfig,
    mut on_seed: F,
) {
    let q_len = query.sequence().len();
    let q_sa = query.sa();
    let target_id = TargetId(target_idx);
    let interval = query.seed_interval();
    let q_start = interval.start;
    let q_end = interval.end;
    let min_len = query.min_seed_len();
    let max_len = q_end.saturating_sub(q_start);
    let seq_len = target.seq_len;

    // Build C-style pre-pruned query SA, then pad for unchecked lookup.
    // Query SAs are small, so this per-query setup is cheap and removes large
    // portions of the recursion tree in mismatch-heavy searches.
    let partial_q_sa = build_partial_query_sa(query.sequence(), q_sa, q_start, q_end, min_len);
    let q_sa_real_len = partial_q_sa.len();
    let q_sa_start = partial_q_sa
        .iter()
        .position(|&packed| PackedSaEntry::from(packed).pos() != q_len)
        .unwrap_or(q_sa_real_len);
    if q_sa_start == q_sa_real_len {
        return;
    }
    let mut padded_q_sa = Vec::with_capacity(q_sa_real_len + SA_CHAR_PADDING);
    padded_q_sa.extend_from_slice(&partial_q_sa);
    padded_q_sa.resize(q_sa_real_len + SA_CHAR_PADDING, 0u64);
    let mut padded_q_seq = Vec::with_capacity(q_len + SA_CHAR_PADDING);
    padded_q_seq.extend_from_slice(query.sequence());
    padded_q_seq.resize(q_len + SA_CHAR_PADDING, Base::Gap);

    let searcher = SeedSearcher::new(
        &padded_q_sa,
        &padded_q_seq,
        q_sa_start,
        q_sa_real_len,
        target.combined_sa,
        target.combined_seq,
        target.sa_real_len,
        config,
    );
    let mut group_id: u32 = 0;
    searcher.for_each_length_range(min_len, max_len, |m| {
        group_id = group_id.wrapping_add(1);
        let this_group = group_id;
        let seed_len = m.seed_len;
        let seed_len_typed = SeedLen::new(seed_len)
            .expect("seed length from search must be positive and fit in u16");
        for &q_packed in &padded_q_sa[m.query_interval.start..m.query_interval.end] {
            let q_pos = PackedSaEntry::from(q_packed).pos();
            if q_pos + seed_len > q_len {
                continue;
            }
            if q_pos < q_start || q_pos + seed_len > q_end {
                continue;
            }
            if has_n_in_range(query, q_pos, seed_len) {
                continue;
            }

            for &t_packed in &target.combined_sa[m.target_interval.start..m.target_interval.end] {
                let t_pos = PackedSaEntry::from(t_packed).pos();

                // Check maximality before doing anything else
                // if !is_maximal(
                //     query.sequence(),
                //     target.combined_seq,
                //     q_pos,
                //     t_pos,
                //     seed_len,
                //     pair_table,
                // ) {
                //     continue;
                // }

                // Determine strand from position in combined sequence.
                // t_pos == seq_len is the Gap separator — skip it.
                let (strand, target_start) = if t_pos < seq_len {
                    if t_pos + seed_len > seq_len {
                        continue;
                    }
                    (Strand::Reverse, t_pos)
                } else if t_pos > seq_len {
                    let rc_pos = t_pos - seq_len - 1;
                    if rc_pos + seed_len > seq_len {
                        continue;
                    }
                    (Strand::Forward, rc_pos)
                } else {
                    continue; // Gap separator position
                };

                on_seed(SeedHit {
                    group_id: this_group,
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
