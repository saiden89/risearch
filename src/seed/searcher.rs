//! Parallel Suffix Array seed search.
//!
//! Recursive traversal of two suffix arrays with base-pairing constraints.
//! Singleton fast-paths avoid partition overhead when one or both SA intervals
//! have a single entry.

use crate::config::SeedConfig;
use crate::index::store::{TargetStore, TargetView};
use crate::registry::QueryRegistry;
use crate::types::{Base, SeedLen, TargetId};
use std::ops::Range;

use super::SeedHit;

/// The four matchable RNA bases indexed alongside SLOTS.
const BASES: [Base; 4] = [Base::A, Base::G, Base::C, Base::U];

/// Partition slot for each matchable base in BASES.
/// `partition()` returns `[A=0, G=1, C=2, U=3, N=4, end]`; this mapping
/// addresses only searchable A/G/C/U slots.
const SLOTS: [usize; 4] = [0, 1, 2, 3];

pub(crate) trait SeedSaView: Copy {
    fn sa_real_len(&self) -> usize;
    fn sa_suffix_pos(&self, sa_idx: usize) -> usize;
    fn sa_base(&self, sa_idx: usize, offset: usize) -> Base;
}

impl SeedSaView for (&[u64], &[Base], usize) {
    #[inline(always)]
    fn sa_real_len(&self) -> usize {
        self.2
    }

    #[inline(always)]
    fn sa_suffix_pos(&self, sa_idx: usize) -> usize {
        unsafe { *self.0.get_unchecked(sa_idx) as usize }
    }

    #[inline(always)]
    fn sa_base(&self, sa_idx: usize, offset: usize) -> Base {
        let suffix_pos = self.sa_suffix_pos(sa_idx);
        unsafe { *self.1.get_unchecked(suffix_pos + offset) }
    }
}

/// A seed match found by parallel SA search.
#[derive(Debug, Clone)]
struct SeedMatch {
    query_interval: Range<usize>,
    target_interval: Range<usize>,
    seed_len: usize,
}

/// Find all seeds across all queries × all targets in a single SA traversal.
///
/// Returns non-empty per-query seed buckets ready for parallel extension.
pub fn collect(
    queries: &QueryRegistry,
    targets: &TargetStore,
    config: &SeedConfig,
) -> Vec<(u32, Vec<SeedHit>)> {
    let qview = queries.view();
    let target = targets.view();
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

    let q = (qview.combined_sa, qview.combined_seed_seq, qview.len);
    let t = (target.combined_sa, target.combined_seq, target.sa_real_len);

    let mut seeds_by_query: Vec<Vec<SeedHit>> = (0..queries.len()).map(|_| Vec::new()).collect();
    let mut ctx = SeedingContext {
        q,
        t,
        min_len: global_min,
        max_len: global_max,
        max_mm: config.mismatch.max_mismatches,
        min_prefix: config.mismatch.min_prefix_matches,
        min_suffix: config.mismatch.min_suffix_matches,
        on_match: &mut |m| {
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

                for &t_sa_pos in &target.combined_sa[m.target_interval.start..m.target_interval.end]
                {
                    let Some((ti, t_local_pos)) = remap(target.offsets, t_sa_pos as usize) else {
                        continue;
                    };
                    let Some((strand, target_start)) = TargetView::map_target_pos(
                        t_local_pos,
                        target.seq_lens[ti] as usize,
                        seed_len,
                    ) else {
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
        },
    };

    let q = 0..ctx.q.sa_real_len();
    let s = 0..ctx.t.sa_real_len();
    if config.seed_wobble {
        recurse::<_, true>(&mut ctx, q, s, 0, 0, 0);
    } else {
        recurse::<_, false>(&mut ctx, q, s, 0, 0, 0);
    }

    seeds_by_query
        .into_iter()
        .enumerate()
        .filter_map(|(qi, seeds)| (!seeds.is_empty()).then_some((qi as u32, seeds)))
        .collect()
}

#[cfg(test)]
mod tests;

/// Check if any position in range [start, start+len) contains 'N'.
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

struct SeedingContext<'a, F: FnMut(SeedMatch)> {
    q: (&'a [u64], &'a [Base], usize),
    t: (&'a [u64], &'a [Base], usize),
    min_len: usize,
    max_len: usize,
    max_mm: usize,
    min_prefix: usize,
    min_suffix: usize,
    on_match: &'a mut F,
}

impl<F: FnMut(SeedMatch)> SeedingContext<'_, F> {
    #[inline(always)]
    fn should_emit(&self, depth: usize, match_streak: usize, mm_count: usize) -> bool {
        depth >= self.min_len
            && depth <= self.max_len
            && (mm_count == 0 || (match_streak >= self.min_suffix && match_streak < self.min_len))
    }

    #[inline(always)]
    fn can_reach_suffix(&self, depth: usize, match_streak: usize, mm_count: usize) -> bool {
        mm_count == 0
            || self.min_suffix == 0
            || match_streak + (self.max_len - depth) >= self.min_suffix
    }

    #[inline(always)]
    fn can_mismatch_next(&self, depth: usize, match_streak: usize, mm_count: usize) -> bool {
        let next_depth = depth + 1;
        self.max_mm > 0
            && mm_count < self.max_mm
            && next_depth > self.min_prefix
            && match_streak < self.min_len
            && self.max_len - next_depth >= self.min_suffix
    }
}

fn recurse<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut SeedingContext<'_, F>,
    q: Range<usize>,
    s: Range<usize>,
    depth: usize,
    match_streak: usize,
    mm_count: usize,
) {
    let (ql, qr) = (q.start, q.end);
    let (sl, sr) = (s.start, s.end);

    if ctx.should_emit(depth, match_streak, mm_count) {
        (ctx.on_match)(SeedMatch {
            query_interval: ql..qr,
            target_interval: sl..sr,
            seed_len: depth,
        });
    }

    if depth >= ctx.max_len {
        return;
    }

    if !ctx.can_reach_suffix(depth, match_streak, mm_count) {
        return;
    }

    if qr - ql == 1 && sr - sl == 1 {
        recurse_singleton::<F, WOBBLE>(ctx, ql, sl, depth, match_streak, mm_count);
        return;
    }
    if qr - ql == 1 {
        recurse_half_singleton::<F, WOBBLE, true>(ctx, ql, sl..sr, depth, match_streak, mm_count);
        return;
    }
    if sr - sl == 1 {
        recurse_half_singleton::<F, WOBBLE, false>(ctx, sl, ql..qr, depth, match_streak, mm_count);
        return;
    }

    let qi = partition(ctx.q, ql, qr, depth);
    let si = partition(ctx.t, sl, sr, depth);

    if qi[0] == qr || si[0] == sr {
        return;
    }

    let d1 = depth + 1;
    let ms = match_streak + 1;
    let can_mm = ctx.can_mismatch_next(depth, match_streak, mm_count);

    for i in 0..4 {
        let qs = SLOTS[i];
        if qi[qs] >= qi[qs + 1] {
            continue;
        }

        for j in 0..4 {
            let ts = SLOTS[j];
            if si[ts] >= si[ts + 1] {
                continue;
            }

            if BASES[i].pair_type(BASES[j]).is_match(WOBBLE) {
                recurse::<F, WOBBLE>(
                    ctx,
                    qi[qs]..qi[qs + 1],
                    si[ts]..si[ts + 1],
                    d1,
                    ms,
                    mm_count,
                );
            } else if can_mm {
                recurse::<F, WOBBLE>(
                    ctx,
                    qi[qs]..qi[qs + 1],
                    si[ts]..si[ts + 1],
                    d1,
                    0,
                    mm_count + 1,
                );
            }
        }
    }
}

const LINEAR_PARTITION_CUTOFF: usize = 1024;

/// Base discriminants that mark partition boundaries (sorted by SA order).
const BOUNDARIES: [u8; 5] = [
    Base::A as u8,
    Base::G as u8,
    Base::C as u8,
    Base::U as u8,
    Base::N as u8,
];

/// Partition a sorted SA interval by base character at `depth`.
///
/// Returns 6 boundary positions `[A, C, G, N, U, end]` in suffix-array order.
/// Sub-interval for slot `k` is `bounds[k]..bounds[k+1]`.
/// Callers that search only matchable RNA bases should skip the `N` bucket.
#[inline(always)]
fn partition<V: SeedSaView>(view: V, start: usize, end: usize, depth: usize) -> [usize; 6] {
    if start >= end {
        return [start; 6];
    }

    let mut out = [0usize; 6];
    if end - start <= LINEAR_PARTITION_CUTOFF {
        let mut i = start;
        for (slot, &target) in BOUNDARIES.iter().enumerate() {
            while i < end && (view.sa_base(i, depth) as u8) < target {
                i += 1;
            }
            out[slot] = i;
        }
    } else {
        for (slot, &target) in BOUNDARIES.iter().enumerate() {
            out[slot] = binary_search(view, start, end, depth, target);
        }
    }
    out[5] = end;
    out
}

/// Binary search for leftmost position where character at `depth` >= `target`.
#[inline(always)]
fn binary_search<V: SeedSaView>(
    view: V,
    mut start: usize,
    mut end: usize,
    depth: usize,
    target: u8,
) -> usize {
    let mut half = (end - start) >> 1;
    while start < end {
        let mid = start + half;
        if (view.sa_base(mid, depth) as u8) >= target {
            end = start + half;
        } else {
            start += if half != 0 { half } else { 1 };
        }
        half >>= 1;
    }
    start
}

/// Half-singleton: one SA interval has a single entry, the other has multiple.
/// `Q_SINGLETON=true` → query is the single entry, partition target.
/// `Q_SINGLETON=false` → target is the single entry, partition query.
#[inline(always)]
fn recurse_half_singleton<F: FnMut(SeedMatch), const WOBBLE: bool, const Q_SINGLETON: bool>(
    ctx: &mut SeedingContext<'_, F>,
    singleton_idx: usize,
    multi: Range<usize>,
    depth: usize,
    match_streak: usize,
    mm_count: usize,
) {
    let (single, multi_view) = if Q_SINGLETON {
        (ctx.q, ctx.t)
    } else {
        (ctx.t, ctx.q)
    };

    let single_base = single.sa_base(singleton_idx, depth);
    if !single_base.is_matchable() {
        return;
    }

    let part = partition(multi_view, multi.start, multi.end, depth);

    let d1 = depth + 1;
    let can_mm = ctx.can_mismatch_next(depth, match_streak, mm_count);
    let single_range = singleton_idx..singleton_idx + 1;
    let ms = match_streak + 1;
    let mm1 = mm_count + 1;

    for k in 0..4 {
        let ps = SLOTS[k];
        if part[ps] >= part[ps + 1] {
            continue;
        }
        let part_range = part[ps]..part[ps + 1];
        let (q, s) = if Q_SINGLETON {
            (single_range.clone(), part_range)
        } else {
            (part_range, single_range.clone())
        };
        if single_base.pair_type(BASES[k]).is_match(WOBBLE) {
            recurse::<F, WOBBLE>(ctx, q, s, d1, ms, mm_count);
        } else if can_mm {
            recurse::<F, WOBBLE>(ctx, q, s, d1, 0, mm1);
        }
    }
}

/// Both-singleton fast path: tight linear scan when both SA intervals have one entry.
///
/// Assumes the caller already handled emission at the current depth, then advances
/// linearly and emits only for newly reached depths.
#[inline(always)]
fn recurse_singleton<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut SeedingContext<'_, F>,
    q_idx: usize,
    s_idx: usize,
    mut depth: usize,
    mut match_streak: usize,
    mut mm_count: usize,
) {
    loop {
        if depth >= ctx.max_len {
            return;
        }

        if !ctx.can_reach_suffix(depth, match_streak, mm_count) {
            return;
        }

        let d1 = depth + 1;
        let can_mm = ctx.can_mismatch_next(depth, match_streak, mm_count);

        let q_base = ctx.q.sa_base(q_idx, depth);
        let s_base = ctx.t.sa_base(s_idx, depth);
        if !q_base.is_matchable() || !s_base.is_matchable() {
            return;
        }

        if q_base.pair_type(s_base).is_match(WOBBLE) {
            depth = d1;
            match_streak += 1;
        } else if can_mm {
            depth = d1;
            mm_count += 1;
            match_streak = 0;
        } else {
            return;
        }

        if ctx.should_emit(depth, match_streak, mm_count) {
            (ctx.on_match)(SeedMatch {
                query_interval: q_idx..q_idx + 1,
                target_interval: s_idx..s_idx + 1,
                seed_len: depth,
            });
        }
    }
}
