use std::ops::Range;

use super::{partition, recurse, SeedMatch, SeedSaView, SeedingContext, BASES, SLOTS};

/// Half-singleton: one SA interval has a single entry, the other has multiple.
/// `Q_SINGLETON=true` → query is the single entry, partition target.
/// `Q_SINGLETON=false` → target is the single entry, partition query.
#[inline(always)]
pub(super) fn recurse_half_singleton<
    F: FnMut(SeedMatch),
    const WOBBLE: bool,
    const Q_SINGLETON: bool,
>(
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
pub(super) fn recurse_singleton<F: FnMut(SeedMatch), const WOBBLE: bool>(
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
