use std::ops::Range;

use crate::types::Base;

use super::{is_valid_base, partition, recurse, sa_char, sa_suffix_pos, SeedMatch, SeedingContext,
    BASES, SLOTS};

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
    let (single_sa, single_seq, multi_sa, multi_seq) = if Q_SINGLETON {
        (ctx.q_sa, ctx.q_seq, ctx.t_sa, ctx.t_seq)
    } else {
        (ctx.t_sa, ctx.t_seq, ctx.q_sa, ctx.q_seq)
    };

    let single_char = sa_char(single_sa, single_seq, singleton_idx, depth);
    if !is_valid_base(single_char) {
        return;
    }

    let part = partition(multi_sa, multi_seq, multi.start, multi.end, depth);

    let d1 = depth + 1;
    let can_mm = ctx.max_mm > 0
        && mm_count < ctx.max_mm
        && d1 > ctx.min_prefix
        && match_streak < ctx.min_len
        && ctx.max_len - d1 >= ctx.min_suffix;
    let single_base = Base::from_idx(single_char as usize);
    let single_range = singleton_idx..singleton_idx + 1;
    let ms = match_streak + 1;
    let mm1 = mm_count + 1;

    for k in 0..4 {
        let ps = SLOTS[k];
        if part[ps] >= part[ps + 1] {
            continue;
        }
        let part_range = part[ps]..part[ps + 1];
        // pair_type is symmetric — argument order doesn't affect result.
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
/// Avoids all partition overhead — reads characters directly and advances in a loop.
#[inline(always)]
pub(super) fn recurse_singleton<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut SeedingContext<'_, F>,
    q_idx: usize,
    s_idx: usize,
    mut depth: usize,
    mut match_streak: usize,
    mut mm_count: usize,
    mut emit_current: bool,
) {
    let q_suffix_pos = sa_suffix_pos(ctx.q_sa, q_idx);
    let s_suffix_pos = sa_suffix_pos(ctx.t_sa, s_idx);

    loop {
        if emit_current
            && depth >= ctx.min_len
            && depth <= ctx.max_len
            && (mm_count == 0 || (match_streak >= ctx.min_suffix && match_streak < ctx.min_len))
        {
            (ctx.on_match)(SeedMatch {
                query_interval: q_idx..q_idx + 1,
                target_interval: s_idx..s_idx + 1,
                seed_len: depth,
            });
        }
        emit_current = true;

        if depth >= ctx.max_len {
            return;
        }

        if mm_count > 0 && ctx.min_suffix > 0 {
            let max_possible = match_streak + (ctx.max_len - depth);
            if max_possible < ctx.min_suffix {
                return;
            }
        }

        let d1 = depth + 1;
        let can_mm = ctx.max_mm > 0
            && mm_count < ctx.max_mm
            && d1 > ctx.min_prefix
            && match_streak < ctx.min_len
            && ctx.max_len - d1 >= ctx.min_suffix;

        // SAFETY: SA_CHAR_PADDING sentinels guarantee in-bounds access.
        // Read Base directly — #[repr(u8)] means same layout, zero-cost.
        let q_base = unsafe { *ctx.q_seq.get_unchecked(q_suffix_pos + depth) };
        let s_base = unsafe { *ctx.t_seq.get_unchecked(s_suffix_pos + depth) };
        if !q_base.is_matchable() || !s_base.is_matchable() {
            return;
        }

        if q_base.pair_type(s_base).is_match(WOBBLE) {
            depth = d1;
            match_streak += 1;
            continue;
        }
        if can_mm {
            depth = d1;
            mm_count += 1;
            match_streak = 0;
            continue;
        }
        return;
    }
}
