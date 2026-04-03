use std::ops::Range;

use crate::types::Base;

use super::{is_valid_base, partition_interval_into, recurse, sa_char, sa_suffix_pos, SeedMatch,
    SeedingContext, BASES, SLOTS};

/// Q-singleton: one query SA entry, multiple target SA entries (partitioned).
#[inline(always)]
pub(super) fn recurse_q_singleton<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut SeedingContext<'_, F>,
    q_idx: usize,
    s: Range<usize>,
    depth: usize,
    match_streak: usize,
    mm_count: usize,
) {
    let q_char = sa_char(ctx.q_sa, ctx.q_seq, q_idx, depth);
    if !is_valid_base(q_char) {
        return;
    }

    let mut sint = [0usize; 6];
    partition_interval_into(ctx.t_sa, ctx.t_seq, s.start, s.end, depth, &mut sint);

    let d1 = depth + 1;
    let can_mm = ctx.max_mm > 0
        && mm_count < ctx.max_mm
        && d1 > ctx.min_prefix
        && match_streak < ctx.min_len
        && ctx.max_len - d1 >= ctx.min_suffix;
    let q_base = Base::from_idx(q_char as usize);
    let q1 = q_idx..q_idx + 1;
    let ms = match_streak + 1;
    let mm1 = mm_count + 1;

    for j in 0..4 {
        let ts = SLOTS[j];
        if sint[ts] >= sint[ts + 1] {
            continue;
        }
        if q_base.pair_type(BASES[j]).is_match(WOBBLE) {
            recurse::<F, WOBBLE>(ctx, q1.clone(), sint[ts]..sint[ts + 1], d1, ms, mm_count);
        } else if can_mm {
            recurse::<F, WOBBLE>(ctx, q1.clone(), sint[ts]..sint[ts + 1], d1, 0, mm1);
        }
    }
}

/// S-singleton: one target SA entry, multiple query SA entries (partitioned).
#[inline(always)]
pub(super) fn recurse_s_singleton<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut SeedingContext<'_, F>,
    q: Range<usize>,
    s_idx: usize,
    depth: usize,
    match_streak: usize,
    mm_count: usize,
) {
    let s_char = sa_char(ctx.t_sa, ctx.t_seq, s_idx, depth);
    if !is_valid_base(s_char) {
        return;
    }

    let mut qint = [0usize; 6];
    partition_interval_into(ctx.q_sa, ctx.q_seq, q.start, q.end, depth, &mut qint);

    let d1 = depth + 1;
    let can_mm = ctx.max_mm > 0
        && mm_count < ctx.max_mm
        && d1 > ctx.min_prefix
        && match_streak < ctx.min_len
        && ctx.max_len - d1 >= ctx.min_suffix;
    let s_base = Base::from_idx(s_char as usize);
    let s1 = s_idx..s_idx + 1;
    let ms = match_streak + 1;
    let mm1 = mm_count + 1;

    for i in 0..4 {
        let qs = SLOTS[i];
        if qint[qs] >= qint[qs + 1] {
            continue;
        }
        if BASES[i].pair_type(s_base).is_match(WOBBLE) {
            recurse::<F, WOBBLE>(ctx, qint[qs]..qint[qs + 1], s1.clone(), d1, ms, mm_count);
        } else if can_mm {
            recurse::<F, WOBBLE>(ctx, qint[qs]..qint[qs + 1], s1.clone(), d1, 0, mm1);
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
