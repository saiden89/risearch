use std::ops::Range;

use super::{
    partition_interval_into, recurse_if_nonempty, sa_char, sa_suffix_pos, SeedMatch,
    SeedingContext, BASE_A, BASE_C, BASE_G, BASE_U,
};

/// Check if a query base pairs with a target base (in complement-transformed space).
#[inline(always)]
fn is_match_pair<const WOBBLE: bool>(q_char: u8, s_char: u8) -> bool {
    match q_char {
        BASE_A => s_char == BASE_U,
        BASE_C => s_char == BASE_G,
        BASE_G => s_char == BASE_C || (WOBBLE && s_char == BASE_U),
        BASE_U => s_char == BASE_A || (WOBBLE && s_char == BASE_G),
        _ => false,
    }
}

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
    if !(BASE_A..=BASE_U).contains(&q_char) {
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

    let q1 = q_idx..q_idx + 1;
    let ms = match_streak + 1;

    // Match branches — direct dispatch on q_char.
    match q_char {
        BASE_A => recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[4]..sint[5], d1, ms, mm_count),
        BASE_C => recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[2]..sint[3], d1, ms, mm_count),
        BASE_G => {
            recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[1]..sint[2], d1, ms, mm_count);
            if WOBBLE {
                recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[4]..sint[5], d1, ms, mm_count);
            }
        }
        BASE_U => {
            recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[0]..sint[1], d1, ms, mm_count);
            if WOBBLE {
                recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[2]..sint[3], d1, ms, mm_count);
            }
        }
        _ => {} // N: no matches
    }

    if !can_mm {
        return;
    }

    let mm1 = mm_count + 1;

    // Mismatch branches — all target slots not paired with q_char.
    match q_char {
        BASE_A => {
            // A matches U: mismatches A, C, G
            recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[0]..sint[1], d1, 0, mm1);
            recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[1]..sint[2], d1, 0, mm1);
            recurse_if_nonempty::<F, WOBBLE>(ctx, q1, sint[2]..sint[3], d1, 0, mm1);
        }
        BASE_C => {
            // C matches G: mismatches A, C, U
            recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[0]..sint[1], d1, 0, mm1);
            recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[1]..sint[2], d1, 0, mm1);
            recurse_if_nonempty::<F, WOBBLE>(ctx, q1, sint[4]..sint[5], d1, 0, mm1);
        }
        BASE_G => {
            // G matches C (+ wobble U): mismatches A, G, and U when no wobble
            recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[0]..sint[1], d1, 0, mm1);
            recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[2]..sint[3], d1, 0, mm1);
            if !WOBBLE {
                recurse_if_nonempty::<F, WOBBLE>(ctx, q1, sint[4]..sint[5], d1, 0, mm1);
            }
        }
        BASE_U => {
            // U matches A (+ wobble G): mismatches G when no wobble, C, U
            if !WOBBLE {
                recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[2]..sint[3], d1, 0, mm1);
            }
            recurse_if_nonempty::<F, WOBBLE>(ctx, q1.clone(), sint[1]..sint[2], d1, 0, mm1);
            recurse_if_nonempty::<F, WOBBLE>(ctx, q1, sint[4]..sint[5], d1, 0, mm1);
        }
        _ => {} // N: no mismatch handling
    }
}

/// S-singleton: one target SA entry, multiple query SA entries (partitioned).
///
/// When the target base is N, it cannot match any query base but still counts
/// as a mismatch opportunity for all non-empty query slots.
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
    if !(BASE_A..=BASE_U).contains(&s_char) {
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

    let s1 = s_idx..s_idx + 1;
    let ms = match_streak + 1;

    // Match branches — dispatch on s_char to select the matching query slot(s).
    if qint[0] < qint[1] && is_match_pair::<WOBBLE>(BASE_A, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qint[0]..qint[1], s1.clone(), d1, ms, mm_count);
    }
    if qint[1] < qint[2] && is_match_pair::<WOBBLE>(BASE_C, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qint[1]..qint[2], s1.clone(), d1, ms, mm_count);
    }
    if qint[2] < qint[3] && is_match_pair::<WOBBLE>(BASE_G, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qint[2]..qint[3], s1.clone(), d1, ms, mm_count);
    }
    if qint[4] < qint[5] && is_match_pair::<WOBBLE>(BASE_U, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qint[4]..qint[5], s1.clone(), d1, ms, mm_count);
    }

    if !can_mm {
        return;
    }

    let mm1 = mm_count + 1;

    // Mismatch branches.
    if qint[0] < qint[1] && !is_match_pair::<WOBBLE>(BASE_A, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qint[0]..qint[1], s1.clone(), d1, 0, mm1);
    }
    if qint[1] < qint[2] && !is_match_pair::<WOBBLE>(BASE_C, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qint[1]..qint[2], s1.clone(), d1, 0, mm1);
    }
    if qint[2] < qint[3] && !is_match_pair::<WOBBLE>(BASE_G, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qint[2]..qint[3], s1.clone(), d1, 0, mm1);
    }
    if qint[4] < qint[5] && !is_match_pair::<WOBBLE>(BASE_U, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qint[4]..qint[5], s1.clone(), d1, 0, mm1);
    }
}

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

        // SAFETY: SA_CHAR_PADDING sentinels guarantee in-bounds access
        let q_char = unsafe { *ctx.q_seq.get_unchecked(q_suffix_pos + depth) as u8 };
        let s_char = unsafe { *ctx.t_seq.get_unchecked(s_suffix_pos + depth) as u8 };
        if !(BASE_A..=BASE_U).contains(&q_char) || !(BASE_A..=BASE_U).contains(&s_char) {
            return;
        }

        if is_match_pair::<WOBBLE>(q_char, s_char) {
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
