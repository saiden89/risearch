use crate::types::Interval;

use super::{
    partition_interval_into, recurse_if_nonempty, sa_char, sa_suffix_pos, SeedMatch,
    SeedingContext,
    BASE_A, BASE_C, BASE_G, BASE_U,
};

/// Check if a query base pairs with a target base (in complement-transformed space).
#[inline(always)]
fn is_match_pair<const WOBBLE: bool>(q_char: u8, s_char: u8) -> bool {
    match q_char {
        BASE_A => s_char == BASE_U,
        BASE_G => s_char == BASE_C || (WOBBLE && s_char == BASE_U),
        BASE_C => s_char == BASE_G,
        BASE_U => s_char == BASE_A || (WOBBLE && s_char == BASE_G),
        _ => false,
    }
}

#[inline(always)]
pub(super) fn recurse_q_singleton<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut SeedingContext<'_, F>,
    q_idx: usize,
    sl: usize,
    sr: usize,
    depth: usize,
    match_streak: usize,
    mm_count: usize,
) {
    let q_char = sa_char(ctx.q_sa, ctx.q_seq, q_idx, depth);
    if !(BASE_A..=BASE_U).contains(&q_char) {
        return;
    }

    let mut sint = [0usize; 6];
    partition_interval_into(ctx.t_sa, ctx.t_seq, sl, sr, depth, &mut sint);

    let d1 = depth + 1;
    let can_mm = ctx.max_mm > 0
        && mm_count < ctx.max_mm
        && d1 > ctx.min_prefix
        && ctx.max_len - d1 >= ctx.min_suffix;

    let (sa_lo, sa_hi) = (sint[0], sint[1]);
    let (sc_lo, sc_hi) = (sint[1], sint[2]);
    let (sg_lo, sg_hi) = (sint[2], sint[3]);
    let (su_lo, su_hi) = (sint[4], sint[5]);

    match q_char {
        BASE_A => recurse_if_nonempty::<F, WOBBLE>(
            ctx,
            q_idx,
            q_idx + 1,
            su_lo,
            su_hi,
            d1,
            match_streak + 1,
            mm_count,
        ),
        BASE_G => {
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                q_idx,
                q_idx + 1,
                sc_lo,
                sc_hi,
                d1,
                match_streak + 1,
                mm_count,
            );
            if WOBBLE {
                recurse_if_nonempty::<F, WOBBLE>(
                    ctx,
                    q_idx,
                    q_idx + 1,
                    su_lo,
                    su_hi,
                    d1,
                    match_streak + 1,
                    mm_count,
                );
            }
        }
        BASE_C => recurse_if_nonempty::<F, WOBBLE>(
            ctx,
            q_idx,
            q_idx + 1,
            sg_lo,
            sg_hi,
            d1,
            match_streak + 1,
            mm_count,
        ),
        BASE_U => {
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                q_idx,
                q_idx + 1,
                sa_lo,
                sa_hi,
                d1,
                match_streak + 1,
                mm_count,
            );
            if WOBBLE {
                recurse_if_nonempty::<F, WOBBLE>(
                    ctx,
                    q_idx,
                    q_idx + 1,
                    sg_lo,
                    sg_hi,
                    d1,
                    match_streak + 1,
                    mm_count,
                );
            }
        }
        _ => {}
    }

    if !can_mm {
        return;
    }

    match q_char {
        BASE_A => {
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                q_idx,
                q_idx + 1,
                sa_lo,
                sa_hi,
                d1,
                0,
                mm_count + 1,
            );
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                q_idx,
                q_idx + 1,
                sc_lo,
                sc_hi,
                d1,
                0,
                mm_count + 1,
            );
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                q_idx,
                q_idx + 1,
                sg_lo,
                sg_hi,
                d1,
                0,
                mm_count + 1,
            );
        }
        BASE_G => {
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                q_idx,
                q_idx + 1,
                sa_lo,
                sa_hi,
                d1,
                0,
                mm_count + 1,
            );
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                q_idx,
                q_idx + 1,
                sg_lo,
                sg_hi,
                d1,
                0,
                mm_count + 1,
            );
            if !WOBBLE {
                recurse_if_nonempty::<F, WOBBLE>(
                    ctx,
                    q_idx,
                    q_idx + 1,
                    su_lo,
                    su_hi,
                    d1,
                    0,
                    mm_count + 1,
                );
            }
        }
        BASE_C => {
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                q_idx,
                q_idx + 1,
                sa_lo,
                sa_hi,
                d1,
                0,
                mm_count + 1,
            );
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                q_idx,
                q_idx + 1,
                sc_lo,
                sc_hi,
                d1,
                0,
                mm_count + 1,
            );
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                q_idx,
                q_idx + 1,
                su_lo,
                su_hi,
                d1,
                0,
                mm_count + 1,
            );
        }
        BASE_U => {
            if !WOBBLE {
                recurse_if_nonempty::<F, WOBBLE>(
                    ctx,
                    q_idx,
                    q_idx + 1,
                    sg_lo,
                    sg_hi,
                    d1,
                    0,
                    mm_count + 1,
                );
            }
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                q_idx,
                q_idx + 1,
                sc_lo,
                sc_hi,
                d1,
                0,
                mm_count + 1,
            );
            recurse_if_nonempty::<F, WOBBLE>(
                ctx,
                q_idx,
                q_idx + 1,
                su_lo,
                su_hi,
                d1,
                0,
                mm_count + 1,
            );
        }
        _ => {}
    }
}

/// S-singleton: one target SA entry, multiple query SA entries (partitioned).
///
/// Uses `is_match_pair` predicate rather than `match s_char` dispatch (as in
/// `recurse_q_singleton`) because iterating query base classes against a known
/// target base is the natural direction when the partition is on the query side.
#[inline(always)]
pub(super) fn recurse_s_singleton<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut SeedingContext<'_, F>,
    ql: usize,
    qr: usize,
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
    partition_interval_into(ctx.q_sa, ctx.q_seq, ql, qr, depth, &mut qint);

    let d1 = depth + 1;
    let can_mm = ctx.max_mm > 0
        && mm_count < ctx.max_mm
        && d1 > ctx.min_prefix
        && ctx.max_len - d1 >= ctx.min_suffix;

    let (qa_lo, qa_hi) = (qint[0], qint[1]);
    let (qc_lo, qc_hi) = (qint[1], qint[2]);
    let (qg_lo, qg_hi) = (qint[2], qint[3]);
    let (qu_lo, qu_hi) = (qint[4], qint[5]);

    if qa_lo < qa_hi && is_match_pair::<WOBBLE>(BASE_A, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(
            ctx,
            qa_lo,
            qa_hi,
            s_idx,
            s_idx + 1,
            d1,
            match_streak + 1,
            mm_count,
        );
    }
    if qc_lo < qc_hi && is_match_pair::<WOBBLE>(BASE_C, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(
            ctx,
            qc_lo,
            qc_hi,
            s_idx,
            s_idx + 1,
            d1,
            match_streak + 1,
            mm_count,
        );
    }
    if qg_lo < qg_hi && is_match_pair::<WOBBLE>(BASE_G, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(
            ctx,
            qg_lo,
            qg_hi,
            s_idx,
            s_idx + 1,
            d1,
            match_streak + 1,
            mm_count,
        );
    }
    if qu_lo < qu_hi && is_match_pair::<WOBBLE>(BASE_U, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(
            ctx,
            qu_lo,
            qu_hi,
            s_idx,
            s_idx + 1,
            d1,
            match_streak + 1,
            mm_count,
        );
    }

    if !can_mm {
        return;
    }

    if qa_lo < qa_hi && !is_match_pair::<WOBBLE>(BASE_A, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qa_lo, qa_hi, s_idx, s_idx + 1, d1, 0, mm_count + 1);
    }
    if qc_lo < qc_hi && !is_match_pair::<WOBBLE>(BASE_C, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qc_lo, qc_hi, s_idx, s_idx + 1, d1, 0, mm_count + 1);
    }
    if qg_lo < qg_hi && !is_match_pair::<WOBBLE>(BASE_G, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qg_lo, qg_hi, s_idx, s_idx + 1, d1, 0, mm_count + 1);
    }
    if qu_lo < qu_hi && !is_match_pair::<WOBBLE>(BASE_U, s_char) {
        recurse_if_nonempty::<F, WOBBLE>(ctx, qu_lo, qu_hi, s_idx, s_idx + 1, d1, 0, mm_count + 1);
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
) {
    let q_suffix_pos = sa_suffix_pos(ctx.q_sa, q_idx);
    let s_suffix_pos = sa_suffix_pos(ctx.t_sa, s_idx);

    loop {
        if depth >= ctx.min_len
            && depth <= ctx.max_len
            && (mm_count == 0 || (match_streak >= ctx.min_suffix && match_streak < depth))
        {
            (ctx.on_match)(SeedMatch {
                query_interval: Interval::new(q_idx, q_idx + 1),
                target_interval: Interval::new(s_idx, s_idx + 1),
                seed_len: depth,
            });
        }

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
