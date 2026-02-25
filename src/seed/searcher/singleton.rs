use crate::types::{Interval, PackedSaEntry};

use super::{partition_interval_into, recurse, RecurseCtx, SeedMatch};

#[inline(always)]
fn is_match_pair<const WOBBLE: bool>(q_char: u8, s_char: u8) -> bool {
    match q_char {
        1 => s_char == 4,                            // A ↔ U
        2 => s_char == 3 || (WOBBLE && s_char == 4), // G ↔ C or wobble U
        3 => s_char == 2,                            // C ↔ G
        4 => s_char == 1 || (WOBBLE && s_char == 2), // U ↔ A or wobble G
        _ => false,
    }
}

#[inline(always)]
pub(super) fn recurse_q_singleton<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut RecurseCtx<'_, F>,
    q_idx: usize,
    sl: usize,
    sr: usize,
    depth: usize,
    msm: usize,
    mc: usize,
) {
    let q_suffix_pos =
        (unsafe { *ctx.q_sa.get_unchecked(q_idx) } & PackedSaEntry::POS_MASK) as usize;
    let q_char = unsafe { *ctx.q_seq.get_unchecked(q_suffix_pos + depth) } as u8;
    if !(1..=4).contains(&q_char) {
        return;
    }

    let mut sint = [0usize; 6];
    partition_interval_into(ctx.t_sa, ctx.t_seq, sl, sr, depth, &mut sint);

    let d1 = depth + 1;
    let can_mm = ctx.max_mm > 0
        && mc < ctx.max_mm
        && d1 > ctx.min_prefix
        && ctx.max_len - d1 >= ctx.min_suffix;

    let (sa_lo, sa_hi) = (sint[0], sint[1]);
    let (sc_lo, sc_hi) = (sint[1], sint[2]);
    let (sg_lo, sg_hi) = (sint[2], sint[3]);
    let (su_lo, su_hi) = (sint[4], sint[5]);

    macro_rules! rec_s {
        ($s_lo:expr, $s_hi:expr, $next_msm:expr, $next_mc:expr) => {
            if $s_lo < $s_hi {
                recurse::<F, WOBBLE>(ctx, q_idx, q_idx + 1, $s_lo, $s_hi, d1, $next_msm, $next_mc);
            }
        };
    }

    match q_char {
        1 => rec_s!(su_lo, su_hi, msm + 1, mc), // A-U
        2 => {
            rec_s!(sc_lo, sc_hi, msm + 1, mc); // G-C
            if WOBBLE {
                rec_s!(su_lo, su_hi, msm + 1, mc); // G-U wobble
            }
        }
        3 => rec_s!(sg_lo, sg_hi, msm + 1, mc), // C-G
        4 => {
            rec_s!(sa_lo, sa_hi, msm + 1, mc); // U-A
            if WOBBLE {
                rec_s!(sg_lo, sg_hi, msm + 1, mc); // U-G wobble
            }
        }
        _ => {}
    }

    if !can_mm {
        return;
    }

    match q_char {
        1 => {
            rec_s!(sa_lo, sa_hi, 0, mc + 1);
            rec_s!(sc_lo, sc_hi, 0, mc + 1);
            rec_s!(sg_lo, sg_hi, 0, mc + 1);
        }
        2 => {
            rec_s!(sa_lo, sa_hi, 0, mc + 1);
            rec_s!(sg_lo, sg_hi, 0, mc + 1);
            if !WOBBLE {
                rec_s!(su_lo, su_hi, 0, mc + 1);
            }
        }
        3 => {
            rec_s!(sa_lo, sa_hi, 0, mc + 1);
            rec_s!(sc_lo, sc_hi, 0, mc + 1);
            rec_s!(su_lo, su_hi, 0, mc + 1);
        }
        4 => {
            if !WOBBLE {
                rec_s!(sg_lo, sg_hi, 0, mc + 1);
            }
            rec_s!(sc_lo, sc_hi, 0, mc + 1);
            rec_s!(su_lo, su_hi, 0, mc + 1);
        }
        _ => {}
    }
}

#[inline(always)]
pub(super) fn recurse_s_singleton<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut RecurseCtx<'_, F>,
    ql: usize,
    qr: usize,
    s_idx: usize,
    depth: usize,
    msm: usize,
    mc: usize,
) {
    let s_suffix_pos =
        (unsafe { *ctx.t_sa.get_unchecked(s_idx) } & PackedSaEntry::POS_MASK) as usize;
    let s_char = unsafe { *ctx.t_seq.get_unchecked(s_suffix_pos + depth) } as u8;
    if !(1..=4).contains(&s_char) {
        return;
    }

    let mut qint = [0usize; 6];
    partition_interval_into(ctx.q_sa, ctx.q_seq, ql, qr, depth, &mut qint);

    let d1 = depth + 1;
    let can_mm = ctx.max_mm > 0
        && mc < ctx.max_mm
        && d1 > ctx.min_prefix
        && ctx.max_len - d1 >= ctx.min_suffix;

    let (qa_lo, qa_hi) = (qint[0], qint[1]);
    let (qc_lo, qc_hi) = (qint[1], qint[2]);
    let (qg_lo, qg_hi) = (qint[2], qint[3]);
    let (qu_lo, qu_hi) = (qint[4], qint[5]);

    macro_rules! rec_q {
        ($q_lo:expr, $q_hi:expr, $next_msm:expr, $next_mc:expr) => {
            if $q_lo < $q_hi {
                recurse::<F, WOBBLE>(ctx, $q_lo, $q_hi, s_idx, s_idx + 1, d1, $next_msm, $next_mc);
            }
        };
    }

    if qa_lo < qa_hi && is_match_pair::<WOBBLE>(1, s_char) {
        rec_q!(qa_lo, qa_hi, msm + 1, mc);
    }
    if qc_lo < qc_hi && is_match_pair::<WOBBLE>(3, s_char) {
        rec_q!(qc_lo, qc_hi, msm + 1, mc);
    }
    if qg_lo < qg_hi && is_match_pair::<WOBBLE>(2, s_char) {
        rec_q!(qg_lo, qg_hi, msm + 1, mc);
    }
    if qu_lo < qu_hi && is_match_pair::<WOBBLE>(4, s_char) {
        rec_q!(qu_lo, qu_hi, msm + 1, mc);
    }

    if !can_mm {
        return;
    }

    if qa_lo < qa_hi && !is_match_pair::<WOBBLE>(1, s_char) {
        rec_q!(qa_lo, qa_hi, 0, mc + 1);
    }
    if qc_lo < qc_hi && !is_match_pair::<WOBBLE>(3, s_char) {
        rec_q!(qc_lo, qc_hi, 0, mc + 1);
    }
    if qg_lo < qg_hi && !is_match_pair::<WOBBLE>(2, s_char) {
        rec_q!(qg_lo, qg_hi, 0, mc + 1);
    }
    if qu_lo < qu_hi && !is_match_pair::<WOBBLE>(4, s_char) {
        rec_q!(qu_lo, qu_hi, 0, mc + 1);
    }
}

#[inline(always)]
pub(super) fn recurse_singleton<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut RecurseCtx<'_, F>,
    q_idx: usize,
    s_idx: usize,
    mut depth: usize,
    mut msm: usize,
    mut mc: usize,
) {
    let q_suffix_pos =
        (unsafe { *ctx.q_sa.get_unchecked(q_idx) } & PackedSaEntry::POS_MASK) as usize;
    let s_suffix_pos =
        (unsafe { *ctx.t_sa.get_unchecked(s_idx) } & PackedSaEntry::POS_MASK) as usize;

    loop {
        if depth >= ctx.min_len
            && depth <= ctx.max_len
            && (mc == 0 || (msm >= ctx.min_suffix && msm < depth))
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

        if mc > 0 && ctx.min_suffix > 0 {
            let max_possible = msm + (ctx.max_len - depth);
            if max_possible < ctx.min_suffix {
                return;
            }
        }

        let d1 = depth + 1;
        let can_mm = ctx.max_mm > 0
            && mc < ctx.max_mm
            && d1 > ctx.min_prefix
            && ctx.max_len - d1 >= ctx.min_suffix;

        let q_char = unsafe { *ctx.q_seq.get_unchecked(q_suffix_pos + depth) } as u8;
        let s_char = unsafe { *ctx.t_seq.get_unchecked(s_suffix_pos + depth) } as u8;
        if !(1..=4).contains(&q_char) || !(1..=4).contains(&s_char) {
            return;
        }

        if is_match_pair::<WOBBLE>(q_char, s_char) {
            depth = d1;
            msm += 1;
            continue;
        }
        if can_mm {
            depth = d1;
            mc += 1;
            msm = 0;
            continue;
        }
        return;
    }
}
