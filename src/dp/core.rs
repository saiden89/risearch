use std::cmp::max;

use crate::dsm::{DirectionalDsm, GAP};

use super::{add_e, max3, BestScore, DpCell, DpGrid};

/// Threshold where the row-profile-heavy kernel starts to amortize better.
pub(super) const LONG_KERNEL_T_LEN_THRESHOLD: usize = 24;

#[cfg_attr(feature = "prof", inline(never))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn dp_main_loop_generic(
    q_ptr: *const usize,
    t_ptr: *const usize,
    grid: &mut DpGrid,
    q_len: usize,
    t_len: usize,
    dsm: &DirectionalDsm<'_>,
    best: &mut BestScore,
) {
    if t_len > LONG_KERNEL_T_LEN_THRESHOLD {
        dp_main_loop_long(q_ptr, t_ptr, grid, q_len, t_len, dsm, best);
    } else {
        dp_main_loop_short(q_ptr, t_ptr, grid, q_len, t_len, dsm, best);
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn dp_main_loop_short(
    q_ptr: *const usize,
    t_ptr: *const usize,
    grid: &mut DpGrid,
    q_len: usize,
    t_len: usize,
    dsm: &DirectionalDsm<'_>,
    best: &mut BestScore,
) {
    // SAFETY invariants:
    // - q_ptr/t_ptr valid for indices [0, q_len) / [0, t_len)
    // - grid sized at least (q_len+1) x (t_len+1)
    unsafe {
        let ptr = grid.ptr();
        let width = grid.width();
        let gap_gap_profile = dsm.gap_gap_profile();

        for i in 3..q_len {
            let row_i = i * width;
            let row_prev = (i - 1) * width;
            let qi = *q_ptr.add(i);
            let qi_prev = *q_ptr.add(i - 1);

            let mut q_profile = [0i32; 36];
            for t1 in 0..6 {
                for t2 in 0..6 {
                    q_profile[t1 * 6 + t2] = dsm.lookup(qi_prev, qi, t1, t2);
                }
            }

            for j in 3..t_len {
                let diag_idx = row_prev + j - 1;
                let up_idx = row_prev + j;
                let left_idx = row_i + j - 1;
                let curr_idx = row_i + j;
                let tj = *t_ptr.add(j);
                let tj_prev = *t_ptr.add(j - 1);

                let t_stack_idx = tj_prev * 6 + tj;
                let tj_x6 = tj * 6;

                let diag = *ptr.add(diag_idx);
                let s_mm = add_e(diag.m, q_profile[t_stack_idx]);
                let s_mq = add_e(diag.bq, q_profile[tj]);
                let s_mt = add_e(diag.bt, dsm.lookup(GAP, qi, tj_prev, tj));
                let val_m = max3(s_mm, s_mq, s_mt);

                let term = dsm.terminal(qi, tj);
                best.update(val_m, term, i, j);

                let up = *ptr.add(up_idx);
                let s_qm = add_e(up.m, q_profile[tj_x6]);
                let s_qq = add_e(up.bq, q_profile[0]); // GAP * 6 + GAP = 0

                let left = *ptr.add(left_idx);
                let s_tm = add_e(left.m, dsm.lookup(qi, GAP, tj_prev, tj));
                let s_tt = add_e(left.bt, gap_gap_profile[t_stack_idx]);

                *ptr.add(curr_idx) = DpCell {
                    m: val_m,
                    bq: max(s_qm, s_qq),
                    bt: max(s_tm, s_tt),
                };
            }
        }
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn dp_main_loop_long(
    q_ptr: *const usize,
    t_ptr: *const usize,
    grid: &mut DpGrid,
    q_len: usize,
    t_len: usize,
    dsm: &DirectionalDsm<'_>,
    best: &mut BestScore,
) {
    // SAFETY invariants:
    // - q_ptr/t_ptr valid for indices [0, q_len) / [0, t_len)
    // - grid sized at least (q_len+1) x (t_len+1)
    unsafe {
        let ptr = grid.ptr();
        let width = grid.width();
        let gap_gap_profile = dsm.gap_gap_profile();

        for i in 3..q_len {
            let row_i = i * width;
            let row_prev = (i - 1) * width;
            let qi = *q_ptr.add(i);
            let qi_prev = *q_ptr.add(i - 1);

            let mut q_profile = [0i32; 36];
            let mut bt_to_m_profile = [0i32; 36];
            let mut m_to_bt_profile = [0i32; 36];
            for t1 in 0..6 {
                for t2 in 0..6 {
                    let idx = t1 * 6 + t2;
                    q_profile[idx] = dsm.lookup(qi_prev, qi, t1, t2);
                    bt_to_m_profile[idx] = dsm.lookup(GAP, qi, t1, t2);
                    m_to_bt_profile[idx] = dsm.lookup(qi, GAP, t1, t2);
                }
            }

            let mut terminal_profile = [0i32; 6];
            for (t, item) in terminal_profile.iter_mut().enumerate() {
                *item = dsm.terminal(qi, t);
            }

            for j in 3..t_len {
                let diag_idx = row_prev + j - 1;
                let up_idx = row_prev + j;
                let left_idx = row_i + j - 1;
                let curr_idx = row_i + j;
                let tj = *t_ptr.add(j);
                let tj_prev = *t_ptr.add(j - 1);

                let t_stack_idx = tj_prev * 6 + tj;
                let tj_x6 = tj * 6;

                let diag = *ptr.add(diag_idx);
                let s_mm = add_e(diag.m, q_profile[t_stack_idx]);
                let s_mq = add_e(diag.bq, q_profile[tj]);
                let s_mt = add_e(diag.bt, bt_to_m_profile[t_stack_idx]);
                let val_m = max3(s_mm, s_mq, s_mt);
                best.update(val_m, terminal_profile[tj], i, j);

                let up = *ptr.add(up_idx);
                let s_qm = add_e(up.m, q_profile[tj_x6]);
                let s_qq = add_e(up.bq, q_profile[0]);

                let left = *ptr.add(left_idx);
                let s_tm = add_e(left.m, m_to_bt_profile[t_stack_idx]);
                let s_tt = add_e(left.bt, gap_gap_profile[t_stack_idx]);

                *ptr.add(curr_idx) = DpCell {
                    m: val_m,
                    bq: max(s_qm, s_qq),
                    bt: max(s_tm, s_tt),
                };
            }
        }
    }
}
