use crate::dsm::dsm_lookup_raw;

use super::{ScoreOnlyGrid, GAP, MIN_SCORE};
use super::init::{add_e, max2, max3, update_best_with_term};

#[cfg_attr(feature = "prof", inline(never))]
pub(super) fn dp_main_loop_left(
    q_ptr: *const usize,
    t_ptr: *const usize,
    m: &mut ScoreOnlyGrid,
    bq: &mut ScoreOnlyGrid,
    bt: &mut ScoreOnlyGrid,
    width: usize,
    q_len: usize,
    t_len: usize,
    max_stack: i32,
    max_terminal: i32,
    best_e: &mut i32,
    best_i: &mut usize,
    best_j: &mut usize,
) {
    // SAFETY invariants:
    // - q_ptr/t_ptr valid for indices [0, q_len) / [0, t_len)
    // - matrices sized at least (q_len+1) x (t_len+1)
    unsafe {
        let m_ptr = m.ptr();
        let bq_ptr = bq.ptr();
        let bt_ptr = bt.ptr();
        for i in 3..q_len {
            let row_i = i * width;
            let row_prev = (i - 1) * width;
            let qi = *q_ptr.add(i);
            let qi_prev = *q_ptr.add(i - 1);

            for j in 3..t_len {
                let diag_idx = row_prev + j - 1;
                let up_idx = row_prev + j;
                let left_idx = row_i + j - 1;
                let curr_idx = row_i + j;
                let tj = *t_ptr.add(j);
                let tj_prev = *t_ptr.add(j - 1);

                let m_diag = *m_ptr.add(diag_idx);
                let bq_diag = *bq_ptr.add(diag_idx);
                let bt_diag = *bt_ptr.add(diag_idx);

                let s_mm = add_e(m_diag, dsm_lookup_raw(qi, qi_prev, tj, tj_prev));
                let s_mq = add_e(bq_diag, dsm_lookup_raw(qi, qi_prev, tj, GAP));
                let s_mt = add_e(bt_diag, dsm_lookup_raw(qi, GAP, tj, tj_prev));
                let val_m = max3(s_mm, s_mq, s_mt);

                if val_m > MIN_SCORE {
                    let remaining = (q_len - i).min(t_len - j) as i32;
                    let upper = val_m + remaining * max_stack + max_terminal;
                    if upper > *best_e {
                        update_best_with_term(
                            best_e,
                            best_i,
                            best_j,
                            val_m,
                            dsm_lookup_raw(GAP, qi, GAP, tj),
                            i,
                            j,
                        );
                    }
                }

                *m_ptr.add(curr_idx) = val_m;

                let m_up = *m_ptr.add(up_idx);
                let bq_up = *bq_ptr.add(up_idx);
                let s_qm = add_e(m_up, dsm_lookup_raw(qi, qi_prev, GAP, tj));
                let s_qq = add_e(bq_up, dsm_lookup_raw(qi, qi_prev, GAP, GAP));
                *bq_ptr.add(curr_idx) = max2(s_qm, s_qq);

                let m_left = *m_ptr.add(left_idx);
                let bt_left = *bt_ptr.add(left_idx);
                let s_tm = add_e(m_left, dsm_lookup_raw(GAP, qi, tj, tj_prev));
                let s_tt = add_e(bt_left, dsm_lookup_raw(GAP, GAP, tj, tj_prev));
                *bt_ptr.add(curr_idx) = max2(s_tm, s_tt);
            }
        }
    }
}

#[cfg_attr(feature = "prof", inline(never))]
pub(super) fn dp_main_loop_right(
    q_ptr: *const usize,
    t_ptr: *const usize,
    m: &mut ScoreOnlyGrid,
    bq: &mut ScoreOnlyGrid,
    bt: &mut ScoreOnlyGrid,
    width: usize,
    q_len: usize,
    t_len: usize,
    max_stack: i32,
    max_terminal: i32,
    best_e: &mut i32,
    best_i: &mut usize,
    best_j: &mut usize,
) {
    // SAFETY invariants:
    // - q_ptr/t_ptr valid for indices [0, q_len) / [0, t_len)
    // - matrices sized at least (q_len+1) x (t_len+1)
    unsafe {
        let m_ptr = m.ptr();
        let bq_ptr = bq.ptr();
        let bt_ptr = bt.ptr();
        for i in 3..q_len {
            let row_i = i * width;
            let row_prev = (i - 1) * width;
            let qi = *q_ptr.add(i);
            let qi_prev = *q_ptr.add(i - 1);

            for j in 3..t_len {
                let diag_idx = row_prev + j - 1;
                let up_idx = row_prev + j;
                let left_idx = row_i + j - 1;
                let curr_idx = row_i + j;
                let tj = *t_ptr.add(j);
                let tj_prev = *t_ptr.add(j - 1);

                let m_diag = *m_ptr.add(diag_idx);
                let bq_diag = *bq_ptr.add(diag_idx);
                let bt_diag = *bt_ptr.add(diag_idx);

                let s_mm = add_e(m_diag, dsm_lookup_raw(qi_prev, qi, tj_prev, tj));
                let s_mq = add_e(bq_diag, dsm_lookup_raw(qi_prev, qi, GAP, tj));
                let s_mt = add_e(bt_diag, dsm_lookup_raw(GAP, qi, tj_prev, tj));
                let val_m = max3(s_mm, s_mq, s_mt);

                if val_m > MIN_SCORE {
                    let remaining = (q_len - i).min(t_len - j) as i32;
                    let upper = val_m + remaining * max_stack + max_terminal;
                    if upper > *best_e {
                        update_best_with_term(
                            best_e,
                            best_i,
                            best_j,
                            val_m,
                            dsm_lookup_raw(qi, GAP, tj, GAP),
                            i,
                            j,
                        );
                    }
                }

                *m_ptr.add(curr_idx) = val_m;

                let m_up = *m_ptr.add(up_idx);
                let bq_up = *bq_ptr.add(up_idx);
                let s_qm = add_e(m_up, dsm_lookup_raw(qi_prev, qi, tj, GAP));
                let s_qq = add_e(bq_up, dsm_lookup_raw(qi_prev, qi, GAP, GAP));
                *bq_ptr.add(curr_idx) = max2(s_qm, s_qq);

                let m_left = *m_ptr.add(left_idx);
                let bt_left = *bt_ptr.add(left_idx);
                let s_tm = add_e(m_left, dsm_lookup_raw(qi, GAP, tj_prev, tj));
                let s_tt = add_e(bt_left, dsm_lookup_raw(GAP, GAP, tj_prev, tj));
                *bt_ptr.add(curr_idx) = max2(s_tm, s_tt);
            }
        }
    }
}
