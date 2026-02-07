use std::cmp::max;

use crate::dsm::{stack_with_penalty, DsmModel};

use super::init::{add_e, max3, update_best_with_term};
use super::{ScoreGrid, MIN_SCORE};

#[allow(clippy::too_many_arguments)]
#[cfg_attr(feature = "prof", inline(never))]
pub(super) fn dp_main_loop_generic<const LEFT: bool, M: DsmModel>(
    q_ptr: *const usize,
    t_ptr: *const usize,
    m: &mut ScoreGrid,
    bq: &mut ScoreGrid,
    bt: &mut ScoreGrid,
    width: usize,
    q_len: usize,
    t_len: usize,
    max_stack: i32,
    max_terminal: i32,
    penalty: i32,
    best_e: &mut i32,
    best_i: &mut usize,
    best_j: &mut usize,
) {
    const GAP: usize = 0; // Base::Gap as usize

    // SAFETY invariants:
    // - q_ptr/t_ptr valid for indices [0, q_len) / [0, t_len)
    // - matrices sized at least (q_len+1) x (t_len+1)
    unsafe {
        let m_ptr = m.ptr();
        let bq_ptr = bq.ptr();
        let bt_ptr = bt.ptr();

        // Precompute GAP-GAP profile (constant across all rows)
        // Store lookup_raw(GAP, GAP, t1, t2) at index [t1 * 6 + t2]
        let mut gap_gap_profile = [0i32; 36];
        for t1 in 0..6 {
            for t2 in 0..6 {
                gap_gap_profile[t1 * 6 + t2] = stack_with_penalty::<M>(GAP, GAP, t1, t2, penalty);
            }
        }

        for i in 3..q_len {
            let row_i = i * width;
            let row_prev = (i - 1) * width;
            let qi = *q_ptr.add(i);
            let qi_prev = *q_ptr.add(i - 1);

            // Precompute q-profile for this row (36 lookups, amortized over t_len iterations)
            // Store lookup_raw(q1, q2, t1, t2) at index [t1 * 6 + t2]
            // LEFT:  q1=qi, q2=qi_prev
            // RIGHT: q1=qi_prev, q2=qi
            let mut q_profile = [0i32; 36];
            for t1 in 0..6 {
                for t2 in 0..6 {
                    q_profile[t1 * 6 + t2] = if LEFT {
                        stack_with_penalty::<M>(qi, qi_prev, t1, t2, penalty)
                    } else {
                        stack_with_penalty::<M>(qi_prev, qi, t1, t2, penalty)
                    };
                }
            }

            for j in 3..t_len {
                let diag_idx = row_prev + j - 1;
                let up_idx = row_prev + j;
                let left_idx = row_i + j - 1;
                let curr_idx = row_i + j;
                let tj = *t_ptr.add(j);
                let tj_prev = *t_ptr.add(j - 1);

                // Unified index calculation - compute once, reuse multiple times
                // This eliminates redundant tj * 6 + tj_prev calculations
                let t_stack_idx = if LEFT {
                    tj * 6 + tj_prev
                } else {
                    tj_prev * 6 + tj
                };
                let tj_x6 = tj * 6; // for tj * 6 + GAP cases (GAP = 0)
                                    // Note: GAP * 6 + tj = tj (no computation needed)
                                    // Note: GAP * 6 + GAP = 0 (constant)

                let m_diag = *m_ptr.add(diag_idx);
                let bq_diag = *bq_ptr.add(diag_idx);
                let bt_diag = *bt_ptr.add(diag_idx);

                // Use precomputed t_stack_idx
                let s_mm = add_e(m_diag, q_profile[t_stack_idx]);

                // Use tj_x6 for LEFT, tj for RIGHT
                let s_mq = if LEFT {
                    add_e(bq_diag, q_profile[tj_x6])
                } else {
                    add_e(bq_diag, q_profile[tj])
                };

                let s_mt = if LEFT {
                    add_e(
                        bt_diag,
                        stack_with_penalty::<M>(qi, GAP, tj, tj_prev, penalty),
                    )
                } else {
                    add_e(
                        bt_diag,
                        stack_with_penalty::<M>(GAP, qi, tj_prev, tj, penalty),
                    )
                };
                let val_m = max3(s_mm, s_mq, s_mt);

                if val_m > MIN_SCORE {
                    let remaining = (q_len - i).min(t_len - j) as i32;
                    let upper = val_m + remaining * max_stack + max_terminal;
                    if upper > *best_e {
                        let term = if LEFT {
                            stack_with_penalty::<M>(GAP, qi, GAP, tj, penalty)
                        } else {
                            stack_with_penalty::<M>(qi, GAP, tj, GAP, penalty)
                        };
                        update_best_with_term(best_e, best_i, best_j, val_m, term, i, j);
                    }
                }

                *m_ptr.add(curr_idx) = val_m;

                let m_up = *m_ptr.add(up_idx);
                let bq_up = *bq_ptr.add(up_idx);

                // Use tj for LEFT, tj_x6 for RIGHT
                let s_qm = if LEFT {
                    add_e(m_up, q_profile[tj])
                } else {
                    add_e(m_up, q_profile[tj_x6])
                };
                let s_qq = add_e(bq_up, q_profile[0]); // GAP * 6 + GAP = 0
                *bq_ptr.add(curr_idx) = max(s_qm, s_qq);

                let m_left = *m_ptr.add(left_idx);
                let bt_left = *bt_ptr.add(left_idx);
                let s_tm = if LEFT {
                    add_e(
                        m_left,
                        stack_with_penalty::<M>(GAP, qi, tj, tj_prev, penalty),
                    )
                } else {
                    add_e(
                        m_left,
                        stack_with_penalty::<M>(qi, GAP, tj_prev, tj, penalty),
                    )
                };

                // Use precomputed t_stack_idx
                let s_tt = add_e(bt_left, gap_gap_profile[t_stack_idx]);
                *bt_ptr.add(curr_idx) = max(s_tm, s_tt);
            }
        }
    }
}
