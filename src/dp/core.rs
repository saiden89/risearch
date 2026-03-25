use std::cmp::max;

use crate::dp::gotoh::Gotoh;

use super::{max3, BestScore, DpCell, DpGrid};

/// Threshold where the row-profile-heavy kernel starts to amortize better.
pub(super) const LONG_KERNEL_T_LEN_THRESHOLD: usize = 24;

impl Gotoh {
    #[cfg_attr(feature = "prof", inline(never))]
    #[allow(clippy::too_many_arguments)]
    pub(super) fn dp_main_loop_generic(
        &self,
        q_ptr: *const usize,
        t_ptr: *const usize,
        grid: &mut DpGrid,
        q_len: usize,
        t_len: usize,
        best: &mut BestScore,
    ) {
        if t_len > LONG_KERNEL_T_LEN_THRESHOLD {
            self.dp_main_loop_long(q_ptr, t_ptr, grid, q_len, t_len, best);
        } else {
            self.dp_main_loop_short(q_ptr, t_ptr, grid, q_len, t_len, best);
        }
    }

    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    fn dp_main_loop_short(
        &self,
        q_ptr: *const usize,
        t_ptr: *const usize,
        grid: &mut DpGrid,
        q_len: usize,
        t_len: usize,
        best: &mut BestScore,
    ) {
        // SAFETY: The entire block relies on these invariants established by
        // Gotoh::extend():
        // - q_ptr valid for reads [0, q_len), values ∈ 0..6
        // - t_ptr valid for reads [0, t_len), values ∈ 0..6
        // - grid allocated as (q_len+1) × (t_len+1), so all (i,j) with
        //   i ∈ 0..q_len, j ∈ 0..t_len are in-bounds
        // - All Gotoh lookups receive args ∈ 0..6, satisfying their contracts
        // - t_stack_idx = tp*6+tc ≤ 35 < 36 for all 36-entry profile slices
        unsafe {
            let ptr = grid.ptr();
            let width = grid.width();
            let bt_ext_profile = self.bt_extend_profile();

            for i in 3..q_len {
                let row_i = i * width;
                let row_prev = (i - 1) * width;
                let qi = *q_ptr.add(i);
                let qi_prev = *q_ptr.add(i - 1);

                let q_profile = self.match_profile(qi_prev, qi);
                let bq_ext_e = self.bq_extend(qi_prev, qi);
                let m_from_bt = self.m_from_bt_profile(qi);
                let bt_open = self.bt_open_profile(qi);

                for j in 3..t_len {
                    let diag_idx = row_prev + j - 1;
                    let up_idx = row_prev + j;
                    let left_idx = row_i + j - 1;
                    let curr_idx = row_i + j;
                    let tj = *t_ptr.add(j);
                    let tj_prev = *t_ptr.add(j - 1);

                    let t_stack_idx = tj_prev * 6 + tj;

                    let diag = *ptr.add(diag_idx);
                    let s_mm = diag.m + *q_profile.get_unchecked(t_stack_idx);
                    let s_mq = diag.bq + self.m_from_bq(qi_prev, qi, tj);
                    let s_mt = diag.bt + *m_from_bt.get_unchecked(t_stack_idx);
                    let val_m = max3(s_mm, s_mq, s_mt);

                    best.update(val_m, self.terminal(qi, tj), i, j);

                    let up = *ptr.add(up_idx);
                    let s_qm = up.m + self.bq_open(qi_prev, qi, tj);
                    let s_qq = up.bq + bq_ext_e;

                    let left = *ptr.add(left_idx);
                    let s_tm = left.m + *bt_open.get_unchecked(t_stack_idx);
                    let s_tt = left.bt + *bt_ext_profile.get_unchecked(t_stack_idx);

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
        &self,
        q_ptr: *const usize,
        t_ptr: *const usize,
        grid: &mut DpGrid,
        q_len: usize,
        t_len: usize,
        best: &mut BestScore,
    ) {
        // SAFETY: Same invariants as dp_main_loop_short. Additionally, 6-entry
        // profile slices (m_from_bq_profile, bq_open_profile, terminal_profile)
        // are indexed by tj ∈ 0..6, which is within the 6-entry slice bounds.
        unsafe {
            let ptr = grid.ptr();
            let width = grid.width();
            let bt_ext_profile = self.bt_extend_profile();

            for i in 3..q_len {
                let row_i = i * width;
                let row_prev = (i - 1) * width;
                let qi = *q_ptr.add(i);
                let qi_prev = *q_ptr.add(i - 1);

                let q_profile = self.match_profile(qi_prev, qi);
                let m_from_bt = self.m_from_bt_profile(qi);
                let bt_open = self.bt_open_profile(qi);
                let bq_ext_e = self.bq_extend(qi_prev, qi);
                let m_from_bq_slice = self.m_from_bq_profile(qi_prev, qi);
                let bq_open_slice = self.bq_open_profile(qi_prev, qi);
                let terminal_slice = self.terminal_profile(qi);

                for j in 3..t_len {
                    let diag_idx = row_prev + j - 1;
                    let up_idx = row_prev + j;
                    let left_idx = row_i + j - 1;
                    let curr_idx = row_i + j;
                    let tj = *t_ptr.add(j);
                    let tj_prev = *t_ptr.add(j - 1);

                    let t_stack_idx = tj_prev * 6 + tj;

                    let diag = *ptr.add(diag_idx);
                    let s_mm = diag.m + *q_profile.get_unchecked(t_stack_idx);
                    let s_mq = diag.bq + *m_from_bq_slice.get_unchecked(tj);
                    let s_mt = diag.bt + *m_from_bt.get_unchecked(t_stack_idx);
                    let val_m = max3(s_mm, s_mq, s_mt);
                    best.update(val_m, *terminal_slice.get_unchecked(tj), i, j);

                    let up = *ptr.add(up_idx);
                    let s_qm = up.m + *bq_open_slice.get_unchecked(tj);
                    let s_qq = up.bq + bq_ext_e;

                    let left = *ptr.add(left_idx);
                    let s_tm = left.m + *bt_open.get_unchecked(t_stack_idx);
                    let s_tt = left.bt + *bt_ext_profile.get_unchecked(t_stack_idx);

                    *ptr.add(curr_idx) = DpCell {
                        m: val_m,
                        bq: max(s_qm, s_qq),
                        bt: max(s_tm, s_tt),
                    };
                }
            }
        }
    }
}
