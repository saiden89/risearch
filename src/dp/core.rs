use std::cmp::max;

use crate::dp::gotoh::Gotoh;

use super::{max3, BestScore, DpCell, DpGrid};

impl Gotoh {
    #[cfg_attr(feature = "prof", inline(never))]
    #[allow(clippy::too_many_arguments)]
    pub(super) fn dp_main_loop(
        &self,
        q_ptr: *const u8,
        t_ptr: *const u8,
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
        unsafe {
            let ptr = grid.ptr();
            let width = grid.width();

            for i in 3..q_len {
                let row_i = i * width;
                let row_prev = (i - 1) * width;
                let qi_prev = *q_ptr.add(i - 1);
                let qi = *q_ptr.add(i);

                let stack = self.stack_row(qi_prev, qi).as_ptr();
                let close_query_gap = self.close_query_gap_row(qi_prev, qi).as_ptr();
                let open_query_gap = self.open_query_gap_row(qi_prev, qi).as_ptr();
                let extend_query_gap = self.extend_query_gap(qi_prev, qi);
                let close_target_gap = self.close_target_gap_row(qi).as_ptr();
                let open_target_gap = self.open_target_gap_row(qi).as_ptr();
                let extend_target_gap = self.extend_target_gap_row().as_ptr();
                let terminal = self.terminal_row(qi).as_ptr();

                let mut diag = *ptr.add(row_prev + 2);
                let mut left = *ptr.add(row_i + 2);
                let mut up_ptr = ptr.add(row_prev + 3);
                let mut curr_ptr = ptr.add(row_i + 3);
                let mut t_curr_ptr = t_ptr.add(3);
                let mut tp = *t_ptr.add(2);

                for j in 3..t_len {
                    let up = *up_ptr;
                    let tc = *t_curr_ptr;
                    let tc_ix = usize::from(tc);
                    let pair_ix = usize::from(tp) * 6 + tc_ix;

                    let m = max3(
                        diag.m + *stack.add(pair_ix),
                        diag.bq + *close_query_gap.add(tc_ix),
                        diag.bt + *close_target_gap.add(pair_ix),
                    );
                    let bq = max(up.m + *open_query_gap.add(tc_ix), up.bq + extend_query_gap);
                    let bt = max(
                        left.m + *open_target_gap.add(pair_ix),
                        left.bt + *extend_target_gap.add(pair_ix),
                    );
                    let curr = DpCell { m, bq, bt };

                    best.update(m, *terminal.add(tc_ix), i, j);
                    *curr_ptr = curr;

                    diag = up;
                    left = curr;
                    up_ptr = up_ptr.add(1);
                    curr_ptr = curr_ptr.add(1);
                    t_curr_ptr = t_curr_ptr.add(1);
                    tp = tc;
                }
            }
        }
    }
}
