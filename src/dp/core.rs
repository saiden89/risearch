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
            let ptr = grid.as_mut_ptr();
            let width = grid.width();
            let table = self.model.table_ptr();

            for i in 3..q_len {
                let row_i = i * width;
                let row_prev = (i - 1) * width;
                let qi_prev = *q_ptr.add(i - 1) as usize;
                let qi = *q_ptr.add(i) as usize;

                let qp_qc = qi_prev * 216 + qi * 36;
                let stack_ptr = table.add(qp_qc);
                let t_close_ptr = table.add(qi * 36);
                let t_open_ptr = table.add(qi * 216);
                let t_extend_ptr = table;

                let extend_query_gap = *stack_ptr; // GAP = 0

                let mut diag = *ptr.add(row_prev + 2);
                let mut left = *ptr.add(row_i + 2);
                let mut up_ptr = ptr.add(row_prev + 3);
                let mut curr_ptr = ptr.add(row_i + 3);
                let mut t_curr_ptr = t_ptr.add(3);
                let mut tp = *t_ptr.add(2) as usize;

                for j in 3..t_len {
                    let up = *up_ptr;
                    let tc = *t_curr_ptr as usize;
                    let pair_ix = tp * 6 + tc;

                    let m = max3(
                        diag.m + *stack_ptr.add(pair_ix),
                        diag.bq + *stack_ptr.add(tc),
                        diag.bt + *t_close_ptr.add(pair_ix),
                    );
                    let bq = max(up.m + *stack_ptr.add(tc * 6), up.bq + extend_query_gap);
                    let bt = max(
                        left.m + *t_open_ptr.add(pair_ix),
                        left.bt + *t_extend_ptr.add(pair_ix),
                    );
                    let curr = DpCell { m, bq, bt };

                    best.update_if_better(m, *t_open_ptr.add(tc * 6), i, j);
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
