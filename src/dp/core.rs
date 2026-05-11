use std::cmp::max;

use crate::dp::gotoh::Gotoh;

use super::{max3, BestScore, DpCell, DpGrid};

impl Gotoh {
    #[cfg_attr(feature = "prof", inline(never))]
    #[allow(clippy::too_many_arguments)]
    pub(super) fn main_loop(
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

            for i in 3..q_len {
                let prev_row_offset = (i - 1) * width;
                let curr_row_offset = i * width;

                let row = self.row_lookup(*q_ptr.add(i - 1), *q_ptr.add(i));

                // 1. Target Sequence Context
                let mut tp = *t_ptr.add(2) as usize; // Target previous
                let mut tc_ptr = t_ptr.add(3);       // Target current

                // 2. Previous Row Grid Context (Reads only)
                let mut diag = *ptr.add(prev_row_offset + 2);
                let mut up_ptr = ptr.add(prev_row_offset + 3);

                // 3. Current Row Grid Context (Reads and Writes)
                let mut left = *ptr.add(curr_row_offset + 2);
                let mut curr_ptr = ptr.add(curr_row_offset + 3);

                for j in 3..t_len {
                    let up = *up_ptr;
                    let tc = *tc_ptr as usize;

                    let m = max3(
                        diag.m + row.stack(tp, tc),
                        diag.bq + row.close_query_gap(tc),
                        diag.bt + row.close_target_gap(tp, tc),
                    );
                    let bq = max(up.m + row.open_query_gap(tc), up.bq + row.ext_qgap);
                    let bt = max(
                        left.m + row.open_target_gap(tp, tc),
                        left.bt + row.extend_target_gap(tp, tc),
                    );
                    let curr = DpCell { m, bq, bt };

                    best.update_if_better(m, row.terminal(tc), i, j);
                    *curr_ptr = curr;

                    diag = up;
                    left = curr;
                    up_ptr = up_ptr.add(1);
                    curr_ptr = curr_ptr.add(1);
                    tc_ptr = tc_ptr.add(1);
                    tp = tc;
                }
            }
        }
    }
}
