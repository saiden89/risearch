use std::cmp::max;

use crate::dsm::{DirectionalDsm, GAP};

use super::{add_e, BestScore, DpCell, DpGrid, NEG_INF};

// =============================================================================
// INIT HELPERS - Reduce code duplication in DP initialization
// =============================================================================

/// Transition helper: choose the best of two predecessor paths.
#[inline(always)]
fn best2(from_m: i32, e_m: i32, from_p: i32, e_p: i32) -> i32 {
    max(add_e(from_m, e_m), add_e(from_p, e_p))
}

/// Get a pointer to a DP cell at `(i, j)` from the grid.
#[inline(always)]
unsafe fn cell(ptr: *mut DpCell, width: usize, i: usize, j: usize) -> *mut DpCell {
    ptr.add(i * width + j)
}

/// Initialize DP boundary cells and bridge into limited row/column setup.
///
/// Returns `true` when the main DP region (`i >= 3`, `j >= 3`) exists.
/// Returns `false` when initialization is complete and no main-loop pass is needed.
#[inline]
#[allow(clippy::too_many_arguments)]
pub(crate) fn init_frontier(
    q_ptr: *const usize,
    t_ptr: *const usize,
    grid: &mut DpGrid,
    q_len: usize,
    t_len: usize,
    dsm: &DirectionalDsm<'_>,
    best: &mut BestScore,
) -> bool {
    let ptr = grid.ptr();
    let width = grid.width();

    debug_assert!(width > t_len, "matrix width too small for t_len");
    unsafe {
        let q0 = *q_ptr.add(0);
        let q1 = *q_ptr.add(1);
        let t0 = *t_ptr.add(0);
        let t1 = *t_ptr.add(1);

        // Corner cells: explicit NA values (like C code)
        *cell(ptr, width, 0, 0) = DpCell {
            m: 0,
            bq: NEG_INF,
            bt: NEG_INF,
        };
        *cell(ptr, width, 0, 1) = DpCell {
            m: NEG_INF,
            bq: NEG_INF,
            bt: dsm.lookup(q0, GAP, t0, t1),
        };
        *cell(ptr, width, 1, 0) = DpCell {
            m: NEG_INF,
            bq: dsm.lookup(q0, q1, t0, GAP),
            bt: NEG_INF,
        };

        let m11 = dsm.lookup(q0, q1, t0, t1);
        *cell(ptr, width, 1, 1) = DpCell {
            m: m11,
            bq: NEG_INF,
            bt: NEG_INF,
        };
        best.update(m11, dsm.terminal(q1, t1), 1, 1);

        // Row 0 (Bt only) and Row 1 (M) - unconditional writes
        for k in 2..t_len {
            let t_prev = *t_ptr.add(k - 1);
            let t_curr = *t_ptr.add(k);
            let prev = (*cell(ptr, width, 0, k - 1)).bt;
            let bt_val = add_e(prev, dsm.lookup(GAP, GAP, t_prev, t_curr));
            let m_val = add_e(prev, dsm.lookup(GAP, q1, t_prev, t_curr));
            *cell(ptr, width, 0, k) = DpCell {
                m: NEG_INF,
                bq: NEG_INF,
                bt: bt_val,
            };
            *cell(ptr, width, 1, k) = DpCell {
                m: m_val,
                bq: NEG_INF,
                bt: NEG_INF,
            };
            best.update(m_val, dsm.terminal(q1, t_curr), 1, k);
        }

        // Col 0 (Bq only) and Col 1 (M) - unconditional writes
        for k in 2..q_len {
            let q_prev = *q_ptr.add(k - 1);
            let q_curr = *q_ptr.add(k);
            let prev = (*cell(ptr, width, k - 1, 0)).bq;
            let bq_val = add_e(prev, dsm.lookup(q_prev, q_curr, GAP, GAP));
            let m_val = add_e(prev, dsm.lookup(q_prev, q_curr, GAP, t1));
            *cell(ptr, width, k, 0) = DpCell {
                m: NEG_INF,
                bq: bq_val,
                bt: NEG_INF,
            };
            *cell(ptr, width, k, 1) = DpCell {
                m: m_val,
                bq: NEG_INF,
                bt: NEG_INF,
            };
            best.update(m_val, dsm.terminal(q_curr, t1), k, 1);
        }
    }

    // Early return when either axis is too small for row/col 2 cells
    if q_len <= 2 || t_len <= 2 {
        return false;
    }

    // Cell (2,2) init: bridge corner to limited rows/cols
    unsafe {
        let q1 = *q_ptr.add(1);
        let q2 = *q_ptr.add(2);
        let t1 = *t_ptr.add(1);
        let t2 = *t_ptr.add(2);
        let m11_val = (*cell(ptr, width, 1, 1)).m;
        let bt12 = add_e(m11_val, dsm.lookup(q1, GAP, t1, t2));
        let bq21 = add_e(m11_val, dsm.lookup(q1, q2, t1, GAP));
        let m22 = add_e(m11_val, dsm.lookup(q1, q2, t1, t2));
        (*cell(ptr, width, 1, 2)).bt = bt12;
        (*cell(ptr, width, 2, 1)).bq = bq21;
        (*cell(ptr, width, 2, 2)).m = m22;
        best.update(m22, dsm.terminal(q2, t2), 2, 2);

        let m12 = (*cell(ptr, width, 1, 2)).m;
        let m21 = (*cell(ptr, width, 2, 1)).m;
        (*cell(ptr, width, 2, 2)).bq = add_e(m12, dsm.lookup(q1, q2, t2, GAP));
        (*cell(ptr, width, 2, 2)).bt = add_e(m21, dsm.lookup(q2, GAP, t1, t2));
    }

    init_limited_rows(q_ptr, t_ptr, grid, q_len, t_len, dsm, best);
    init_limited_cols(q_ptr, t_ptr, grid, q_len, t_len, dsm, best);

    true
}

#[inline]
/// Initialize the limited top rows (`i=1` and `i=2`) for columns `j>=3`.
#[allow(clippy::too_many_arguments)]
pub(super) fn init_limited_rows(
    q_ptr: *const usize,
    t_ptr: *const usize,
    grid: &mut DpGrid,
    q_len: usize,
    t_len: usize,
    dsm: &DirectionalDsm<'_>,
    best: &mut BestScore,
) {
    let ptr = grid.ptr();
    let width = grid.width();

    debug_assert!(
        q_len > 2 && t_len > 2,
        "limited rows require q_len/t_len >= 3"
    );
    unsafe {
        let qi1 = *q_ptr.add(1);
        let qi2 = *q_ptr.add(2);

        let c1_2 = *cell(ptr, width, 1, 2);
        let c2_2 = *cell(ptr, width, 2, 2);
        let mut m1_prev = c1_2.m;
        let mut bt1_prev = c1_2.bt;
        let mut m2_prev = c2_2.m;
        let mut bt2_prev = c2_2.bt;

        for k in 3..t_len {
            let tj = *t_ptr.add(k);
            let tj_prev = *t_ptr.add(k - 1);

            // Primary[1,k] = Bt
            let bt1 = best2(
                m1_prev,
                dsm.lookup(qi1, GAP, tj_prev, tj),
                bt1_prev,
                dsm.lookup(GAP, GAP, tj_prev, tj),
            );
            (*cell(ptr, width, 1, k)).bt = bt1;

            // M[2,k]
            let m2 = best2(
                m1_prev,
                dsm.lookup(qi1, qi2, tj_prev, tj),
                bt1_prev,
                dsm.lookup(GAP, qi2, tj_prev, tj),
            );
            (*cell(ptr, width, 2, k)).m = m2;
            best.update(m2, dsm.terminal(qi2, tj), 2, k);

            // Secondary[2,k] = Bq
            let m1k = (*cell(ptr, width, 1, k)).m;
            (*cell(ptr, width, 2, k)).bq = add_e(m1k, dsm.lookup(qi1, qi2, tj, GAP));

            // Primary[2,k] = Bt
            let bt2 = best2(
                m2_prev,
                dsm.lookup(qi2, GAP, tj_prev, tj),
                bt2_prev,
                dsm.lookup(GAP, GAP, tj_prev, tj),
            );
            (*cell(ptr, width, 2, k)).bt = bt2;

            m1_prev = m1k;
            bt1_prev = bt1;
            m2_prev = m2;
            bt2_prev = bt2;
        }
    }
}

#[inline]
/// Initialize the limited left columns (`j=1` and `j=2`) for rows `i>=3`.
#[allow(clippy::too_many_arguments)]
pub(super) fn init_limited_cols(
    q_ptr: *const usize,
    t_ptr: *const usize,
    grid: &mut DpGrid,
    q_len: usize,
    t_len: usize,
    dsm: &DirectionalDsm<'_>,
    best: &mut BestScore,
) {
    let ptr = grid.ptr();
    let width = grid.width();

    debug_assert!(
        q_len > 2 && t_len > 2,
        "limited cols require q_len/t_len >= 3"
    );
    unsafe {
        let tj1 = *t_ptr.add(1);
        let tj2 = *t_ptr.add(2);

        let c2_1 = *cell(ptr, width, 2, 1);
        let c2_2 = *cell(ptr, width, 2, 2);
        let mut m1_prev = c2_1.m;
        let mut bq1_prev = c2_1.bq;
        let mut m2_prev = c2_2.m;
        let mut bq2_prev = c2_2.bq;

        for k in 3..q_len {
            let qi = *q_ptr.add(k);
            let qi_prev = *q_ptr.add(k - 1);

            // Primary[ k,1 ] = Bq
            let bq1 = best2(
                m1_prev,
                dsm.lookup(qi_prev, qi, tj1, GAP),
                bq1_prev,
                dsm.lookup(qi_prev, qi, GAP, GAP),
            );
            (*cell(ptr, width, k, 1)).bq = bq1;

            // M[ k,2 ]
            let val = best2(
                m1_prev,
                dsm.lookup(qi_prev, qi, tj1, tj2),
                bq1_prev,
                dsm.lookup(qi_prev, qi, GAP, tj2),
            );
            (*cell(ptr, width, k, 2)).m = val;
            best.update(val, dsm.terminal(qi, tj2), k, 2);

            // Secondary[ k,2 ] = Bt
            let m_k1 = (*cell(ptr, width, k, 1)).m;
            (*cell(ptr, width, k, 2)).bt = add_e(m_k1, dsm.lookup(qi, GAP, tj1, tj2));

            // Primary[ k,2 ] = Bq
            let bq2 = best2(
                m2_prev,
                dsm.lookup(qi_prev, qi, tj2, GAP),
                bq2_prev,
                dsm.lookup(qi_prev, qi, GAP, GAP),
            );
            (*cell(ptr, width, k, 2)).bq = bq2;

            m1_prev = m_k1;
            bq1_prev = bq1;
            m2_prev = val;
            bq2_prev = bq2;
        }
    }
}
