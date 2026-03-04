use std::cmp::max;

use crate::dsm::{DsmModel, GAP};

use super::{add_e, BestScore, DpCell, DpGrid, DpView, ExtendDir, MIN_SCORE};

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

/// Direction-aware stack lookup.
#[inline(always)]
fn stack_dir<const LEFT: bool>(
    q1: usize,
    q2: usize,
    t1: usize,
    t2: usize,
    model: &DsmModel,
) -> i32 {
    if LEFT {
        model.lookup(q2, q1, t2, t1)
    } else {
        model.lookup(q1, q2, t1, t2)
    }
}

/// Initialize DP boundary cells and bridge into limited row/column setup.
///
/// Returns `true` when the main DP region (`i >= 3`, `j >= 3`) exists.
/// Returns `false` when initialization is complete and no main-loop pass is needed.
#[inline]
#[allow(clippy::too_many_arguments)]
pub(crate) fn init_frontier(
    view: &DpView<'_>,
    q_ptr: *const usize,
    t_ptr: *const usize,
    grid: &mut DpGrid,
    q_len: usize,
    t_len: usize,
    model: &DsmModel,
    best: &mut BestScore,
) -> bool {
    let ptr = grid.ptr();
    let width = grid.width();

    debug_assert!(width > t_len, "matrix width too small for t_len");
    unsafe {
        // Corner cells: explicit NA values (like C code)
        *cell(ptr, width, 0, 0) = DpCell {
            m: 0,
            bq: MIN_SCORE,
            bt: MIN_SCORE,
        };
        *cell(ptr, width, 0, 1) = DpCell {
            m: MIN_SCORE,
            bq: MIN_SCORE,
            bt: view.bt_open(0, 1, model),
        };
        *cell(ptr, width, 1, 0) = DpCell {
            m: MIN_SCORE,
            bq: view.bq_open(1, 0, model),
            bt: MIN_SCORE,
        };

        let m11 = view.match_e(1, 1, model);
        *cell(ptr, width, 1, 1) = DpCell {
            m: m11,
            bq: MIN_SCORE,
            bt: MIN_SCORE,
        };
        best.update(m11, view.terminal(1, 1, model), 1, 1);

        // Row 0 (Bt only) and Row 1 (M) - unconditional writes
        for k in 2..t_len {
            let prev = (*cell(ptr, width, 0, k - 1)).bt;
            // Always write - use add_e to propagate MIN_SCORE
            let bt_val = add_e(prev, view.bt_ext(k, model));
            let m_val = add_e(prev, view.m_from_bt(1, k, model));
            *cell(ptr, width, 0, k) = DpCell {
                m: MIN_SCORE,
                bq: MIN_SCORE,
                bt: bt_val,
            };
            *cell(ptr, width, 1, k) = DpCell {
                m: m_val,
                bq: MIN_SCORE,
                bt: MIN_SCORE,
            };
            best.update(m_val, view.terminal(1, k, model), 1, k);
        }

        // Col 0 (Bq only) and Col 1 (M) - unconditional writes
        for k in 2..q_len {
            let prev = (*cell(ptr, width, k - 1, 0)).bq;
            let bq_val = add_e(prev, view.bq_ext(k, model));
            let m_val = add_e(prev, view.m_from_bq(k, 1, model));
            *cell(ptr, width, k, 0) = DpCell {
                m: MIN_SCORE,
                bq: bq_val,
                bt: MIN_SCORE,
            };
            *cell(ptr, width, k, 1) = DpCell {
                m: m_val,
                bq: MIN_SCORE,
                bt: MIN_SCORE,
            };
            best.update(m_val, view.terminal(k, 1, model), k, 1);
        }
    }

    // Early return when either axis is too small for row/col 2 cells
    if q_len <= 2 || t_len <= 2 {
        return false;
    }

    // Cell (2,2) init: bridge corner to limited rows/cols
    unsafe {
        let m11_val = (*cell(ptr, width, 1, 1)).m;
        let bt12 = add_e(m11_val, view.bt_open(1, 2, model));
        let bq21 = add_e(m11_val, view.bq_open(2, 1, model));
        let m22 = add_e(m11_val, view.match_e(2, 2, model));
        (*cell(ptr, width, 1, 2)).bt = bt12;
        (*cell(ptr, width, 2, 1)).bq = bq21;
        (*cell(ptr, width, 2, 2)).m = m22;
        best.update(m22, view.terminal(2, 2, model), 2, 2);

        let m12 = (*cell(ptr, width, 1, 2)).m;
        let m21 = (*cell(ptr, width, 2, 1)).m;
        (*cell(ptr, width, 2, 2)).bq = add_e(m12, view.bq_open(2, 2, model));
        (*cell(ptr, width, 2, 2)).bt = add_e(m21, view.bt_open(2, 2, model));
    }

    if view.dir == ExtendDir::Left {
        init_limited_rows::<true>(view, q_ptr, t_ptr, grid, q_len, t_len, model, best);
        init_limited_cols::<true>(view, q_ptr, t_ptr, grid, q_len, t_len, model, best);
    } else {
        init_limited_rows::<false>(view, q_ptr, t_ptr, grid, q_len, t_len, model, best);
        init_limited_cols::<false>(view, q_ptr, t_ptr, grid, q_len, t_len, model, best);
    }

    true
}

#[inline]
/// Initialize the limited top rows (`i=1` and `i=2`) for columns `j>=3`.
///
/// This computes the seed-adjacent band that bridges boundary initialization
/// to the main DP region, while updating global best score candidates.
#[allow(clippy::too_many_arguments)]
pub(super) fn init_limited_rows<const LEFT: bool>(
    view: &DpView<'_>,
    q_ptr: *const usize,
    t_ptr: *const usize,
    grid: &mut DpGrid,
    q_len: usize,
    t_len: usize,
    model: &DsmModel,
    best: &mut BestScore,
) {
    let ptr = grid.ptr();
    let width = grid.width();

    // SAFETY invariants:
    // - q_ptr and t_ptr are valid for indices [0, q_len) and [0, t_len)
    // - q_len and t_len are >= 3 (we index 1 and 2)
    // - ptr points to a grid sized at least (q_len+1) x (t_len+1)
    debug_assert!(
        q_len > 2 && t_len > 2,
        "limited rows require q_len/t_len >= 3"
    );
    unsafe {
        let qi1 = *q_ptr.add(1);
        let qi2 = *q_ptr.add(2);

        // Rolling values for row 1 and 2
        let c1_2 = *cell(ptr, width, 1, 2);
        let c2_2 = *cell(ptr, width, 2, 2);
        let mut m1_prev = c1_2.m; // m(1, k-1)
        let mut bt1_prev = c1_2.bt; // bt(1, k-1)
        let mut m2_prev = c2_2.m; // m(2, k-1)
        let mut bt2_prev = c2_2.bt; // bt(2, k-1)

        for k in 3..t_len {
            let tj = *t_ptr.add(k);
            let tj_prev = *t_ptr.add(k - 1);

            // Primary[1,k] = Bt
            let bt1 = best2(
                m1_prev,
                stack_dir::<LEFT>(qi1, GAP, tj_prev, tj, model),
                bt1_prev,
                stack_dir::<LEFT>(GAP, GAP, tj_prev, tj, model),
            );
            (*cell(ptr, width, 1, k)).bt = bt1;

            // M[2,k]
            let m2 = best2(
                m1_prev,
                stack_dir::<LEFT>(qi1, qi2, tj_prev, tj, model),
                bt1_prev,
                stack_dir::<LEFT>(GAP, qi2, tj_prev, tj, model),
            );
            (*cell(ptr, width, 2, k)).m = m2;
            best.update(m2, view.terminal(2, k, model), 2, k);

            // Secondary[2,k] = Bq
            let m1k = (*cell(ptr, width, 1, k)).m;
            (*cell(ptr, width, 2, k)).bq = add_e(m1k, stack_dir::<LEFT>(qi1, qi2, tj, GAP, model));

            // Primary[2,k] = Bt
            let bt2 = best2(
                m2_prev,
                stack_dir::<LEFT>(qi2, GAP, tj_prev, tj, model),
                bt2_prev,
                stack_dir::<LEFT>(GAP, GAP, tj_prev, tj, model),
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
///
/// Symmetric companion of `init_limited_rows`, using rolling predecessor
/// values to keep pointer traffic low and logic easy to audit.
#[allow(clippy::too_many_arguments)]
pub(super) fn init_limited_cols<const LEFT: bool>(
    view: &DpView<'_>,
    q_ptr: *const usize,
    t_ptr: *const usize,
    grid: &mut DpGrid,
    q_len: usize,
    t_len: usize,
    model: &DsmModel,
    best: &mut BestScore,
) {
    let ptr = grid.ptr();
    let width = grid.width();

    // SAFETY invariants:
    // - q_ptr and t_ptr are valid for indices [0, q_len) and [0, t_len)
    // - q_len and t_len are >= 3 (we index 1 and 2)
    // - ptr points to a grid sized at least (q_len+1) x (t_len+1)
    debug_assert!(
        q_len > 2 && t_len > 2,
        "limited cols require q_len/t_len >= 3"
    );
    unsafe {
        let tj1 = *t_ptr.add(1);
        let tj2 = *t_ptr.add(2);

        // Rolling values for rows (k-1) at columns 1 and 2.
        let c2_1 = *cell(ptr, width, 2, 1);
        let c2_2 = *cell(ptr, width, 2, 2);
        let mut m1_prev = c2_1.m; // m(k-1,1)
        let mut bq1_prev = c2_1.bq; // bq(k-1,1)
        let mut m2_prev = c2_2.m; // m(k-1,2)
        let mut bq2_prev = c2_2.bq; // bq(k-1,2)

        for k in 3..q_len {
            let qi = *q_ptr.add(k);
            let qi_prev = *q_ptr.add(k - 1);

            // Primary[ k,1 ] = Bq
            let bq1 = best2(
                m1_prev,
                stack_dir::<LEFT>(qi_prev, qi, tj1, GAP, model),
                bq1_prev,
                stack_dir::<LEFT>(qi_prev, qi, GAP, GAP, model),
            );
            (*cell(ptr, width, k, 1)).bq = bq1;

            // M[ k,2 ]
            let val = best2(
                m1_prev,
                stack_dir::<LEFT>(qi_prev, qi, tj1, tj2, model),
                bq1_prev,
                stack_dir::<LEFT>(qi_prev, qi, GAP, tj2, model),
            );
            (*cell(ptr, width, k, 2)).m = val;
            best.update(val, view.terminal(k, 2, model), k, 2);

            // Secondary[ k,2 ] = Bt
            let m_k1 = (*cell(ptr, width, k, 1)).m; // Preinitialized in frontier.
            (*cell(ptr, width, k, 2)).bt = add_e(m_k1, stack_dir::<LEFT>(qi, GAP, tj1, tj2, model));

            // Primary[ k,2 ] = Bq
            let bq2 = best2(
                m2_prev,
                stack_dir::<LEFT>(qi_prev, qi, tj2, GAP, model),
                bq2_prev,
                stack_dir::<LEFT>(qi_prev, qi, GAP, GAP, model),
            );
            (*cell(ptr, width, k, 2)).bq = bq2;

            // Advance rolling state.
            m1_prev = m_k1;
            bq1_prev = bq1;
            m2_prev = val;
            bq2_prev = bq2;
        }
    }
}
