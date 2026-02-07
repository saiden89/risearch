use std::cmp::max;

use crate::dsm::{stack_with_penalty, DsmModel};

use super::{DpView, ExtendDir, GAP, MAX_EXT, MIN_SCORE};

// =============================================================================
// INIT HELPERS - Reduce code duplication in DP initialization
// =============================================================================

#[inline(always)]
pub(super) fn update_best_with_term(
    best_e: &mut i32,
    best_i: &mut usize,
    best_j: &mut usize,
    val: i32,
    term: i32,
    i: usize,
    j: usize,
) {
    if val > MIN_SCORE {
        let curr = val + term;
        if curr > *best_e {
            *best_e = curr;
            *best_i = i;
            *best_j = j;
        }
    }
}

/// Helper: add energy if base is valid (not MIN_SCORE).
/// Simple branch - LLVM optimizes to CMOV when beneficial.
#[inline(always)]
pub(super) fn add_e(base: i32, energy: i32) -> i32 {
    if base > MIN_SCORE {
        base + energy
    } else {
        MIN_SCORE
    }
}

/// max of 3 values - branchless
#[inline(always)]
pub(super) fn max3(a: i32, b: i32, c: i32) -> i32 {
    max(max(a, b), c)
}

/// Transition helper: choose the best of two predecessor paths.
#[inline(always)]
fn best2(from_m: i32, e_m: i32, from_p: i32, e_p: i32) -> i32 {
    max(add_e(from_m, e_m), add_e(from_p, e_p))
}

/// Read one matrix cell at `(i, j)` from a raw DP buffer.
#[inline(always)]
unsafe fn cell_get(ptr: *mut i32, width: usize, i: usize, j: usize) -> i32 {
    unsafe { *ptr.add(i * width + j) }
}

/// Write one matrix cell at `(i, j)` into a raw DP buffer.
#[inline(always)]
unsafe fn cell_set(ptr: *mut i32, width: usize, i: usize, j: usize, val: i32) {
    unsafe { *ptr.add(i * width + j) = val };
}

/// Direction-aware stack lookup.
#[inline(always)]
fn stack_dir<const LEFT: bool, M: DsmModel>(
    q1: usize,
    q2: usize,
    t1: usize,
    t2: usize,
    penalty: i32,
) -> i32 {
    if LEFT {
        stack_with_penalty::<M>(q2, q1, t2, t1, penalty)
    } else {
        stack_with_penalty::<M>(q1, q2, t1, t2, penalty)
    }
}

/// Initialize DP boundary cells and bridge into limited row/column setup.
///
/// Returns `true` when the main DP region (`i >= 3`, `j >= 3`) exists.
/// Returns `false` when initialization is complete and no main-loop pass is needed.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(super) fn init_frontier<M: DsmModel>(
    view: &DpView<'_, M>,
    q_ptr: *const usize,
    t_ptr: *const usize,
    m_ptr: *mut i32,
    bq_ptr: *mut i32,
    bt_ptr: *mut i32,
    width: usize,
    q_len: usize,
    t_len: usize,
    best_e: &mut i32,
    best_i: &mut usize,
    best_j: &mut usize,
) -> bool {
    debug_assert!(width > t_len, "matrix width too small for t_len");
    unsafe {
        // Corner cells: explicit NA values (like C code)
        cell_set(m_ptr, width, 0, 0, 0);
        cell_set(bq_ptr, width, 0, 0, MIN_SCORE);
        cell_set(bt_ptr, width, 0, 0, MIN_SCORE);
        cell_set(m_ptr, width, 0, 1, MIN_SCORE);
        cell_set(bq_ptr, width, 0, 1, MIN_SCORE);
        cell_set(m_ptr, width, 1, 0, MIN_SCORE);
        cell_set(bt_ptr, width, 1, 0, MIN_SCORE);
        cell_set(bq_ptr, width, 1, 1, MIN_SCORE);
        cell_set(bt_ptr, width, 1, 1, MIN_SCORE);

        // Valid corner values
        cell_set(bt_ptr, width, 0, 1, view.bt_open(0, 1));
        cell_set(bq_ptr, width, 1, 0, view.bq_open(1, 0));
        let m11 = view.match_e(1, 1);
        cell_set(m_ptr, width, 1, 1, m11);
        update_best_with_term(best_e, best_i, best_j, m11, view.terminal(1, 1), 1, 1);

        // Row 0 (Bt only) and Row 1 (M) - unconditional writes
        for k in 2..t_len {
            let prev = cell_get(bt_ptr, width, 0, k - 1);
            // Always write - use add_e to propagate MIN_SCORE
            let bt_val = add_e(prev, view.bt_ext(k));
            let m_val = add_e(prev, view.m_from_bt(1, k));
            cell_set(bt_ptr, width, 0, k, bt_val);
            cell_set(m_ptr, width, 0, k, MIN_SCORE); // M[0,k] is NA
            cell_set(bq_ptr, width, 0, k, MIN_SCORE); // Bq[0,k] is NA
            cell_set(bq_ptr, width, 1, k, MIN_SCORE); // Bq[1,k] is NA
            cell_set(m_ptr, width, 1, k, m_val);
            update_best_with_term(best_e, best_i, best_j, m_val, view.terminal(1, k), 1, k);
        }

        // Col 0 (Bq only) and Col 1 (M) - unconditional writes
        for k in 2..q_len {
            let prev = cell_get(bq_ptr, width, k - 1, 0);
            let bq_val = add_e(prev, view.bq_ext(k));
            let m_val = add_e(prev, view.m_from_bq(k, 1));
            cell_set(bq_ptr, width, k, 0, bq_val);
            cell_set(m_ptr, width, k, 0, MIN_SCORE); // M[k,0] is NA
            cell_set(bt_ptr, width, k, 0, MIN_SCORE); // Bt[k,0] is NA
            cell_set(bt_ptr, width, k, 1, MIN_SCORE); // Bt[k,1] is NA
            cell_set(m_ptr, width, k, 1, m_val);
            update_best_with_term(best_e, best_i, best_j, m_val, view.terminal(k, 1), k, 1);
        }
    }

    // Early return when either axis is too small for row/col 2 cells
    if q_len <= 2 || t_len <= 2 {
        return false;
    }

    // Cell (2,2) init: bridge corner to limited rows/cols
    unsafe {
        let m11_val = cell_get(m_ptr, width, 1, 1);
        let bt12 = add_e(m11_val, view.bt_open(1, 2));
        let bq21 = add_e(m11_val, view.bq_open(2, 1));
        let m22 = add_e(m11_val, view.match_e(2, 2));
        cell_set(bt_ptr, width, 1, 2, bt12);
        cell_set(bq_ptr, width, 2, 1, bq21);
        cell_set(m_ptr, width, 2, 2, m22);
        update_best_with_term(best_e, best_i, best_j, m22, view.terminal(2, 2), 2, 2);

        let m12 = cell_get(m_ptr, width, 1, 2);
        let m21 = cell_get(m_ptr, width, 2, 1);
        cell_set(bq_ptr, width, 2, 2, add_e(m12, view.bq_open(2, 2)));
        cell_set(bt_ptr, width, 2, 2, add_e(m21, view.bt_open(2, 2)));
    }

    if view.dir == ExtendDir::Left {
        init_limited_rows::<true, M>(
            view, q_ptr, t_ptr, m_ptr, bt_ptr, bq_ptr, width, q_len, t_len, best_e, best_i, best_j,
        );
        init_limited_cols::<true, M>(
            view, q_ptr, t_ptr, m_ptr, bq_ptr, bt_ptr, width, q_len, t_len, best_e, best_i, best_j,
        );
    } else {
        init_limited_rows::<false, M>(
            view, q_ptr, t_ptr, m_ptr, bt_ptr, bq_ptr, width, q_len, t_len, best_e, best_i, best_j,
        );
        init_limited_cols::<false, M>(
            view, q_ptr, t_ptr, m_ptr, bq_ptr, bt_ptr, width, q_len, t_len, best_e, best_i, best_j,
        );
    }

    true
}

#[allow(clippy::too_many_arguments)]
#[inline(always)]
/// Initialize the limited top rows (`i=1` and `i=2`) for columns `j>=3`.
///
/// This computes the seed-adjacent band that bridges boundary initialization
/// to the main DP region, while updating global best score candidates.
pub(super) fn init_limited_rows<const LEFT: bool, M: DsmModel>(
    view: &DpView<'_, M>,
    q_ptr: *const usize,
    t_ptr: *const usize,
    m_ptr: *mut i32,
    bt_ptr: *mut i32,
    bq_ptr: *mut i32,
    width: usize,
    q_len: usize,
    t_len: usize,
    best_e: &mut i32,
    best_i: &mut usize,
    best_j: &mut usize,
) {
    // SAFETY invariants:
    // - q_ptr and t_ptr are valid for indices [0, q_len) and [0, t_len)
    // - q_len and t_len are >= 3 (we index 1 and 2)
    // - m_ptr/bt_ptr/bq_ptr point to matrices sized at least (q_len+1) x (t_len+1)
    debug_assert!(
        q_len > 2 && t_len > 2,
        "limited rows require q_len/t_len >= 3"
    );
    debug_assert!(q_len <= MAX_EXT, "q_len exceeds precomputed index capacity");
    debug_assert!(t_len <= MAX_EXT, "t_len exceeds precomputed index capacity");
    unsafe {
        let qi1 = *q_ptr.add(1);
        let qi2 = *q_ptr.add(2);

        // Rolling values for row 1 and 2
        let mut m1_prev = cell_get(m_ptr, width, 1, 2); // m(1, k-1)
        let mut bt1_prev = cell_get(bt_ptr, width, 1, 2); // bt(1, k-1)
        let mut m2_prev = cell_get(m_ptr, width, 2, 2); // m(2, k-1)
        let mut bt2_prev = cell_get(bt_ptr, width, 2, 2); // bt(2, k-1)

        for k in 3..t_len {
            let tj = *t_ptr.add(k);
            let tj_prev = *t_ptr.add(k - 1);

            // Primary[1,k] = Bt
            let bt1 = best2(
                m1_prev,
                stack_dir::<LEFT, M>(qi1, GAP, tj_prev, tj, view.penalty),
                bt1_prev,
                stack_dir::<LEFT, M>(GAP, GAP, tj_prev, tj, view.penalty),
            );
            cell_set(bt_ptr, width, 1, k, bt1);

            // M[2,k]
            let m2 = best2(
                m1_prev,
                stack_dir::<LEFT, M>(qi1, qi2, tj_prev, tj, view.penalty),
                bt1_prev,
                stack_dir::<LEFT, M>(GAP, qi2, tj_prev, tj, view.penalty),
            );
            cell_set(m_ptr, width, 2, k, m2);
            update_best_with_term(best_e, best_i, best_j, m2, view.terminal(2, k), 2, k);

            // Secondary[2,k] = Bq
            let m1k = cell_get(m_ptr, width, 1, k);
            cell_set(
                bq_ptr,
                width,
                2,
                k,
                add_e(m1k, stack_dir::<LEFT, M>(qi1, qi2, tj, GAP, view.penalty)),
            );

            // Primary[2,k] = Bt
            let bt2 = best2(
                m2_prev,
                stack_dir::<LEFT, M>(qi2, GAP, tj_prev, tj, view.penalty),
                bt2_prev,
                stack_dir::<LEFT, M>(GAP, GAP, tj_prev, tj, view.penalty),
            );
            cell_set(bt_ptr, width, 2, k, bt2);

            m1_prev = cell_get(m_ptr, width, 1, k);
            bt1_prev = bt1;
            m2_prev = m2;
            bt2_prev = bt2;
        }
    }
}

#[allow(clippy::too_many_arguments)]
#[inline(always)]
/// Initialize the limited left columns (`j=1` and `j=2`) for rows `i>=3`.
///
/// Symmetric companion of `init_limited_rows`, using rolling predecessor
/// values to keep pointer traffic low and logic easy to audit.
pub(super) fn init_limited_cols<const LEFT: bool, M: DsmModel>(
    view: &DpView<'_, M>,
    q_ptr: *const usize,
    t_ptr: *const usize,
    m_ptr: *mut i32,
    bq_ptr: *mut i32,
    bt_ptr: *mut i32,
    width: usize,
    q_len: usize,
    t_len: usize,
    best_e: &mut i32,
    best_i: &mut usize,
    best_j: &mut usize,
) {
    // SAFETY invariants:
    // - q_ptr and t_ptr are valid for indices [0, q_len) and [0, t_len)
    // - q_len and t_len are >= 3 (we index 1 and 2)
    // - m_ptr/bq_ptr/bt_ptr point to matrices sized at least (q_len+1) x (t_len+1)
    debug_assert!(
        q_len > 2 && t_len > 2,
        "limited cols require q_len/t_len >= 3"
    );
    debug_assert!(q_len <= MAX_EXT, "q_len exceeds precomputed index capacity");
    debug_assert!(t_len <= MAX_EXT, "t_len exceeds precomputed index capacity");
    unsafe {
        let tj1 = *t_ptr.add(1);
        let tj2 = *t_ptr.add(2);

        // Rolling values for rows (k-1) at columns 1 and 2.
        let mut m1_prev = cell_get(m_ptr, width, 2, 1); // m(k-1,1)
        let mut bq1_prev = cell_get(bq_ptr, width, 2, 1); // bq(k-1,1)
        let mut m2_prev = cell_get(m_ptr, width, 2, 2); // m(k-1,2)
        let mut bq2_prev = cell_get(bq_ptr, width, 2, 2); // bq(k-1,2)

        for k in 3..q_len {
            let qi = *q_ptr.add(k);
            let qi_prev = *q_ptr.add(k - 1);

            // Primary[ k,1 ] = Bq
            let bq1 = best2(
                m1_prev,
                stack_dir::<LEFT, M>(qi_prev, qi, tj1, GAP, view.penalty),
                bq1_prev,
                stack_dir::<LEFT, M>(qi_prev, qi, GAP, GAP, view.penalty),
            );
            cell_set(bq_ptr, width, k, 1, bq1);

            // M[ k,2 ]
            let val = best2(
                m1_prev,
                stack_dir::<LEFT, M>(qi_prev, qi, tj1, tj2, view.penalty),
                bq1_prev,
                stack_dir::<LEFT, M>(qi_prev, qi, GAP, tj2, view.penalty),
            );
            cell_set(m_ptr, width, k, 2, val);
            update_best_with_term(best_e, best_i, best_j, val, view.terminal(k, 2), k, 2);

            // Secondary[ k,2 ] = Bt
            let m_k1 = cell_get(m_ptr, width, k, 1); // Preinitialized in frontier.
            cell_set(
                bt_ptr,
                width,
                k,
                2,
                add_e(m_k1, stack_dir::<LEFT, M>(qi, GAP, tj1, tj2, view.penalty)),
            );

            // Primary[ k,2 ] = Bq
            let bq2 = best2(
                m2_prev,
                stack_dir::<LEFT, M>(qi_prev, qi, tj2, GAP, view.penalty),
                bq2_prev,
                stack_dir::<LEFT, M>(qi_prev, qi, GAP, GAP, view.penalty),
            );
            cell_set(bq_ptr, width, k, 2, bq2);

            // Advance rolling state.
            m1_prev = m_k1;
            bq1_prev = bq1;
            m2_prev = val;
            bq2_prev = bq2;
        }
    }
}
