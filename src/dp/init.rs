use std::cmp::max;

use crate::dp::gotoh::Gotoh;

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
///
/// # Safety
/// Caller must ensure `i * width + j` is within the grid allocation.
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
    gotoh: &Gotoh,
    best: &mut BestScore,
) -> bool {
    let ptr = grid.ptr();
    let width = grid.width();

    debug_assert!(width > t_len, "matrix width too small for t_len");

    // SAFETY: q_ptr valid for [0, q_len), t_ptr valid for [0, t_len), values ∈ 0..6.
    // Grid allocated as (q_len+1) × (t_len+1) by DpExtender::extend().
    // All cell() calls address (i, j) with i < q_len, j < t_len, within the grid.
    // All Gotoh calls receive Base::idx() values ∈ 0..6.
    unsafe {
        let q0 = *q_ptr.add(0);
        let q1 = *q_ptr.add(1);
        let t0 = *t_ptr.add(0);
        let t1 = *t_ptr.add(1);

        *cell(ptr, width, 0, 0) = DpCell {
            m: 0,
            bq: NEG_INF,
            bt: NEG_INF,
        };
        *cell(ptr, width, 0, 1) = DpCell {
            m: NEG_INF,
            bq: NEG_INF,
            bt: gotoh.bt_open(q0, t0, t1),
        };
        *cell(ptr, width, 1, 0) = DpCell {
            m: NEG_INF,
            bq: gotoh.bq_open(q0, q1, t0),
            bt: NEG_INF,
        };

        let m11 = gotoh.match_energy(q0, q1, t0, t1);
        *cell(ptr, width, 1, 1) = DpCell {
            m: m11,
            bq: NEG_INF,
            bt: NEG_INF,
        };
        best.update(m11, gotoh.terminal(q1, t1), 1, 1);

        // Row 0 (Bt only) and Row 1 (M) — k ∈ 2..t_len
        for k in 2..t_len {
            let t_prev = *t_ptr.add(k - 1);
            let t_curr = *t_ptr.add(k);
            let prev = (*cell(ptr, width, 0, k - 1)).bt;
            let bt_val = add_e(prev, gotoh.bt_extend_e(t_prev, t_curr));
            let m_val = add_e(prev, gotoh.m_from_bt(q1, t_prev, t_curr));
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
            best.update(m_val, gotoh.terminal(q1, t_curr), 1, k);
        }

        // Col 0 (Bq only) and Col 1 (M) — k ∈ 2..q_len
        for k in 2..q_len {
            let q_prev = *q_ptr.add(k - 1);
            let q_curr = *q_ptr.add(k);
            let prev = (*cell(ptr, width, k - 1, 0)).bq;
            let bq_val = add_e(prev, gotoh.bq_extend(q_prev, q_curr));
            let m_val = add_e(prev, gotoh.m_from_bq(q_prev, q_curr, t1));
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
            best.update(m_val, gotoh.terminal(q_curr, t1), k, 1);
        }
    }

    if q_len <= 2 || t_len <= 2 {
        return false;
    }

    // SAFETY: Same invariants as above. q_len > 2 and t_len > 2 guaranteed here.
    unsafe {
        let q1 = *q_ptr.add(1);
        let q2 = *q_ptr.add(2);
        let t1 = *t_ptr.add(1);
        let t2 = *t_ptr.add(2);
        let m11_val = (*cell(ptr, width, 1, 1)).m;
        let bt12 = add_e(m11_val, gotoh.bt_open(q1, t1, t2));
        let bq21 = add_e(m11_val, gotoh.bq_open(q1, q2, t1));
        let m22 = add_e(m11_val, gotoh.match_energy(q1, q2, t1, t2));
        (*cell(ptr, width, 1, 2)).bt = bt12;
        (*cell(ptr, width, 2, 1)).bq = bq21;
        (*cell(ptr, width, 2, 2)).m = m22;
        best.update(m22, gotoh.terminal(q2, t2), 2, 2);

        let m12 = (*cell(ptr, width, 1, 2)).m;
        let m21 = (*cell(ptr, width, 2, 1)).m;
        (*cell(ptr, width, 2, 2)).bq = add_e(m12, gotoh.bq_open(q1, q2, t2));
        (*cell(ptr, width, 2, 2)).bt = add_e(m21, gotoh.bt_open(q2, t1, t2));
    }

    init_limited_rows(q_ptr, t_ptr, grid, q_len, t_len, gotoh, best);
    init_limited_cols(q_ptr, t_ptr, grid, q_len, t_len, gotoh, best);

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
    gotoh: &Gotoh,
    best: &mut BestScore,
) {
    let ptr = grid.ptr();
    let width = grid.width();

    debug_assert!(
        q_len > 2 && t_len > 2,
        "limited rows require q_len/t_len >= 3"
    );

    // SAFETY: q_ptr/t_ptr valid for [0, q_len)/[0, t_len), values ∈ 0..6.
    // Grid allocated as (q_len+1) × (t_len+1). All cell accesses at rows 1-2,
    // cols 2..t_len are in-bounds.
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

            let bt1 = best2(
                m1_prev,
                gotoh.bt_open(qi1, tj_prev, tj),
                bt1_prev,
                gotoh.bt_extend_e(tj_prev, tj),
            );
            (*cell(ptr, width, 1, k)).bt = bt1;

            let m2 = best2(
                m1_prev,
                gotoh.match_energy(qi1, qi2, tj_prev, tj),
                bt1_prev,
                gotoh.m_from_bt(qi2, tj_prev, tj),
            );
            (*cell(ptr, width, 2, k)).m = m2;
            best.update(m2, gotoh.terminal(qi2, tj), 2, k);

            let m1k = (*cell(ptr, width, 1, k)).m;
            (*cell(ptr, width, 2, k)).bq = add_e(m1k, gotoh.bq_open(qi1, qi2, tj));

            let bt2 = best2(
                m2_prev,
                gotoh.bt_open(qi2, tj_prev, tj),
                bt2_prev,
                gotoh.bt_extend_e(tj_prev, tj),
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
    gotoh: &Gotoh,
    best: &mut BestScore,
) {
    let ptr = grid.ptr();
    let width = grid.width();

    debug_assert!(
        q_len > 2 && t_len > 2,
        "limited cols require q_len/t_len >= 3"
    );

    // SAFETY: q_ptr/t_ptr valid for [0, q_len)/[0, t_len), values ∈ 0..6.
    // Grid allocated as (q_len+1) × (t_len+1). All cell accesses at cols 1-2,
    // rows 2..q_len are in-bounds.
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

            let bq1 = best2(
                m1_prev,
                gotoh.bq_open(qi_prev, qi, tj1),
                bq1_prev,
                gotoh.bq_extend(qi_prev, qi),
            );
            (*cell(ptr, width, k, 1)).bq = bq1;

            let m2 = best2(
                m1_prev,
                gotoh.match_energy(qi_prev, qi, tj1, tj2),
                bq1_prev,
                gotoh.m_from_bq(qi_prev, qi, tj2),
            );
            (*cell(ptr, width, k, 2)).m = m2;
            best.update(m2, gotoh.terminal(qi, tj2), k, 2);

            let mk1 = (*cell(ptr, width, k, 1)).m;
            (*cell(ptr, width, k, 2)).bt = add_e(mk1, gotoh.bt_open(qi, tj1, tj2));

            let bq2 = best2(
                m2_prev,
                gotoh.bq_open(qi_prev, qi, tj2),
                bq2_prev,
                gotoh.bq_extend(qi_prev, qi),
            );
            (*cell(ptr, width, k, 2)).bq = bq2;

            m1_prev = mk1;
            bq1_prev = bq1;
            m2_prev = m2;
            bq2_prev = bq2;
        }
    }
}
