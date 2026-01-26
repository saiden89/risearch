use crate::dsm::dsm_lookup_raw;

use super::{DpView, ExtendDir, GAP, MAX_EXT, MIN_SCORE};

// =============================================================================
// INIT HELPERS - Reduce code duplication in DP initialization
// =============================================================================

/// Pick best value from two sources.
/// Uses std::cmp::max which LLVM compiles to branchless CMOV.
#[inline(always)]
pub(super) fn pick_best(val_a: i32, val_b: i32) -> i32 {
    std::cmp::max(val_a, val_b)
}

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

/// max of 2 values - branchless via std::cmp::max (compiles to CMOV)
#[inline(always)]
pub(super) fn max2(a: i32, b: i32) -> i32 {
    std::cmp::max(a, b)
}

/// max of 3 values - branchless
#[inline(always)]
pub(super) fn max3(a: i32, b: i32, c: i32) -> i32 {
    std::cmp::max(std::cmp::max(a, b), c)
}

/// Initialize limited rows (t_len axis) - score-only version.
/// All writes are unconditional to avoid reading stale data.
#[inline(always)]
pub(super) fn init_limited_rows(
    view: &DpView<'_>,
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
    let left = view.dir == ExtendDir::Left;
    let stack = |q1: usize, q2: usize, t1: usize, t2: usize| -> i32 {
        if left {
            dsm_lookup_raw(q2, q1, t2, t1)
        } else {
            dsm_lookup_raw(q1, q2, t1, t2)
        }
    };

    debug_assert!(
        q_len > 2 && t_len > 2,
        "limited rows require q_len/t_len >= 3"
    );
    debug_assert!(q_len <= MAX_EXT, "q_len exceeds precomputed index capacity");
    debug_assert!(t_len <= MAX_EXT, "t_len exceeds precomputed index capacity");
    unsafe {
        let idx = |i: usize, j: usize| -> usize { i * width + j };
        let qi1 = *q_ptr.add(1);
        let qi2 = *q_ptr.add(2);

        // Rolling values for row 1 and 2
        let mut m1_prev = *m_ptr.add(idx(1, 2)); // m(1, k-1)
        let mut bt1_prev = *bt_ptr.add(idx(1, 2)); // bt(1, k-1)
        let mut m2_prev = *m_ptr.add(idx(2, 2)); // m(2, k-1)
        let mut bt2_prev = *bt_ptr.add(idx(2, 2)); // bt(2, k-1)

        for k in 3..t_len {
            let tj = *t_ptr.add(k);
            let tj_prev = *t_ptr.add(k - 1);

            // Primary[1,k] = Bt
            let from_m = add_e(m1_prev, stack(qi1, GAP, tj_prev, tj));
            let from_b = add_e(bt1_prev, stack(GAP, GAP, tj_prev, tj));
            let bt1 = pick_best(from_m, from_b);
            *bt_ptr.add(idx(1, k)) = bt1;

            // M[2,k]
            let from_m = add_e(m1_prev, stack(qi1, qi2, tj_prev, tj));
            let from_b = add_e(bt1_prev, stack(GAP, qi2, tj_prev, tj));
            let m2 = pick_best(from_m, from_b);
            *m_ptr.add(idx(2, k)) = m2;
            let term = if left {
                dsm_lookup_raw(GAP, qi2, GAP, tj)
            } else {
                dsm_lookup_raw(qi2, GAP, tj, GAP)
            };
            update_best_with_term(best_e, best_i, best_j, m2, term, 2, k);

            // Secondary[2,k] = Bq
            let m1k = *m_ptr.add(idx(1, k));
            *bq_ptr.add(idx(2, k)) = add_e(m1k, stack(qi1, qi2, tj, GAP));

            // Primary[2,k] = Bt
            let from_m = add_e(m2_prev, stack(qi2, GAP, tj_prev, tj));
            let from_b = add_e(bt2_prev, stack(GAP, GAP, tj_prev, tj));
            let bt2 = pick_best(from_m, from_b);
            *bt_ptr.add(idx(2, k)) = bt2;

            m1_prev = *m_ptr.add(idx(1, k));
            bt1_prev = bt1;
            m2_prev = m2;
            bt2_prev = bt2;
        }
    }
}

/// Initialize limited columns (q_len axis) - score-only version.
/// All writes are unconditional to avoid reading stale data.
#[inline(always)]
pub(super) fn init_limited_cols(
    view: &DpView<'_>,
    q_ptr: *const usize,
    t_ptr: *const usize,
    m_ptr: *mut i32,
    bq_ptr: *mut i32,
    bt_ptr: *mut i32,
    width: usize,
    q_len: usize,
    best_e: &mut i32,
    best_i: &mut usize,
    best_j: &mut usize,
) {
    // SAFETY invariants:
    // - q_ptr and t_ptr are valid for indices [0, q_len) and [0, t_len)
    // - q_len and t_len are >= 3 (we index 1 and 2)
    // - m_ptr/bq_ptr/bt_ptr point to matrices sized at least (q_len+1) x (t_len+1)
    let left = view.dir == ExtendDir::Left;
    let stack = |q1: usize, q2: usize, t1: usize, t2: usize| -> i32 {
        if left {
            dsm_lookup_raw(q2, q1, t2, t1)
        } else {
            dsm_lookup_raw(q1, q2, t1, t2)
        }
    };

    debug_assert!(q_len > 2, "limited cols require q_len >= 3");
    debug_assert!(q_len <= MAX_EXT, "q_len exceeds precomputed index capacity");
    unsafe {
        let idx = |i: usize, j: usize| -> usize { i * width + j };
        let tj1 = *t_ptr.add(1);
        let tj2 = *t_ptr.add(2);

        for k in 3..q_len {
            let qi = *q_ptr.add(k);
            let qi_prev = *q_ptr.add(k - 1);

            // Primary[ k,1 ] = Bq
            let from_m = add_e(*m_ptr.add(idx(k - 1, 1)), stack(qi_prev, qi, tj1, GAP));
            let from_b = add_e(*bq_ptr.add(idx(k - 1, 1)), stack(qi_prev, qi, GAP, GAP));
            *bq_ptr.add(idx(k, 1)) = pick_best(from_m, from_b);

            // M[ k,2 ]
            let from_m = add_e(*m_ptr.add(idx(k - 1, 1)), stack(qi_prev, qi, tj1, tj2));
            let from_b = add_e(*bq_ptr.add(idx(k - 1, 1)), stack(qi_prev, qi, GAP, tj2));
            let val = pick_best(from_m, from_b);
            *m_ptr.add(idx(k, 2)) = val;
            let term = if left {
                dsm_lookup_raw(GAP, qi, GAP, tj2)
            } else {
                dsm_lookup_raw(qi, GAP, tj2, GAP)
            };
            update_best_with_term(best_e, best_i, best_j, val, term, k, 2);

            // Secondary[ k,2 ] = Bt
            let m_k1 = *m_ptr.add(idx(k, 1));
            *bt_ptr.add(idx(k, 2)) = add_e(m_k1, stack(qi, GAP, tj1, tj2));

            // Primary[ k,2 ] = Bq
            let from_m = add_e(*m_ptr.add(idx(k - 1, 2)), stack(qi_prev, qi, tj2, GAP));
            let from_b = add_e(*bq_ptr.add(idx(k - 1, 2)), stack(qi_prev, qi, GAP, GAP));
            *bq_ptr.add(idx(k, 2)) = pick_best(from_m, from_b);
        }
    }
}
