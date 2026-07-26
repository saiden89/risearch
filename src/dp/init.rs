use std::cmp::max;

use crate::dp::gotoh::Gotoh;
use crate::dp::scoring::GotohScoring;

use super::{BestScore, DpCell, DpGrid};

impl<S: GotohScoring> Gotoh<S> {
    /// Initialize DP boundary cells.
    ///
    /// Returns `true` when the main DP region (`i >= 3`, `j >= 3`) exists.
    /// Returns `false` when the matrix is too small for a main-loop pass.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    pub(super) fn init_frontier(
        &self,
        q_ptr: *const u8,
        t_ptr: *const u8,
        grid: &mut DpGrid,
        q_len: usize,
        t_len: usize,
        best: &mut BestScore,
    ) -> bool {
        let ptr = grid.as_mut_ptr();
        let width = grid.width();
        let scoring = &self.scoring;

        debug_assert!(width > t_len, "matrix width too small for t_len");

        // SAFETY: q_ptr valid for [0, q_len), t_ptr for [0, t_len), values in symbol index range.
        // Grid is (q_len+1) x (t_len+1) via Gotoh::extend(). All writes below are
        // to cells (i, j) with i < q_len, j < t_len — within the allocation.
        unsafe {
            let q0 = *q_ptr;
            let q1 = *q_ptr.add(1);
            let t0 = *t_ptr;
            let t1 = *t_ptr.add(1);

            let gt01 = scoring.open_target_gap(q0, t0, t1);
            let gq10 = scoring.open_query_gap(q0, q1, t0);
            let m11 = scoring.r#match(q0, q1, t0, t1);

            *ptr = DpCell {
                m: 0,
                ..DpCell::EMPTY
            };
            *ptr.add(1) = DpCell {
                gap_t: gt01,
                ..DpCell::EMPTY
            };
            *ptr.add(width) = DpCell {
                gap_q: gq10,
                ..DpCell::EMPTY
            };
            *ptr.add(width + 1) = DpCell {
                m: m11,
                ..DpCell::EMPTY
            };
            best.update_if_better(m11, scoring.boundary(q1, t1), 1, 1);

            // Tiny-matrix path: no row/col 2 exists
            if q_len <= 2 || t_len <= 2 {
                let mut gt_prev = gt01;
                let mut t_prev = t1;
                for j in 2..t_len {
                    let tc = *t_ptr.add(j);
                    let gt = gt_prev + scoring.extend_target_gap(t_prev, tc);
                    let m1 = gt_prev + scoring.close_target_gap(q1, t_prev, tc);
                    *ptr.add(j) = DpCell {
                        gap_t: gt,
                        ..DpCell::EMPTY
                    };
                    *ptr.add(width + j) = DpCell {
                        m: m1,
                        ..DpCell::EMPTY
                    };
                    best.update_if_better(m1, scoring.boundary(q1, tc), 1, j);
                    gt_prev = gt;
                    t_prev = tc;
                }

                let mut gq_prev = gq10;
                let mut q_prev = q1;
                for i in 2..q_len {
                    let qc = *q_ptr.add(i);
                    let gq = gq_prev + scoring.extend_query_gap(q_prev, qc);
                    let m1 = gq_prev + scoring.close_query_gap(q_prev, qc, t1);
                    *ptr.add(i * width) = DpCell {
                        gap_q: gq,
                        ..DpCell::EMPTY
                    };
                    *ptr.add(i * width + 1) = DpCell {
                        m: m1,
                        ..DpCell::EMPTY
                    };
                    best.update_if_better(m1, scoring.boundary(qc, t1), i, 1);
                    gq_prev = gq;
                    q_prev = qc;
                }

                return false;
            }

            // Full path: q_len >= 3 and t_len >= 3
            let q2 = *q_ptr.add(2);
            let t2 = *t_ptr.add(2);

            // 3x3 corner extension (5 remaining cells)
            let gt02 = gt01 + scoring.extend_target_gap(t1, t2);
            let m12 = gt01 + scoring.close_target_gap(q1, t1, t2);
            let gt12 = m11 + scoring.open_target_gap(q1, t1, t2);

            let gq20 = gq10 + scoring.extend_query_gap(q1, q2);
            let m21 = gq10 + scoring.close_query_gap(q1, q2, t1);
            let gq21 = m11 + scoring.open_query_gap(q1, q2, t1);

            let m22 = m11 + scoring.r#match(q1, q2, t1, t2);
            let gq22 = m12 + scoring.open_query_gap(q1, q2, t2);
            let gt22 = m21 + scoring.open_target_gap(q2, t1, t2);

            *ptr.add(2) = DpCell {
                gap_t: gt02,
                ..DpCell::EMPTY
            };
            *ptr.add(width + 2) = DpCell {
                m: m12,
                gap_t: gt12,
                ..DpCell::EMPTY
            };
            *ptr.add(2 * width) = DpCell {
                gap_q: gq20,
                ..DpCell::EMPTY
            };
            *ptr.add(2 * width + 1) = DpCell {
                m: m21,
                gap_q: gq21,
                ..DpCell::EMPTY
            };
            *ptr.add(2 * width + 2) = DpCell {
                m: m22,
                gap_q: gq22,
                gap_t: gt22,
            };

            best.update_if_better(m12, scoring.boundary(q1, t2), 1, 2);
            best.update_if_better(m21, scoring.boundary(q2, t1), 2, 1);
            best.update_if_better(m22, scoring.boundary(q2, t2), 2, 2);

            // Top-rows fused loop: rows 0, 1, 2 for j >= 3
            let mut gt0_prev = gt02;
            let mut m1_prev = m12;
            let mut gt1_prev = gt12;
            let mut m2_prev = m22;
            let mut gt2_prev = gt22;
            let mut t_prev = t2;

            for k in 3..t_len {
                let tj = *t_ptr.add(k);
                let ext_t = scoring.extend_target_gap(t_prev, tj);

                let gt0 = gt0_prev + ext_t;
                let m1 = gt0_prev + scoring.close_target_gap(q1, t_prev, tj);
                let gt1 = max(
                    m1_prev + scoring.open_target_gap(q1, t_prev, tj),
                    gt1_prev + ext_t,
                );
                let m2 = max(
                    m1_prev + scoring.r#match(q1, q2, t_prev, tj),
                    gt1_prev + scoring.close_target_gap(q2, t_prev, tj),
                );
                let gq2 = m1 + scoring.open_query_gap(q1, q2, tj);
                let gt2 = max(
                    m2_prev + scoring.open_target_gap(q2, t_prev, tj),
                    gt2_prev + ext_t,
                );

                *ptr.add(k) = DpCell {
                    gap_t: gt0,
                    ..DpCell::EMPTY
                };
                *ptr.add(width + k) = DpCell {
                    m: m1,
                    gap_t: gt1,
                    ..DpCell::EMPTY
                };
                *ptr.add(2 * width + k) = DpCell {
                    m: m2,
                    gap_q: gq2,
                    gap_t: gt2,
                };

                best.update_if_better(m1, scoring.boundary(q1, tj), 1, k);
                best.update_if_better(m2, scoring.boundary(q2, tj), 2, k);

                gt0_prev = gt0;
                m1_prev = m1;
                gt1_prev = gt1;
                m2_prev = m2;
                gt2_prev = gt2;
                t_prev = tj;
            }

            // Left-cols fused loop: cols 0, 1, 2 for i >= 3
            let mut gq0_prev = gq20;
            let mut m1_prev = m21;
            let mut gq1_prev = gq21;
            let mut m2_prev = m22;
            let mut gq2_prev = gq22;
            let mut q_prev = q2;

            for k in 3..q_len {
                let qi = *q_ptr.add(k);
                let ext_q = scoring.extend_query_gap(q_prev, qi);

                let gq0 = gq0_prev + ext_q;
                let m1 = gq0_prev + scoring.close_query_gap(q_prev, qi, t1);
                let gq1 = max(
                    m1_prev + scoring.open_query_gap(q_prev, qi, t1),
                    gq1_prev + ext_q,
                );
                let m2 = max(
                    m1_prev + scoring.r#match(q_prev, qi, t1, t2),
                    gq1_prev + scoring.close_query_gap(q_prev, qi, t2),
                );
                let gt2 = m1 + scoring.open_target_gap(qi, t1, t2);
                let gq2 = max(
                    m2_prev + scoring.open_query_gap(q_prev, qi, t2),
                    gq2_prev + ext_q,
                );

                *ptr.add(k * width) = DpCell {
                    gap_q: gq0,
                    ..DpCell::EMPTY
                };
                *ptr.add(k * width + 1) = DpCell {
                    m: m1,
                    gap_q: gq1,
                    ..DpCell::EMPTY
                };
                *ptr.add(k * width + 2) = DpCell {
                    m: m2,
                    gap_q: gq2,
                    gap_t: gt2,
                };

                best.update_if_better(m1, scoring.boundary(qi, t1), k, 1);
                best.update_if_better(m2, scoring.boundary(qi, t2), k, 2);

                gq0_prev = gq0;
                m1_prev = m1;
                gq1_prev = gq1;
                m2_prev = m2;
                gq2_prev = gq2;
                q_prev = qi;
            }
        }

        true
    }
}
