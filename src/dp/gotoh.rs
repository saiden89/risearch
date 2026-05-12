//! Gotoh 3-state DP recurrence.
//!
//! This module implements the core DP engine for directional extension.
//! The algorithm is agnostic to the underlying scoring data, interacting
//! only through the [`GotohScoring`] trait.
//!
//! The DP core identifies three transition families:
//! - `M`: continue the paired region, close a query gap, or close a target gap
//! - `Bq`: open or extend a query gap
//! - `Bt`: open or extend a target gap
//!
//! A mismatch is not its own transition family; it is simply an unfavorable
//! paired transition score.
//!
//! # Index safety
//!
//! All index arguments to point-lookup and slice-accessor methods must be in
//! the valid symbol index range of the scoring source. `Gotoh::extend()`
//! materializes those dense indices once from the semantic `DpView` before
//! entering the hot DP kernels.

use crate::dp::scoring::GotohScoring;
use smallvec::SmallVec;

use log::trace;

use super::{is_valid_score, BestScore, DpGrid, DpView, TraceOp, MAX_EXT};

#[inline(always)]
fn is_transition(val: i32, pred: i32, energy: i32) -> bool {
    is_valid_score(pred) && val == pred + energy
}

/// A Gotoh 3-state DP engine.
///
/// Direction-agnostic: each instance corresponds to one orientation of the
/// underlying scoring source.
pub struct Gotoh<S: GotohScoring> {
    pub(crate) scoring: S,
}

impl<S: GotohScoring> Gotoh<S> {
    /// Create a new Gotoh engine for the given scoring source.
    pub fn new(scoring: &S) -> Self
    where
        S: Clone,
    {
        Self {
            scoring: scoring.clone(),
        }
    }

    /// Boundary penalty for the given symbol pair.
    #[inline(always)]
    pub fn boundary(&self, qc: u8, tc: u8) -> i32 {
        self.scoring.boundary(qc, tc)
    }

    #[cfg_attr(feature = "prof", inline(never))]
    /// DP forward pass over `view`, reusing caller-provided grid and buffers.
    pub fn extend(
        &self,
        view: &DpView<'_>,
        grid: &mut DpGrid,
        q_buf: &mut [u8],
        t_buf: &mut [u8],
    ) -> BestScore {
        let (q_len, t_len) = (view.q_len.min(MAX_EXT), view.t_len.min(MAX_EXT));

        trace!("{} q_len={} t_len={}", view.dir, q_len, t_len);
        debug_assert!(
            q_len > 0 && t_len > 0,
            "DP view must define a non-empty window"
        );

        view.fill_buffers(q_buf, t_buf);

        let q_ptr = q_buf.as_ptr();
        let t_ptr = t_buf.as_ptr();

        let q0 = unsafe { *q_ptr };
        let t0 = unsafe { *t_ptr };
        let mut best = BestScore::new(self.scoring.boundary(q0, t0));

        if q_len <= 1 || t_len <= 1 {
            return best;
        }

        grid.resize(t_len + 1, q_len + 1);

        let has_main_region = self.init_frontier(q_ptr, t_ptr, grid, q_len, t_len, &mut best);
        if !has_main_region {
            return best;
        }

        self.main_loop(q_ptr, t_ptr, grid, q_len, t_len, &mut best);

        trace!(
            "{} result: energy={} q_idx={} t_idx={}",
            view.dir,
            best.energy,
            best.q_idx,
            best.t_idx,
        );

        best
    }

    /// Walk the filled grid backward from `(end_i, end_j)` to recover the
    /// state-transition path.
    pub fn traceback(
        &self,
        q_buf: &[u8],
        t_buf: &[u8],
        grid: &DpGrid,
        end_i: usize,
        end_j: usize,
    ) -> SmallVec<[TraceOp; 64]> {
        let (mut i, mut j) = (end_i, end_j);
        let mut state = TraceOp::Paired;
        let mut out = SmallVec::new();

        while i > 0 || j > 0 {
            let next = match state {
                TraceOp::Paired if i > 0 && j > 0 => {
                    out.push(TraceOp::Paired);

                    let c = grid.get(i, j);
                    let diag = grid.get(i - 1, j - 1);

                    let qi_prev = q_buf[i - 1];
                    let qi = q_buf[i];
                    let tj_prev = t_buf[j - 1];
                    let tj = t_buf[j];

                    let r#match = self.scoring.r#match(qi_prev, qi, tj_prev, tj);
                    let close_query_gap = self.scoring.close_query_gap(qi_prev, qi, tj);
                    let close_target_gap = self.scoring.close_target_gap(qi, tj_prev, tj);

                    i -= 1;
                    j -= 1;

                    if is_transition(c.m, diag.m, r#match) {
                        TraceOp::Paired
                    } else if is_transition(c.m, diag.bq, close_query_gap) {
                        TraceOp::GapQ
                    } else if is_transition(c.m, diag.bt, close_target_gap) {
                        TraceOp::GapT
                    } else {
                        debug_assert!(
                            false,
                            "traceback: no valid predecessor at ({i},{j}) in Paired state"
                        );
                        break;
                    }
                }
                TraceOp::GapQ if i > 0 => {
                    out.push(TraceOp::GapQ);

                    let c = grid.get(i, j);
                    let up = grid.get(i - 1, j);

                    let qi_prev = q_buf[i - 1];
                    let qi = q_buf[i];
                    let tj = t_buf[j];

                    let open_query_gap = self.scoring.open_query_gap(qi_prev, qi, tj);
                    let extend_query_gap = self.scoring.extend_query_gap(qi_prev, qi);

                    i -= 1;

                    if is_transition(c.bq, up.m, open_query_gap) {
                        TraceOp::Paired
                    } else if is_transition(c.bq, up.bq, extend_query_gap) {
                        TraceOp::GapQ
                    } else {
                        debug_assert!(
                            false,
                            "traceback: no valid predecessor at ({i},{j}) in GapQ state"
                        );
                        break;
                    }
                }
                TraceOp::GapT if j > 0 => {
                    out.push(TraceOp::GapT);

                    let c = grid.get(i, j);
                    let left = grid.get(i, j - 1);

                    let qi = q_buf[i];
                    let tj_prev = t_buf[j - 1];
                    let tj = t_buf[j];

                    let open_target_gap = self.scoring.open_target_gap(qi, tj_prev, tj);
                    let extend_target_gap = self.scoring.extend_target_gap(tj_prev, tj);

                    j -= 1;

                    if is_transition(c.bt, left.m, open_target_gap) {
                        TraceOp::Paired
                    } else if is_transition(c.bt, left.bt, extend_target_gap) {
                        TraceOp::GapT
                    } else {
                        debug_assert!(
                            false,
                            "traceback: no valid predecessor at ({i},{j}) in GapT state"
                        );
                        break;
                    }
                }
                _ => break,
            };
            state = next;
        }

        out
    }
}

#[cfg(test)]
mod tests {
    // Note: Integration tests for the DP algorithm (parity with legacy RIsearch)
    // are located in the /tests directory.
}
