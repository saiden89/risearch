//! Gotoh 3-state DP recurrence: precomputed transition tables.
//!
//! All transitions are materialized at construction time from a `ScoringModel`.
//! The DP core asks rank-byte transition questions, not arbitrary 4-base tensor
//! lookups:
//!
//! - `M`: continue the paired stack, close a query gap, or close a target gap
//! - `Bq`: open or extend a query gap
//! - `Bt`: open or extend a target gap
//!
//! A mismatch is not its own transition family; it is simply an unfavorable
//! `stack` score.
//!
//! # Index safety
//!
//! All index arguments to point-lookup and slice-accessor methods must be in
//! `0..6` (valid DP lookup indices derived from the `#[repr(u8)]` `Base` enum).
//! `Gotoh::extend()` materializes those dense indices once from the semantic
//! `DpView` before entering the hot DP kernels.

use smallvec::SmallVec;

use crate::dsm::{RowLookup, ScoringModel};
use log::trace;

use super::{is_valid_score, BestScore, DpGrid, DpView, ExtendDir, TraceOp, MAX_EXT};

#[inline(always)]
fn is_transition(val: i32, pred: i32, energy: i32) -> bool {
    is_valid_score(pred) && val == pred + energy
}

/// Precomputed transition tables for the Gotoh 3-state DP.
///
/// Direction-agnostic: each instance corresponds to one `ScoringModel` orientation.
pub struct Gotoh {
    pub(crate) model: ScoringModel,
}

impl Gotoh {
    /// Materialize all transition tables from a `ScoringModel`.
    pub fn new(table: &ScoringModel) -> Self {
        Self {
            model: table.clone(),
        }
    }

    /// Build a row-local lookup for the given (qp, qc) pair.
    ///
    /// # Safety
    ///
    /// `qp` and `qc` must be in `0..BASE_COUNT` (valid base rank indices).
    #[inline(always)]
    pub(crate) fn row_lookup(&self, qp: u8, qc: u8) -> RowLookup {
        self.model.row_lookup(qp, qc)
    }

    /// M ← M transition score: continue the paired stack.
    #[inline(always)]
    pub(crate) fn stack(&self, qp: u8, qc: u8, tp: u8, tc: u8) -> i32 {
        self.model.stack_score(qp, qc, tp, tc)
    }

    /// M ← Bq transition score: close a query gap and return to the stack.
    #[inline(always)]
    pub(crate) fn close_query_gap(&self, qp: u8, qc: u8, tc: u8) -> i32 {
        self.model.close_query_gap_score(qp, qc, tc)
    }

    /// M ← Bt transition score: close a target gap and return to the stack.
    #[inline(always)]
    pub(crate) fn close_target_gap(&self, qc: u8, tp: u8, tc: u8) -> i32 {
        self.model.close_target_gap_score(qc, tp, tc)
    }

    /// Bq ← M transition score: open a query gap.
    #[inline(always)]
    pub(crate) fn open_query_gap(&self, qp: u8, qc: u8, tc: u8) -> i32 {
        self.model.open_query_gap_score(qp, qc, tc)
    }

    /// Bq ← Bq transition score: extend a query gap.
    #[inline(always)]
    pub(crate) fn extend_query_gap(&self, qp: u8, qc: u8) -> i32 {
        self.model.extend_query_gap_score(qp, qc)
    }

    /// Bt ← M transition score: open a target gap.
    #[inline(always)]
    pub(crate) fn open_target_gap(&self, qc: u8, tp: u8, tc: u8) -> i32 {
        self.model.open_target_gap_score(qc, tp, tc)
    }

    /// Bt ← Bt transition score: extend a target gap.
    #[inline(always)]
    pub(crate) fn extend_target_gap(&self, tp: u8, tc: u8) -> i32 {
        self.model.extend_target_gap_score(tp, tc)
    }

    /// Terminal (boundary) penalty.
    #[inline(always)]
    pub(crate) fn terminal(&self, qc: u8, tc: u8) -> i32 {
        self.model.terminal_penalty(qc, tc)
    }

    #[cfg_attr(feature = "prof", inline(never))]
    /// Run DP forward pass over `view`, reusing caller-provided grid and buffers.
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

        let q_ptr = q_buf.as_mut_ptr();
        let t_ptr = t_buf.as_mut_ptr();

        match view.dir {
            ExtendDir::Left => {
                for i in 0..q_len {
                    unsafe { *q_ptr.add(i) = view.query.get_unchecked(view.q_anchor - i).as_u8() };
                }
                for j in 0..t_len {
                    unsafe { *t_ptr.add(j) = view.target.get_unchecked(view.t_anchor + j).as_u8() };
                }
            }
            ExtendDir::Right => {
                for i in 0..q_len {
                    unsafe { *q_ptr.add(i) = view.query.get_unchecked(view.q_anchor + i).as_u8() };
                }
                for j in 0..t_len {
                    unsafe { *t_ptr.add(j) = view.target.get_unchecked(view.t_anchor - j).as_u8() };
                }
            }
        }

        let q0 = unsafe { *q_ptr };
        let t0 = unsafe { *t_ptr };
        let mut best = BestScore::new(self.terminal(q0, t0));

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

                    let stack = self.stack(qi_prev, qi, tj_prev, tj);
                    let close_query_gap = self.close_query_gap(qi_prev, qi, tj);
                    let close_target_gap = self.close_target_gap(qi, tj_prev, tj);

                    i -= 1;
                    j -= 1;

                    if is_transition(c.m, diag.m, stack) {
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

                    let open_query_gap = self.open_query_gap(qi_prev, qi, tj);
                    let extend_query_gap = self.extend_query_gap(qi_prev, qi);

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

                    let open_target_gap = self.open_target_gap(qi, tj_prev, tj);
                    let extend_target_gap = self.extend_target_gap(tj_prev, tj);

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
