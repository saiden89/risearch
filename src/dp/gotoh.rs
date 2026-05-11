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

use crate::dsm::ScoringModel;
use crate::types::GAP;
use log::trace;

use super::{BestScore, DpGrid, DpView, ExtendDir, MAX_EXT};

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

    /// M ← M transition score: continue the paired stack.
    #[inline(always)]
    pub(crate) fn stack(&self, qp: u8, qc: u8, tp: u8, tc: u8) -> i32 {
        self.model.transition_score(qp, qc, tp, tc)
    }

    /// M ← Bq transition score: close a query gap and return to the stack.
    #[inline(always)]
    pub(crate) fn close_query_gap(&self, qp: u8, qc: u8, tc: u8) -> i32 {
        self.model.transition_score(qp, qc, GAP, tc)
    }

    /// M ← Bt transition score: close a target gap and return to the stack.
    #[inline(always)]
    pub(crate) fn close_target_gap(&self, qc: u8, tp: u8, tc: u8) -> i32 {
        self.model.transition_score(GAP, qc, tp, tc)
    }

    /// Bq ← M transition score: open a query gap.
    #[inline(always)]
    pub(crate) fn open_query_gap(&self, qp: u8, qc: u8, tc: u8) -> i32 {
        self.model.transition_score(qp, qc, tc, GAP)
    }

    /// Bq ← Bq transition score: extend a query gap.
    #[inline(always)]
    pub(crate) fn extend_query_gap(&self, qp: u8, qc: u8) -> i32 {
        self.model.transition_score(qp, qc, GAP, GAP)
    }

    /// Bt ← M transition score: open a target gap.
    #[inline(always)]
    pub(crate) fn open_target_gap(&self, qc: u8, tp: u8, tc: u8) -> i32 {
        self.model.transition_score(qc, GAP, tp, tc)
    }

    /// Bt ← Bt transition score: extend a target gap.
    #[inline(always)]
    pub(crate) fn extend_target_gap(&self, tp: u8, tc: u8) -> i32 {
        self.model.transition_score(GAP, GAP, tp, tc)
    }

    /// Terminal (boundary) penalty.
    #[inline(always)]
    pub(crate) fn terminal(&self, qc: u8, tc: u8) -> i32 {
        self.model.transition_score(qc, GAP, tc, GAP)
    }

    #[cfg_attr(feature = "prof", inline(never))]
    /// Run DP forward pass over `view`, reusing the caller-provided grid.
    pub fn extend(&self, view: &DpView<'_>, grid: &mut DpGrid) -> BestScore {
        let (q_len, t_len) = (view.q_len.min(MAX_EXT), view.t_len.min(MAX_EXT));

        trace!("{} q_len={} t_len={}", view.dir, q_len, t_len);
        debug_assert!(
            q_len > 0 && t_len > 0,
            "DP view must define a non-empty window"
        );

        let mut q_idx = std::mem::MaybeUninit::<[u8; MAX_EXT]>::uninit();
        let mut t_idx = std::mem::MaybeUninit::<[u8; MAX_EXT]>::uninit();

        let q_ptr = q_idx.as_mut_ptr() as *mut u8;
        let t_ptr = t_idx.as_mut_ptr() as *mut u8;

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

        // Initial score: terminal penalty for seed boundary
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

        self.dp_main_loop(q_ptr, t_ptr, grid, q_len, t_len, &mut best);

        trace!(
            "{} result: energy={} q_idx={} t_idx={}",
            view.dir,
            best.energy,
            best.q_idx,
            best.t_idx,
        );

        best
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dsm::DsmRegistry;
    use crate::types::{DsmId, Energy, BASE_COUNT};

    #[test]
    fn gotoh_matches_scoring_table() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let table = ScoringModel::new(&source, init, Energy::from_kcal(0.005));
        let gotoh = Gotoh::new(&table);

        for q1 in 0..BASE_COUNT as u8 {
            for q2 in 0..BASE_COUNT as u8 {
                for t1 in 0..BASE_COUNT as u8 {
                    for t2 in 0..BASE_COUNT as u8 {
                        assert_eq!(
                            gotoh.stack(q1, q2, t1, t2),
                            table.transition_score(q1, q2, t1, t2),
                            "stack mismatch at ({},{},{},{})",
                            q1,
                            q2,
                            t1,
                            t2
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn gotoh_semantic_transitions_match_scoring_model() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let table = ScoringModel::new(&source, init, Energy::from_kcal(0.005));
        let gotoh = Gotoh::new(&table);

        for qp in 0..BASE_COUNT as u8 {
            for qc in 0..BASE_COUNT as u8 {
                for tc in 0..BASE_COUNT as u8 {
                    assert_eq!(
                        gotoh.close_query_gap(qp, qc, tc),
                        table.transition_score(qp, qc, GAP, tc),
                    );
                    assert_eq!(
                        gotoh.open_query_gap(qp, qc, tc),
                        table.transition_score(qp, qc, tc, GAP),
                    );
                }
                assert_eq!(
                    gotoh.extend_query_gap(qp, qc),
                    table.transition_score(qp, qc, GAP, GAP),
                );
            }
        }

        for qc in 0..BASE_COUNT as u8 {
            for tp in 0..BASE_COUNT as u8 {
                for tc in 0..BASE_COUNT as u8 {
                    assert_eq!(
                        gotoh.close_target_gap(qc, tp, tc),
                        table.transition_score(GAP, qc, tp, tc),
                    );
                    assert_eq!(
                        gotoh.open_target_gap(qc, tp, tc),
                        table.transition_score(qc, GAP, tp, tc),
                    );
                }
            }
            for tc in 0..BASE_COUNT as u8 {
                assert_eq!(
                    gotoh.terminal(qc, tc),
                    table.transition_score(qc, GAP, tc, GAP),
                );
            }
        }

        for tp in 0..BASE_COUNT as u8 {
            for tc in 0..BASE_COUNT as u8 {
                assert_eq!(
                    gotoh.extend_target_gap(tp, tc),
                    table.transition_score(GAP, GAP, tp, tc),
                );
            }
        }
    }
}
