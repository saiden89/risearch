//! Gotoh 3-state DP recurrence: precomputed transition tables.
//!
//! All transitions are materialized at construction time from a `ScoringModel`.
//! The DP core asks semantic transition questions, not arbitrary 4-base tensor
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

use crate::dsm::{ScoringModel, GAP};
use crate::types::{Base, BASE_COUNT};
use log::trace;

use super::{BestScore, DpGrid, DpView, ExtendDir, GotohResult, MAX_EXT};

const BC: usize = BASE_COUNT; // 6

/// Precomputed transition tables for the Gotoh 3-state DP.
///
/// Direction-agnostic: each instance corresponds to one `ScoringModel` orientation.
pub struct Gotoh {
    /// M ← M: `stack[qp][qc][tp][tc]` — full 4D paired-stack continuation.
    ///
    /// Sized as 6⁴ rather than 5⁴ so that all 36-entry profile slices share the
    /// uniform stride `tp*6 + tc`. This lets the inner loop compute one index
    /// and reuse it across all target-pair-dependent transitions.
    stack: [i32; 1296],
    /// M ← Bq: `close_query_gap[qp][qc][tc]` — re-entry from a query gap.
    close_query_gap: [i32; 216],
    /// M ← Bt: `close_target_gap[qc][tp][tc]` — re-entry from a target gap.
    close_target_gap: [i32; 216],
    /// Bq ← M: `open_query_gap[qp][qc][tc]` — open a query gap.
    open_query_gap: [i32; 216],
    /// Bq ← Bq: `extend_query_gap[qp][qc]` — extend a query gap.
    extend_query_gap: [i32; 36],
    /// Bt ← M: `open_target_gap[qc][tp][tc]` — open a target gap.
    open_target_gap: [i32; 216],
    /// Bt ← Bt: `extend_target_gap[tp][tc]` — extend a target gap.
    extend_target_gap: [i32; 36],
    /// Terminal: `terminal[qc][tc]` — boundary penalty (36 entries)
    terminal: [i32; 36],
}

impl Gotoh {
    /// Materialize all transition tables from a `ScoringModel`.
    pub fn new(table: &ScoringModel) -> Self {
        let mut stack = [0i32; 1296];
        let mut close_query_gap = [0i32; 216];
        let mut close_target_gap = [0i32; 216];
        let mut open_query_gap = [0i32; 216];
        let mut extend_query_gap = [0i32; 36];
        let mut open_target_gap = [0i32; 216];
        let mut extend_target_gap = [0i32; 36];
        let mut terminal = [0i32; 36];

        for qp in 0..BC {
            for qc in 0..BC {
                let qp_qc_36 = qp * BC * BC * BC + qc * BC * BC;
                let qp_qc_6 = qp * BC + qc;
                for tp in 0..BC {
                    for tc in 0..BC {
                        stack[qp_qc_36 + tp * BC + tc] = table.transition_energy(qp, qc, tp, tc);
                    }
                }
                for tc in 0..BC {
                    close_query_gap[qp * BC * BC + qc * BC + tc] =
                        table.transition_energy(qp, qc, GAP, tc);
                    open_query_gap[qp * BC * BC + qc * BC + tc] =
                        table.transition_energy(qp, qc, tc, GAP);
                }
                extend_query_gap[qp_qc_6] = table.transition_energy(qp, qc, GAP, GAP);
            }
        }

        for qc in 0..BC {
            for tp in 0..BC {
                for tc in 0..BC {
                    close_target_gap[qc * BC * BC + tp * BC + tc] =
                        table.transition_energy(GAP, qc, tp, tc);
                    open_target_gap[qc * BC * BC + tp * BC + tc] =
                        table.transition_energy(qc, GAP, tp, tc);
                }
            }
            for tc in 0..BC {
                terminal[qc * BC + tc] = table.transition_energy(qc, GAP, tc, GAP);
            }
        }

        for tp in 0..BC {
            for tc in 0..BC {
                extend_target_gap[tp * BC + tc] = table.transition_energy(GAP, GAP, tp, tc);
            }
        }

        Self {
            stack,
            close_query_gap,
            close_target_gap,
            open_query_gap,
            extend_query_gap,
            open_target_gap,
            extend_target_gap,
            terminal,
        }
    }

    // =========================================================================
    // SEMANTIC TRANSITION QUERIES — used by traceback and boundary init
    // =========================================================================

    /// M ← M transition energy: continue the paired stack.
    #[inline(always)]
    pub(crate) fn stack(&self, qp: u8, qc: u8, tp: u8, tc: u8) -> i32 {
        let qp = usize::from(qp);
        let qc = usize::from(qc);
        let tp = usize::from(tp);
        let tc = usize::from(tc);
        unsafe { *self.stack.get_unchecked(qp * 216 + qc * 36 + tp * 6 + tc) }
    }

    /// M ← Bq transition energy: close a query gap and return to the stack.
    #[inline(always)]
    pub(crate) fn close_query_gap(&self, qp: u8, qc: u8, tc: u8) -> i32 {
        let qp = usize::from(qp);
        let qc = usize::from(qc);
        let tc = usize::from(tc);
        unsafe { *self.close_query_gap.get_unchecked(qp * 36 + qc * 6 + tc) }
    }

    /// M ← Bt transition energy: close a target gap and return to the stack.
    #[inline(always)]
    pub(crate) fn close_target_gap(&self, qc: u8, tp: u8, tc: u8) -> i32 {
        let qc = usize::from(qc);
        let tp = usize::from(tp);
        let tc = usize::from(tc);
        unsafe { *self.close_target_gap.get_unchecked(qc * 36 + tp * 6 + tc) }
    }

    /// Bq ← M transition energy: open a query gap.
    #[inline(always)]
    pub(crate) fn open_query_gap(&self, qp: u8, qc: u8, tc: u8) -> i32 {
        let qp = usize::from(qp);
        let qc = usize::from(qc);
        let tc = usize::from(tc);
        unsafe { *self.open_query_gap.get_unchecked(qp * 36 + qc * 6 + tc) }
    }

    /// Bq ← Bq transition energy: extend a query gap.
    #[inline(always)]
    pub(crate) fn extend_query_gap(&self, qp: u8, qc: u8) -> i32 {
        let qp = usize::from(qp);
        let qc = usize::from(qc);
        unsafe { *self.extend_query_gap.get_unchecked(qp * 6 + qc) }
    }

    /// Bt ← M transition energy: open a target gap.
    #[inline(always)]
    pub(crate) fn open_target_gap(&self, qc: u8, tp: u8, tc: u8) -> i32 {
        let qc = usize::from(qc);
        let tp = usize::from(tp);
        let tc = usize::from(tc);
        unsafe { *self.open_target_gap.get_unchecked(qc * 36 + tp * 6 + tc) }
    }

    /// Bt ← Bt transition energy: extend a target gap.
    #[inline(always)]
    pub(crate) fn extend_target_gap(&self, tp: u8, tc: u8) -> i32 {
        let tp = usize::from(tp);
        let tc = usize::from(tc);
        unsafe { *self.extend_target_gap.get_unchecked(tp * 6 + tc) }
    }

    /// Terminal (boundary) penalty.
    #[inline(always)]
    pub(crate) fn terminal(&self, qc: u8, tc: u8) -> i32 {
        let qc = usize::from(qc);
        let tc = usize::from(tc);
        unsafe { *self.terminal.get_unchecked(qc * 6 + tc) }
    }

    /// Terminal (boundary) penalty for semantic bases.
    #[inline(always)]
    pub(crate) fn terminal_bases(&self, q: Base, t: Base) -> i32 {
        self.terminal(q as u8, t as u8)
    }

    /// 36-entry row for `stack` with fixed `(qp, qc)`: `&[tp*6+tc]`.
    #[inline(always)]
    pub(super) fn stack_row(&self, qp: u8, qc: u8) -> &[i32] {
        let qp = usize::from(qp);
        let qc = usize::from(qc);
        let start = qp * 216 + qc * 36;
        unsafe { self.stack.get_unchecked(start..start + 36) }
    }

    /// 36-entry row for `close_target_gap` with fixed `qc`: `&[tp*6+tc]`.
    #[inline(always)]
    pub(super) fn close_target_gap_row(&self, qc: u8) -> &[i32] {
        let qc = usize::from(qc);
        let start = qc * 36;
        unsafe { self.close_target_gap.get_unchecked(start..start + 36) }
    }

    /// 36-entry row for `open_target_gap` with fixed `qc`: `&[tp*6+tc]`.
    #[inline(always)]
    pub(super) fn open_target_gap_row(&self, qc: u8) -> &[i32] {
        let qc = usize::from(qc);
        let start = qc * 36;
        unsafe { self.open_target_gap.get_unchecked(start..start + 36) }
    }

    /// Full `extend_target_gap` table (36 entries): `&[tp*6+tc]`.
    #[inline(always)]
    pub(super) fn extend_target_gap_row(&self) -> &[i32; 36] {
        &self.extend_target_gap
    }

    /// 6-entry row for `close_query_gap` with fixed `(qp, qc)`: `&[tc]`.
    #[inline(always)]
    pub(super) fn close_query_gap_row(&self, qp: u8, qc: u8) -> &[i32] {
        let qp = usize::from(qp);
        let qc = usize::from(qc);
        let start = qp * 36 + qc * 6;
        unsafe { self.close_query_gap.get_unchecked(start..start + 6) }
    }

    /// 6-entry row for `open_query_gap` with fixed `(qp, qc)`: `&[tc]`.
    #[inline(always)]
    pub(super) fn open_query_gap_row(&self, qp: u8, qc: u8) -> &[i32] {
        let qp = usize::from(qp);
        let qc = usize::from(qc);
        let start = qp * 36 + qc * 6;
        unsafe { self.open_query_gap.get_unchecked(start..start + 6) }
    }

    /// 6-entry row for `terminal` with fixed `qc`: `&[tc]`.
    #[inline(always)]
    pub(super) fn terminal_row(&self, qc: u8) -> &[i32] {
        let qc = usize::from(qc);
        let start = qc * 6;
        unsafe { self.terminal.get_unchecked(start..start + 6) }
    }

    #[cfg_attr(feature = "prof", inline(never))]
    /// Run DP forward pass over `view`, reusing the caller-provided grid.
    pub fn extend(&self, view: &DpView<'_>, grid: &mut DpGrid) -> GotohResult {
        let (q_len, t_len) = (view.q_len.min(MAX_EXT), view.t_len.min(MAX_EXT));

        trace!("{} q_len={} t_len={}", view.dir, q_len, t_len);
        debug_assert!(
            q_len > 0 && t_len > 0,
            "DP view must define a non-empty window"
        );

        // =====================================================================
        // PRECOMPUTE Q/T LOOKUP INDICES (used by init AND main loop)
        // =====================================================================

        let mut q_idx = std::mem::MaybeUninit::<[u8; MAX_EXT]>::uninit();
        let mut t_idx = std::mem::MaybeUninit::<[u8; MAX_EXT]>::uninit();

        let q_ptr = q_idx.as_mut_ptr() as *mut u8;
        let t_ptr = t_idx.as_mut_ptr() as *mut u8;

        match view.dir {
            ExtendDir::Left => {
                for i in 0..q_len {
                    unsafe { *q_ptr.add(i) = (*view.query.get_unchecked(view.q_anchor - i)) as u8 };
                }
                for j in 0..t_len {
                    unsafe { *t_ptr.add(j) = (*view.target.get_unchecked(view.t_anchor + j)) as u8 };
                }
            }
            ExtendDir::Right => {
                for i in 0..q_len {
                    unsafe { *q_ptr.add(i) = (*view.query.get_unchecked(view.q_anchor + i)) as u8 };
                }
                for j in 0..t_len {
                    unsafe { *t_ptr.add(j) = (*view.target.get_unchecked(view.t_anchor - j)) as u8 };
                }
            }
        }

        // Initial score: terminal penalty for seed boundary
        let q0 = unsafe { *q_ptr };
        let t0 = unsafe { *t_ptr };
        let mut best = BestScore::new(self.terminal(q0, t0));

        if q_len <= 1 || t_len <= 1 {
            return GotohResult {
                score: best.score,
                end_i: 0,
                end_j: 0,
            };
        }

        grid.resize(t_len + 1, q_len + 1);

        let has_main_region = self.init_frontier(q_ptr, t_ptr, grid, q_len, t_len, &mut best);
        if !has_main_region {
            return GotohResult {
                score: best.score,
                end_i: best.i,
                end_j: best.j,
            };
        }

        self.dp_main_loop(q_ptr, t_ptr, grid, q_len, t_len, &mut best);

        trace!(
            "{} result: score={} end_i={} end_j={}",
            view.dir,
            best.score,
            best.i,
            best.j,
        );

        GotohResult {
            score: best.score,
            end_i: best.i,
            end_j: best.j,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Matrix;

    #[test]
    fn gotoh_matches_scoring_table() {
        let table = ScoringModel::new(Matrix::T04, 50);
        let gotoh = Gotoh::new(&table);

        for q1 in 0u8..6 {
            for q2 in 0u8..6 {
                for t1 in 0u8..6 {
                    for t2 in 0u8..6 {
                        assert_eq!(
                            gotoh.stack(q1, q2, t1, t2),
                            table.transition_energy(
                                usize::from(q1),
                                usize::from(q2),
                                usize::from(t1),
                                usize::from(t2),
                            ),
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
        let table = ScoringModel::new(Matrix::T04, 50);
        let gotoh = Gotoh::new(&table);

        for qp in 0u8..6 {
            for qc in 0u8..6 {
                for tc in 0u8..6 {
                    assert_eq!(
                        gotoh.close_query_gap(qp, qc, tc),
                        table.transition_energy(usize::from(qp), usize::from(qc), GAP, usize::from(tc)),
                    );
                    assert_eq!(
                        gotoh.open_query_gap(qp, qc, tc),
                        table.transition_energy(usize::from(qp), usize::from(qc), usize::from(tc), GAP),
                    );
                }
                assert_eq!(
                    gotoh.extend_query_gap(qp, qc),
                    table.transition_energy(usize::from(qp), usize::from(qc), GAP, GAP),
                );
            }
        }

        for qc in 0u8..6 {
            for tp in 0u8..6 {
                for tc in 0u8..6 {
                    assert_eq!(
                        gotoh.close_target_gap(qc, tp, tc),
                        table.transition_energy(GAP, usize::from(qc), usize::from(tp), usize::from(tc)),
                    );
                    assert_eq!(
                        gotoh.open_target_gap(qc, tp, tc),
                        table.transition_energy(usize::from(qc), GAP, usize::from(tp), usize::from(tc)),
                    );
                }
            }
            for tc in 0u8..6 {
                assert_eq!(
                    gotoh.terminal(qc, tc),
                    table.transition_energy(usize::from(qc), GAP, usize::from(tc), GAP),
                );
            }
        }

        for tp in 0u8..6 {
            for tc in 0u8..6 {
                assert_eq!(
                    gotoh.extend_target_gap(tp, tc),
                    table.transition_energy(GAP, GAP, usize::from(tp), usize::from(tc)),
                );
            }
        }
    }
}
