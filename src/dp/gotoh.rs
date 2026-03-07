//! Gotoh 3-state DP recurrence: precomputed transition tables.
//!
//! All 7 transitions are materialized at construction time from a `ScoringTable`.
//! The inner DP loop uses zero-cost slice accesses into these tables instead of
//! computing GAP-pattern lookups on the fly.
//!
//! # Index safety
//!
//! All index arguments to point-lookup and slice-accessor methods must be in
//! `0..6` (valid `Base::idx()` values). Callers source indices from
//! `DpView::q()`/`DpView::t()`, which return `Base::idx()` on a `#[repr(u8)]`
//! enum with variants 0–5.

use crate::dsm::{ScoringTable, GAP};
use crate::types::BASE_COUNT;
use log::trace;

use super::{BestScore, DpGrid, DpView, ExtendResult, MAX_EXT};

const BC: usize = BASE_COUNT; // 6

/// Precomputed transition tables for the Gotoh 3-state DP.
///
/// Direction-agnostic: each instance corresponds to one `ScoringTable` orientation.
pub struct Gotoh {
    /// M ← M: `match_mm[qp][qc][tp][tc]` — full 4D stacking (1296 entries).
    ///
    /// Sized as 6⁴ rather than 5⁴ so that all 36-entry profile slices share the
    /// uniform stride `tp*6 + tc`. This lets the inner loop compute one index
    /// (`t_stack_idx`) and reuse it across match_mm, m_from_bt, bt_open, and
    /// bt_extend lookups. GAP-indexed entries (where any index is 0) are never
    /// accessed in the M←M path — they exist solely to preserve stride uniformity.
    match_mm: [i32; 1296],
    /// M ← Bq: `m_from_bq[qp][qc][tc]` — re-entry from query bulge (216 entries)
    m_from_bq: [i32; 216],
    /// M ← Bt: `m_from_bt[qc][tp][tc]` — re-entry from target bulge (216 entries)
    m_from_bt: [i32; 216],
    /// Bq ← M: `bq_open[qp][qc][tc]` — open query bulge (216 entries)
    bq_open: [i32; 216],
    /// Bq ← Bq: `bq_extend[qp][qc]` — extend query bulge (36 entries)
    bq_extend: [i32; 36],
    /// Bt ← M: `bt_open[qc][tp][tc]` — open target bulge (216 entries)
    bt_open: [i32; 216],
    /// Bt ← Bt: `bt_extend[tp][tc]` — extend target bulge (36 entries)
    bt_extend: [i32; 36],
    /// Terminal: `terminal[qc][tc]` — boundary penalty (36 entries)
    terminal: [i32; 36],
}

impl Gotoh {
    /// Materialize all transition tables from a `ScoringTable`.
    pub fn new(table: &ScoringTable) -> Self {
        let mut match_mm = [0i32; 1296];
        let mut m_from_bq = [0i32; 216];
        let mut m_from_bt = [0i32; 216];
        let mut bq_open = [0i32; 216];
        let mut bq_extend = [0i32; 36];
        let mut bt_open = [0i32; 216];
        let mut bt_extend = [0i32; 36];
        let mut terminal = [0i32; 36];

        for qp in 0..BC {
            for qc in 0..BC {
                let qp_qc_36 = qp * BC * BC * BC + qc * BC * BC; // offset into match_mm
                let qp_qc_6 = qp * BC + qc; // offset into bq_extend
                for tp in 0..BC {
                    for tc in 0..BC {
                        match_mm[qp_qc_36 + tp * BC + tc] = table.lookup(qp, qc, tp, tc);
                    }
                }
                for tc in 0..BC {
                    m_from_bq[qp * BC * BC + qc * BC + tc] = table.lookup(qp, qc, GAP, tc);
                    bq_open[qp * BC * BC + qc * BC + tc] = table.lookup(qp, qc, tc, GAP);
                }
                bq_extend[qp_qc_6] = table.lookup(qp, qc, GAP, GAP);
            }
        }

        for qc in 0..BC {
            for tp in 0..BC {
                for tc in 0..BC {
                    m_from_bt[qc * BC * BC + tp * BC + tc] = table.lookup(GAP, qc, tp, tc);
                    bt_open[qc * BC * BC + tp * BC + tc] = table.lookup(qc, GAP, tp, tc);
                }
            }
            for tc in 0..BC {
                terminal[qc * BC + tc] = table.lookup(qc, GAP, tc, GAP);
            }
        }

        for tp in 0..BC {
            for tc in 0..BC {
                bt_extend[tp * BC + tc] = table.lookup(GAP, GAP, tp, tc);
            }
        }

        Self {
            match_mm,
            m_from_bq,
            m_from_bt,
            bq_open,
            bq_extend,
            bt_open,
            bt_extend,
            terminal,
        }
    }

    // =========================================================================
    // POINT LOOKUPS — for init and traceback
    // =========================================================================

    /// M ← M transition energy.
    #[inline(always)]
    pub fn match_energy(&self, qp: usize, qc: usize, tp: usize, tc: usize) -> i32 {
        // SAFETY: All args are Base::idx() ∈ 0..6.
        // Max index = 5*216 + 5*36 + 5*6 + 5 = 1295 < 1296.
        unsafe {
            *self
                .match_mm
                .get_unchecked(qp * 216 + qc * 36 + tp * 6 + tc)
        }
    }

    /// M ← Bq transition energy.
    #[inline(always)]
    pub fn m_from_bq(&self, qp: usize, qc: usize, tc: usize) -> i32 {
        // SAFETY: All args ∈ 0..6. Max index = 5*36 + 5*6 + 5 = 215 < 216.
        unsafe { *self.m_from_bq.get_unchecked(qp * 36 + qc * 6 + tc) }
    }

    /// M ← Bt transition energy.
    #[inline(always)]
    pub fn m_from_bt(&self, qc: usize, tp: usize, tc: usize) -> i32 {
        // SAFETY: All args ∈ 0..6. Max index = 5*36 + 5*6 + 5 = 215 < 216.
        unsafe { *self.m_from_bt.get_unchecked(qc * 36 + tp * 6 + tc) }
    }

    /// Bq ← M (open query bulge) transition energy.
    #[inline(always)]
    pub fn bq_open(&self, qp: usize, qc: usize, tc: usize) -> i32 {
        // SAFETY: All args ∈ 0..6. Max index = 5*36 + 5*6 + 5 = 215 < 216.
        unsafe { *self.bq_open.get_unchecked(qp * 36 + qc * 6 + tc) }
    }

    /// Bq ← Bq (extend query bulge) transition energy.
    #[inline(always)]
    pub fn bq_extend(&self, qp: usize, qc: usize) -> i32 {
        // SAFETY: All args ∈ 0..6. Max index = 5*6 + 5 = 35 < 36.
        unsafe { *self.bq_extend.get_unchecked(qp * 6 + qc) }
    }

    /// Bt ← M (open target bulge) transition energy.
    #[inline(always)]
    pub fn bt_open(&self, qc: usize, tp: usize, tc: usize) -> i32 {
        // SAFETY: All args ∈ 0..6. Max index = 5*36 + 5*6 + 5 = 215 < 216.
        unsafe { *self.bt_open.get_unchecked(qc * 36 + tp * 6 + tc) }
    }

    /// Bt ← Bt (extend target bulge) transition energy.
    #[inline(always)]
    pub fn bt_extend_e(&self, tp: usize, tc: usize) -> i32 {
        // SAFETY: All args ∈ 0..6. Max index = 5*6 + 5 = 35 < 36.
        unsafe { *self.bt_extend.get_unchecked(tp * 6 + tc) }
    }

    /// Terminal (boundary) penalty.
    #[inline(always)]
    pub fn terminal(&self, qc: usize, tc: usize) -> i32 {
        // SAFETY: All args ∈ 0..6. Max index = 5*6 + 5 = 35 < 36.
        unsafe { *self.terminal.get_unchecked(qc * 6 + tc) }
    }

    // =========================================================================
    // ROW-PROFILE SLICE ACCESSORS — for core loop
    // =========================================================================

    /// 36-entry slice of match_mm for fixed (qp, qc): `&[tp*6+tc]`.
    #[inline(always)]
    pub fn match_profile(&self, qp: usize, qc: usize) -> &[i32] {
        let start = qp * 216 + qc * 36;
        // SAFETY: qp, qc ∈ 0..6. Max start = 5*216 + 5*36 = 1260, end = 1296 = len.
        unsafe { self.match_mm.get_unchecked(start..start + 36) }
    }

    /// 36-entry slice of m_from_bt for fixed qc: `&[tp*6+tc]`.
    #[inline(always)]
    pub fn m_from_bt_profile(&self, qc: usize) -> &[i32] {
        let start = qc * 36;
        // SAFETY: qc ∈ 0..6. Max start = 180, end = 216 = len.
        unsafe { self.m_from_bt.get_unchecked(start..start + 36) }
    }

    /// 36-entry slice of bt_open for fixed qc: `&[tp*6+tc]`.
    #[inline(always)]
    pub fn bt_open_profile(&self, qc: usize) -> &[i32] {
        let start = qc * 36;
        // SAFETY: qc ∈ 0..6. Max start = 180, end = 216 = len.
        unsafe { self.bt_open.get_unchecked(start..start + 36) }
    }

    /// Full bt_extend table (36 entries): `&[tp*6+tc]`.
    #[inline(always)]
    pub fn bt_extend_profile(&self) -> &[i32; 36] {
        &self.bt_extend
    }

    /// 6-entry slice of m_from_bq for fixed (qp, qc): `&[tc]`.
    #[inline(always)]
    pub fn m_from_bq_profile(&self, qp: usize, qc: usize) -> &[i32] {
        let start = qp * 36 + qc * 6;
        // SAFETY: qp, qc ∈ 0..6. Max start = 5*36 + 5*6 = 210, end = 216 = len.
        unsafe { self.m_from_bq.get_unchecked(start..start + 6) }
    }

    /// 6-entry slice of bq_open for fixed (qp, qc): `&[tc]`.
    #[inline(always)]
    pub fn bq_open_profile(&self, qp: usize, qc: usize) -> &[i32] {
        let start = qp * 36 + qc * 6;
        // SAFETY: qp, qc ∈ 0..6. Max start = 5*36 + 5*6 = 210, end = 216 = len.
        unsafe { self.bq_open.get_unchecked(start..start + 6) }
    }

    /// 6-entry slice of terminal for fixed qc: `&[tc]`.
    #[inline(always)]
    pub fn terminal_profile(&self, qc: usize) -> &[i32] {
        let start = qc * 6;
        // SAFETY: qc ∈ 0..6. Max start = 30, end = 36 = len.
        unsafe { self.terminal.get_unchecked(start..start + 6) }
    }

    #[cfg_attr(feature = "prof", inline(never))]
    /// Run DP forward pass over `view`, reusing the caller-provided grid.
    /// Call `.traceback()` on the result if alignment is needed.
    pub fn extend<'a>(&'a self, view: &DpView<'_>, grid: &'a mut DpGrid) -> ExtendResult<'a> {
        let (q_len, t_len) = (view.q_len.min(MAX_EXT), view.t_len.min(MAX_EXT));

        trace!("{} q_len={} t_len={}", view.dir, q_len, t_len);

        // Initial score: terminal penalty for seed boundary
        let mut best = BestScore::new(self.terminal(view.q(0), view.t(0)));

        // Early return
        if q_len <= 1 || t_len <= 1 {
            return ExtendResult {
                grid,
                gotoh: self,
                score: best.score,
                q_len: 0,
                t_len: 0,
            };
        }

        // Resize grid
        grid.resize(t_len + 1, q_len + 1);

        // =====================================================================
        // PRECOMPUTE Q/T BASE INDICES (used by init AND main loop)
        // =====================================================================

        let mut q_idx = std::mem::MaybeUninit::<[usize; MAX_EXT]>::uninit();
        let mut t_idx = std::mem::MaybeUninit::<[usize; MAX_EXT]>::uninit();

        let q_ptr = q_idx.as_mut_ptr() as *mut usize;
        let t_ptr = t_idx.as_mut_ptr() as *mut usize;

        for i in 0..q_len {
            // SAFETY: i < q_len ≤ MAX_EXT, so q_ptr.add(i) is within the
            // MaybeUninit allocation. view.q(i) returns Base::idx() ∈ 0..6.
            unsafe { *q_ptr.add(i) = view.q(i) };
        }
        for j in 0..t_len {
            // SAFETY: j < t_len ≤ MAX_EXT, so t_ptr.add(j) is within the
            // MaybeUninit allocation. view.t(j) returns Base::idx() ∈ 0..6.
            unsafe { *t_ptr.add(j) = view.t(j) };
        }

        // =====================================================================
        // INITIALIZATION - Unconditional writes to avoid stale data
        // =====================================================================

        let has_main_region = self.init_frontier(q_ptr, t_ptr, grid, q_len, t_len, &mut best);
        if !has_main_region {
            return ExtendResult {
                grid,
                gotoh: self,
                score: best.score,
                q_len: best.i,
                t_len: best.j,
            };
        }

        // =======================================================================
        // MAIN DP LOOP (i >= 3, j >= 3)
        // =======================================================================

        self.dp_main_loop_generic(q_ptr, t_ptr, grid, q_len, t_len, &mut best);

        trace!(
            "{} result: score={} q_len={} t_len={}",
            view.dir,
            best.score,
            best.i,
            best.j,
        );

        ExtendResult {
            grid,
            gotoh: self,
            score: best.score,
            q_len: best.i,
            t_len: best.j,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Matrix;

    #[test]
    fn gotoh_matches_scoring_table() {
        let table = ScoringTable::new(Matrix::T04, 50, true);
        let gotoh = Gotoh::new(&table);

        for q1 in 0..6 {
            for q2 in 0..6 {
                for t1 in 0..6 {
                    for t2 in 0..6 {
                        assert_eq!(
                            gotoh.match_energy(q1, q2, t1, t2),
                            table.lookup(q1, q2, t1, t2),
                            "match_energy mismatch at ({},{},{},{})",
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
    fn gotoh_transitions_match_gap_patterns() {
        let table = ScoringTable::new(Matrix::T04, 50, true);
        let gotoh = Gotoh::new(&table);

        for qp in 0..6 {
            for qc in 0..6 {
                for tc in 0..6 {
                    assert_eq!(gotoh.m_from_bq(qp, qc, tc), table.lookup(qp, qc, GAP, tc),);
                    assert_eq!(gotoh.bq_open(qp, qc, tc), table.lookup(qp, qc, tc, GAP),);
                }
                assert_eq!(gotoh.bq_extend(qp, qc), table.lookup(qp, qc, GAP, GAP),);
            }
        }

        for qc in 0..6 {
            for tp in 0..6 {
                for tc in 0..6 {
                    assert_eq!(gotoh.m_from_bt(qc, tp, tc), table.lookup(GAP, qc, tp, tc),);
                    assert_eq!(gotoh.bt_open(qc, tp, tc), table.lookup(qc, GAP, tp, tc),);
                }
            }
            for tc in 0..6 {
                assert_eq!(gotoh.terminal(qc, tc), table.lookup(qc, GAP, tc, GAP),);
            }
        }

        for tp in 0..6 {
            for tc in 0..6 {
                assert_eq!(gotoh.bt_extend_e(tp, tc), table.lookup(GAP, GAP, tp, tc),);
            }
        }
    }
}
