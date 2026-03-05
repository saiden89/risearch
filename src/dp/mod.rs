use crate::alignment::PairClass;
use crate::config::{ExtendConfig, ScoreConfig};
use crate::dsm::{DirectionalDsm, DsmModel, GAP};
use crate::types::Base;
use log::trace;
use smallvec::SmallVec;

mod core;
mod init;
mod traceback;
use std::cmp::max;

pub(crate) use self::core::dp_main_loop_generic;
pub(crate) use self::init::init_frontier;
use self::traceback::traceback;

/// Maximum extension length for precomputed index arrays.
/// Matches the typical max_ext parameter (100-200 bases).
const MAX_EXT: usize = 256;

/// DP runtime configuration derived from high-level search configs.
#[derive(Clone, Copy, Debug)]
pub struct DpConfig {
    max_extension: usize,
    penalty_raw: i32,
}

impl DpConfig {
    #[inline(always)]
    pub const fn max_extension(self) -> usize {
        self.max_extension
    }

    #[inline(always)]
    pub const fn penalty_raw(self) -> i32 {
        self.penalty_raw
    }
}

impl From<(&ScoreConfig, &ExtendConfig)> for DpConfig {
    fn from((score, extend): (&ScoreConfig, &ExtendConfig)) -> Self {
        Self {
            max_extension: usize::from(extend.max_extension).min(MAX_EXT),
            penalty_raw: (score.penalty * 100.0).round() as i32,
        }
    }
}

/// Extension direction - determines sequence indexing polarity.
///
/// After canonical DSM construction, direction only affects `q()`/`t()` index
/// arithmetic. All stacking order is resolved by the canonical table.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ExtendDir {
    Left,
    Right,
}

impl std::fmt::Display for ExtendDir {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Left => write!(f, "[EXT_LEFT]"),
            Self::Right => write!(f, "[EXT_RIGHT]"),
        }
    }
}

// =============================================================================
// DP VIEW - Direction-aware sequence access + canonical DSM
// =============================================================================

/// View into sequences for DP extension.
///
/// Stores a `DirectionalDsm` handle that absorbs all direction-dependent stacking
/// order. The DP loop uses uniform calls like `view.e(q_prev, q_curr, t_prev, t_curr)`
/// and the canonical table provides the correct energy for both directions.
pub struct DpView<'a> {
    query: &'a [Base],
    target_transformed: &'a [Base],
    q_anchor: usize,
    t_anchor: usize,
    pub(super) dir: ExtendDir,
    pub q_len: usize,
    pub t_len: usize,
    dsm: DirectionalDsm<'a>,
}

impl<'a> DpView<'a> {
    #[inline(always)]
    fn base_or_gap(seq: &[Base], pos: usize) -> Base {
        if pos >= seq.len() {
            Base::Gap
        } else {
            // SAFETY: bounds checked above
            unsafe { *seq.get_unchecked(pos) }
        }
    }

    #[inline(always)]
    fn left_base(seq: &[Base], anchor: usize, offset: usize) -> Base {
        if offset > anchor {
            Base::Gap
        } else {
            // SAFETY: offset <= anchor, and anchor < seq.len() by construction
            unsafe { *seq.get_unchecked(anchor - offset) }
        }
    }

    #[inline(always)]
    fn right_base(seq: &[Base], anchor: usize, offset: usize) -> Base {
        Self::base_or_gap(seq, anchor + offset)
    }

    /// Create a left extension view (query toward 5', target toward 3')
    pub fn left(
        query: &'a [Base],
        target_transformed: &'a [Base],
        q_start: usize,
        t_start: usize,
        max_ext: usize,
        model: &'a DsmModel,
    ) -> Self {
        Self {
            query,
            target_transformed,
            q_anchor: q_start,
            t_anchor: t_start,
            dir: ExtendDir::Left,
            q_len: (q_start + 1).min(max_ext),
            t_len: (target_transformed.len() - t_start).min(max_ext),
            dsm: model.left(),
        }
    }

    /// Create a right extension view (query toward 3', target toward 5')
    pub fn right(
        query: &'a [Base],
        target_transformed: &'a [Base],
        q_end: usize,
        t_end: usize,
        max_ext: usize,
        model: &'a DsmModel,
    ) -> Self {
        Self {
            query,
            target_transformed,
            q_anchor: q_end,
            t_anchor: t_end,
            dir: ExtendDir::Right,
            q_len: (query.len() - q_end).min(max_ext),
            t_len: (t_end + 1).min(max_ext),
            dsm: model.right(),
        }
    }

    /// Get query base index at DP position i (0 = anchor)
    #[inline(always)]
    pub fn q(&self, i: usize) -> usize {
        match self.dir {
            ExtendDir::Left => Self::left_base(self.query, self.q_anchor, i).idx(),
            ExtendDir::Right => Self::right_base(self.query, self.q_anchor, i).idx(),
        }
    }

    /// Get target base index at DP position j (0 = anchor).
    /// Bases are complemented on-the-fly from the transformed index.
    #[inline(always)]
    pub fn t(&self, j: usize) -> usize {
        match self.dir {
            ExtendDir::Left => Self::right_base(self.target_transformed, self.t_anchor, j)
                .complement()
                .idx(),
            ExtendDir::Right => Self::left_base(self.target_transformed, self.t_anchor, j)
                .complement()
                .idx(),
        }
    }

    /// Canonical stacking energy lookup.
    ///
    /// Arguments are always in (prev, curr, prev, curr) order.
    /// The canonical table handles direction-specific reordering.
    #[inline(always)]
    pub fn e(&self, q1: usize, q2: usize, t1: usize, t2: usize) -> i32 {
        self.dsm.lookup(q1, q2, t1, t2)
    }

    /// Terminal penalty at position (i, j) - stacking with gap boundary.
    #[inline(always)]
    pub fn terminal(&self, i: usize, j: usize) -> i32 {
        self.dsm.terminal(self.q(i), self.t(j))
    }

    /// Match/mismatch energy at (i, j) from diagonal (i-1, j-1)
    #[inline(always)]
    pub fn match_e(&self, i: usize, j: usize) -> i32 {
        self.e(self.q(i - 1), self.q(i), self.t(j - 1), self.t(j))
    }

    // =========================================================================
    // DP TRANSITION HELPERS - exactly match the DP recurrence semantics
    // =========================================================================

    /// M[i,j] from Bq[i-1,j-1]: re-entry to match from query bulge
    #[inline(always)]
    pub fn m_from_bq(&self, i: usize, j: usize) -> i32 {
        self.e(self.q(i - 1), self.q(i), GAP, self.t(j))
    }

    /// M[i,j] from Bt[i-1,j-1]: re-entry to match from target bulge
    #[inline(always)]
    pub fn m_from_bt(&self, i: usize, j: usize) -> i32 {
        self.e(GAP, self.q(i), self.t(j - 1), self.t(j))
    }

    /// Bq[i,j] from M[i-1,j]: open query bulge (gap in target)
    #[inline(always)]
    pub fn bq_open(&self, i: usize, j: usize) -> i32 {
        self.e(self.q(i - 1), self.q(i), self.t(j), GAP)
    }

    /// Bq[i,j] from Bq[i-1,j]: extend query bulge
    #[inline(always)]
    pub fn bq_ext(&self, i: usize) -> i32 {
        self.e(self.q(i - 1), self.q(i), GAP, GAP)
    }

    /// Bt[i,j] from M[i,j-1]: open target bulge (gap in query)
    #[inline(always)]
    pub fn bt_open(&self, i: usize, j: usize) -> i32 {
        self.e(self.q(i), GAP, self.t(j - 1), self.t(j))
    }

    /// Bt[i,j] from Bt[i,j-1]: extend target bulge
    #[inline(always)]
    pub fn bt_ext(&self, j: usize) -> i32 {
        self.e(GAP, GAP, self.t(j - 1), self.t(j))
    }

    /// Canonical DSM accessor for core loop and init code.
    #[inline(always)]
    pub(super) fn dsm(&self) -> &DirectionalDsm<'a> {
        &self.dsm
    }
}

/// Negative infinity for the (max, +) semiring over DP scores.
///
/// Must satisfy two invariants (enforced by compile-time assert below):
/// 1. Invalid scores can never drift into valid range through accumulated adds
/// 2. No i32 underflow from accumulated negative energy
pub(super) const NEG_INF: i32 = -1_000_000_000;

/// Conservative upper bound on |energy| from a single DSM lookup.
/// Source tables are i16 (max 32767); penalty adds modest overhead.
/// Real values are ~300-400 (0.01 kcal/mol units), but we bound generously.
/// Enforced at runtime in DsmModel::new.
const MAX_ENERGY: i64 = 40_000;

// Compile-time proof that NEG_INF arithmetic is safe for MAX_EXT.
const _: () = {
    // Longest path through MAX_EXT × MAX_EXT grid
    let max_drift = 2 * MAX_EXT as i64 * MAX_ENERGY;
    let neg_inf_abs = -(NEG_INF as i64);

    // Invalid scores must stay below valid range after max positive drift
    assert!(
        neg_inf_abs > 2 * max_drift,
        "NEG_INF too close to zero: invalid scores could enter valid range"
    );
    // NEG_INF minus max negative drift must not underflow i32
    assert!(
        neg_inf_abs + max_drift < (i32::MAX as i64 + 1),
        "NEG_INF too close to i32::MIN: arithmetic could overflow"
    );
};

/// Tracks the best scoring position found during DP extension.
#[derive(Clone, Copy)]
pub(super) struct BestScore {
    pub(super) score: i32,
    pub(super) i: usize,
    pub(super) j: usize,
}

impl BestScore {
    pub(super) fn new(score: i32) -> Self {
        Self { score, i: 0, j: 0 }
    }

    /// Update if `val + term` exceeds current best.
    #[inline(always)]
    pub(super) fn update(&mut self, val: i32, term: i32, i: usize, j: usize) {
        let curr = val + term;
        if curr > self.score {
            self.score = curr;
            self.i = i;
            self.j = j;
        }
    }
}

/// Add stacking energy to a DP score.
///
/// NEG_INF propagates naturally: NEG_INF + energy ≈ NEG_INF,
/// which can never reach the valid score range (see NEG_INF docs).
#[inline(always)]
pub(super) fn add_e(base: i32, energy: i32) -> i32 {
    base + energy
}

/// max of 3 values - branchless
#[inline(always)]
pub(super) fn max3(a: i32, b: i32, c: i32) -> i32 {
    max(max(a, b), c)
}

// =============================================================================
// DP GRID - Interleaved AoS layout for cache-friendly cell access
// =============================================================================

/// Single DP cell: all three state scores packed together.
///
/// `repr(C)` guarantees field order and no padding (3 × i32 = 12 bytes).
/// Every access in the DP touches multiple states at the same (i,j),
/// so interleaving them maximizes cache line utilization.
#[repr(C)]
#[derive(Clone, Copy)]
pub(super) struct DpCell {
    pub(super) m: i32,  // Match/mismatch state
    pub(super) bq: i32, // Query bulge (gap in target)
    pub(super) bt: i32, // Target bulge (gap in query)
}

impl DpCell {
    pub(super) const EMPTY: Self = Self {
        m: NEG_INF,
        bq: NEG_INF,
        bt: NEG_INF,
    };
}

/// Interleaved DP grid storing `DpCell` per position.
///
/// Row-major layout: cell (i, j) is at index `i * width + j`.
/// Reused across extensions (resized, not reallocated).
pub(super) struct DpGrid {
    data: Vec<DpCell>,
    width: usize,
}

impl DpGrid {
    fn new(width: usize, height: usize) -> Self {
        Self {
            data: vec![DpCell::EMPTY; width * height],
            width,
        }
    }

    /// Resize without clearing. Caller must ensure all accessed cells are initialized.
    /// This is critical for performance (avoids O(n²) overhead per extension).
    #[inline]
    pub(super) fn resize(&mut self, width: usize, height: usize) {
        let new_len = width * height;
        if self.data.len() < new_len {
            self.data.resize(new_len, DpCell::EMPTY);
        }
        self.width = width;
    }

    #[inline(always)]
    pub(super) fn get(&self, i: usize, j: usize) -> DpCell {
        let idx = i * self.width + j;
        debug_assert!(
            idx < self.data.len(),
            "DpGrid::get out of bounds: ({}, {}) idx={} len={}",
            i,
            j,
            idx,
            self.data.len()
        );
        unsafe { *self.data.get_unchecked(idx) }
    }

    #[inline(always)]
    pub(super) fn ptr(&mut self) -> *mut DpCell {
        self.data.as_mut_ptr()
    }

    #[inline(always)]
    pub(super) fn width(&self) -> usize {
        self.width
    }
}

/// Stateful DP extender with reusable grid.
///
/// Direction-agnostic: the canonical DSM in `DpView` resolves all
/// direction-dependent stacking order at view construction time.
pub struct DpExtender {
    grid: DpGrid,
}

/// Guard holding a reference to the DP grid after forward pass.
/// While this exists, the extender cannot be used for another extension
/// (borrow checker enforces this via the lifetime on `grid`).
pub struct ExtendResult<'a> {
    grid: &'a DpGrid,
    pub score: i32,
    pub q_len: usize,
    pub t_len: usize,
}

impl<'a> ExtendResult<'a> {
    /// Run traceback to reconstruct alignment as Pairings.
    /// Only call when alignment output is needed (skip for Minimal format).
    pub fn traceback(&self, view: &DpView<'_>) -> SmallVec<[PairClass; 64]> {
        let mut out = SmallVec::new();
        traceback(view, self.grid, self.q_len, self.t_len, &mut out);
        out
    }
}

impl DpExtender {
    pub fn new(max_extension: usize) -> Self {
        let side = max_extension.min(MAX_EXT).saturating_add(1).max(1);
        Self {
            grid: DpGrid::new(side, side),
        }
    }

    #[cfg_attr(feature = "prof", inline(never))]
    /// Run DP forward pass and return an ExtendResult guard.
    /// Call `.traceback()` on the result if alignment is needed.
    pub fn extend(&mut self, view: &DpView<'_>) -> ExtendResult<'_> {
        let (q_len, t_len) = (view.q_len.min(MAX_EXT), view.t_len.min(MAX_EXT));

        trace!("{} q_len={} t_len={}", view.dir, q_len, t_len);

        // Initial score: terminal penalty for seed boundary
        let mut best = BestScore::new(view.terminal(0, 0));

        // Early return
        if q_len <= 1 || t_len <= 1 {
            return ExtendResult {
                grid: &self.grid,
                score: best.score,
                q_len: 0,
                t_len: 0,
            };
        }

        // Resize grid
        self.grid.resize(t_len + 1, q_len + 1);

        // =====================================================================
        // PRECOMPUTE Q/T BASE INDICES (used by init AND main loop)
        // =====================================================================

        let mut q_idx = std::mem::MaybeUninit::<[usize; MAX_EXT]>::uninit();
        let mut t_idx = std::mem::MaybeUninit::<[usize; MAX_EXT]>::uninit();

        // SAFETY: We only read indices we explicitly write below.
        let q_ptr = q_idx.as_mut_ptr() as *mut usize;
        let t_ptr = t_idx.as_mut_ptr() as *mut usize;

        for i in 0..q_len {
            unsafe { *q_ptr.add(i) = view.q(i) };
        }
        for j in 0..t_len {
            unsafe { *t_ptr.add(j) = view.t(j) };
        }

        // =====================================================================
        // INITIALIZATION - Unconditional writes to avoid stale data
        // =====================================================================

        let has_main_region = init_frontier(
            q_ptr,
            t_ptr,
            &mut self.grid,
            q_len,
            t_len,
            view.dsm(),
            &mut best,
        );
        if !has_main_region {
            return ExtendResult {
                grid: &self.grid,
                score: best.score,
                q_len: best.i,
                t_len: best.j,
            };
        }

        // =======================================================================
        // MAIN DP LOOP (i >= 3, j >= 3)
        // =======================================================================

        dp_main_loop_generic(
            q_ptr,
            t_ptr,
            &mut self.grid,
            q_len,
            t_len,
            view.dsm(),
            &mut best,
        );

        trace!(
            "{} result: score={} q_len={} t_len={}",
            view.dir,
            best.score,
            best.i,
            best.j,
        );

        ExtendResult {
            grid: &self.grid,
            score: best.score,
            q_len: best.i,
            t_len: best.j,
        }
    }
}
