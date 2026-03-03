use crate::alignment::PairClass;
use crate::config::{ExtendConfig, ScoreConfig};
use crate::dsm::{build_penalty_adjusted_flat, dsm_flat_idx, DsmModel, DSM_FLAT_SIZE};
use crate::types::Base;
use log::trace;
use smallvec::SmallVec;

mod core;
mod init;
mod traceback;

use self::core::dp_main_loop_generic;
use self::init::init_frontier;
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

/// Extension direction - determines terminal stacking order
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ExtendDir {
    Left,  // Terminal: Gap→Q, Gap→T (extending into sequence from gap)
    Right, // Terminal: Q→Gap, T→Gap (extending out of sequence into gap)
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
// DP VIEW - Direction-aware sequence access abstraction
// =============================================================================

/// View into sequences for DP extension.
///
/// Abstracts away index arithmetic and DSM stacking order differences
/// between left and right extensions. The DP loop uses uniform calls
/// like `view.e(q_prev, q_curr, t_prev, t_curr)` and the view handles
/// direction-specific reordering internally.
///
/// # Stacking Order
///
/// RNA stacking energy is directional (5'→3'). For nearest-neighbor model:
/// - **Left extension** (toward 5'): DSM[curr, prev, curr, prev]
/// - **Right extension** (toward 3'): DSM[prev, curr, prev, curr]
///
/// The `e()` method always takes arguments in (prev, curr, prev, curr) order
/// and internally reorders for left extension.
pub struct DpView<'a, M: DsmModel> {
    query: &'a [Base],
    target_transformed: &'a [Base],
    q_anchor: usize,
    t_anchor: usize,
    pub dir: ExtendDir,
    pub q_len: usize,
    pub t_len: usize,
    _model: std::marker::PhantomData<M>,
}

/// Gap base index constant
const GAP: usize = Base::Gap as usize;

impl<'a, M: DsmModel> DpView<'a, M> {
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
    ) -> Self {
        Self {
            query,
            target_transformed,
            q_anchor: q_start,
            t_anchor: t_start,
            dir: ExtendDir::Left,
            q_len: (q_start + 1).min(max_ext),
            t_len: (target_transformed.len() - t_start).min(max_ext),
            _model: std::marker::PhantomData,
        }
    }

    /// Create a right extension view (query toward 3', target toward 5')
    pub fn right(
        query: &'a [Base],
        target_transformed: &'a [Base],
        q_end: usize,
        t_end: usize,
        max_ext: usize,
    ) -> Self {
        Self {
            query,
            target_transformed,
            q_anchor: q_end,
            t_anchor: t_end,
            dir: ExtendDir::Right,
            q_len: (query.len() - q_end).min(max_ext),
            t_len: (t_end + 1).min(max_ext),
            _model: std::marker::PhantomData,
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

    /// Stacking energy with direction-aware ordering.
    ///
    /// Always called with arguments in (prev, curr, prev, curr) order.
    /// Left extension internally swaps to (curr, prev, curr, prev).
    #[inline(always)]
    pub fn e(&self, q1: usize, q2: usize, t1: usize, t2: usize, dsm: &[i32; DSM_FLAT_SIZE]) -> i32 {
        match self.dir {
            // Left: DSM[curr, prev, curr, prev] - stacking toward 5'
            ExtendDir::Left => dsm[dsm_flat_idx(q2, q1, t2, t1)],
            // Right: DSM[prev, curr, prev, curr] - stacking toward 3'
            ExtendDir::Right => dsm[dsm_flat_idx(q1, q2, t1, t2)],
        }
    }

    /// Terminal penalty at position (i, j) - stacking with gap boundary.
    ///
    /// Terminal energy is direction-specific (cannot use view.e() swap):
    /// - LEFT: GAP at 5' side, base at 3' → DSM[GAP, q, GAP, t]
    /// - RIGHT: base at 5' side, GAP at 3' → DSM[q, GAP, t, GAP]
    #[inline(always)]
    pub fn terminal(&self, i: usize, j: usize, dsm: &[i32; DSM_FLAT_SIZE]) -> i32 {
        match self.dir {
            ExtendDir::Left => dsm[dsm_flat_idx(GAP, self.q(i), GAP, self.t(j))],
            ExtendDir::Right => dsm[dsm_flat_idx(self.q(i), GAP, self.t(j), GAP)],
        }
    }

    /// Match/mismatch energy at (i, j) from diagonal (i-1, j-1)
    #[inline(always)]
    pub fn match_e(&self, i: usize, j: usize, dsm: &[i32; DSM_FLAT_SIZE]) -> i32 {
        self.e(self.q(i - 1), self.q(i), self.t(j - 1), self.t(j), dsm)
    }

    // =========================================================================
    // DP TRANSITION HELPERS - exactly match the DP recurrence semantics
    // =========================================================================
    // All methods take arguments in (prev, curr) order; view.e() swaps for LEFT.
    // These match the actual energy lookups in the DP main loop.

    /// M[i,j] from Bq[i-1,j-1]: re-entry to match from query bulge
    #[inline(always)]
    pub fn m_from_bq(&self, i: usize, j: usize, dsm: &[i32; DSM_FLAT_SIZE]) -> i32 {
        self.e(self.q(i - 1), self.q(i), GAP, self.t(j), dsm)
    }

    /// M[i,j] from Bt[i-1,j-1]: re-entry to match from target bulge
    #[inline(always)]
    pub fn m_from_bt(&self, i: usize, j: usize, dsm: &[i32; DSM_FLAT_SIZE]) -> i32 {
        self.e(GAP, self.q(i), self.t(j - 1), self.t(j), dsm)
    }

    /// Bq[i,j] from M[i-1,j]: open query bulge (gap in target)
    #[inline(always)]
    pub fn bq_open(&self, i: usize, j: usize, dsm: &[i32; DSM_FLAT_SIZE]) -> i32 {
        self.e(self.q(i - 1), self.q(i), self.t(j), GAP, dsm)
    }

    /// Bq[i,j] from Bq[i-1,j]: extend query bulge
    #[inline(always)]
    pub fn bq_ext(&self, i: usize, dsm: &[i32; DSM_FLAT_SIZE]) -> i32 {
        self.e(self.q(i - 1), self.q(i), GAP, GAP, dsm)
    }

    /// Bt[i,j] from M[i,j-1]: open target bulge (gap in query)
    #[inline(always)]
    pub fn bt_open(&self, i: usize, j: usize, dsm: &[i32; DSM_FLAT_SIZE]) -> i32 {
        self.e(self.q(i), GAP, self.t(j - 1), self.t(j), dsm)
    }

    /// Bt[i,j] from Bt[i,j-1]: extend target bulge
    #[inline(always)]
    pub fn bt_ext(&self, j: usize, dsm: &[i32; DSM_FLAT_SIZE]) -> i32 {
        self.e(GAP, GAP, self.t(j - 1), self.t(j), dsm)
    }
}

/// Sentinel value for invalid/uninitialized score (~ -1 billion).
/// Safe for arithmetic: MIN_SCORE + energy (±2000) will not overflow/underflow i32.
const MIN_SCORE: i32 = -1_000_000_000;

/// Tracks the best scoring position found during DP extension.
#[derive(Clone, Copy)]
struct BestScore {
    score: i32,
    i: usize,
    j: usize,
}

impl BestScore {
    fn new(score: i32) -> Self {
        Self { score, i: 0, j: 0 }
    }

    /// Update if `val + term` exceeds current best.
    #[inline(always)]
    fn update(&mut self, val: i32, term: i32, i: usize, j: usize) {
        if val > MIN_SCORE {
            let curr = val + term;
            if curr > self.score {
                self.score = curr;
                self.i = i;
                self.j = j;
            }
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
    use std::cmp::max;
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
        m: MIN_SCORE,
        bq: MIN_SCORE,
        bt: MIN_SCORE,
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

/// Stateful DP extender with reusable grid
pub struct DpExtender<M: DsmModel> {
    grid: DpGrid,
    dsm_adjusted: [i32; DSM_FLAT_SIZE],
    gap_gap_profile: [i32; 36],
    _model: std::marker::PhantomData<M>,
}

/// Guard holding a reference to the DP grid after forward pass.
/// While this exists, the extender cannot be used for another extension
/// (borrow checker enforces this via the lifetime on `grid`).
pub struct ExtendResult<'a, M: DsmModel> {
    grid: &'a DpGrid,
    dsm_adjusted: &'a [i32; DSM_FLAT_SIZE],
    pub score: i32,
    pub q_len: usize,
    pub t_len: usize,
    _model: std::marker::PhantomData<M>,
}

impl<M: DsmModel> ExtendResult<'_, M> {
    /// Run traceback to reconstruct alignment as Pairings.
    /// Only call when alignment output is needed (skip for Minimal format).
    pub fn traceback(&self, view: &DpView<'_, M>) -> SmallVec<[PairClass; 64]> {
        let mut out = SmallVec::new();
        traceback::<M>(
            view,
            self.grid,
            self.dsm_adjusted,
            self.q_len,
            self.t_len,
            &mut out,
        );
        out
    }
}

impl<M: DsmModel> DpExtender<M> {
    pub fn new() -> Self {
        Self::with_penalty(0)
    }

    pub fn from_config(cfg: DpConfig) -> Self {
        Self::with_penalty(cfg.penalty_raw)
    }

    pub fn with_penalty(penalty: i32) -> Self {
        let dsm_adjusted = build_penalty_adjusted_flat::<M>(penalty);
        let mut gap_gap_profile = [0i32; 36];
        for t1 in 0..6 {
            for t2 in 0..6 {
                gap_gap_profile[t1 * 6 + t2] = dsm_adjusted[dsm_flat_idx(GAP, GAP, t1, t2)];
            }
        }

        Self {
            grid: DpGrid::new(200, 200),
            dsm_adjusted,
            gap_gap_profile,
            _model: std::marker::PhantomData,
        }
    }

    #[cfg_attr(feature = "prof", inline(never))]
    /// Run DP forward pass and return an ExtendResult guard.
    /// Call `.traceback()` on the result if alignment is needed.
    pub fn extend(&mut self, view: &DpView<'_, M>) -> ExtendResult<'_, M> {
        let (q_len, t_len) = (view.q_len.min(MAX_EXT), view.t_len.min(MAX_EXT));

        trace!("{} q_len={} t_len={}", view.dir, q_len, t_len);

        // Initial score: terminal penalty for seed boundary
        let mut best = BestScore::new(view.terminal(0, 0, &self.dsm_adjusted));

        // Early return
        if q_len <= 1 || t_len <= 1 {
            return ExtendResult {
                grid: &self.grid,
                dsm_adjusted: &self.dsm_adjusted,
                score: best.score,
                q_len: 0,
                t_len: 0,
                _model: std::marker::PhantomData,
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

        let has_main_region = init_frontier::<M>(
            view,
            q_ptr,
            t_ptr,
            &mut self.grid,
            q_len,
            t_len,
            &self.dsm_adjusted,
            &mut best,
        );
        if !has_main_region {
            return ExtendResult {
                grid: &self.grid,
                dsm_adjusted: &self.dsm_adjusted,
                score: best.score,
                q_len: best.i,
                t_len: best.j,
                _model: std::marker::PhantomData,
            };
        }

        // =======================================================================
        // MAIN DP LOOP (i >= 3, j >= 3)
        // =======================================================================

        if view.dir == ExtendDir::Left {
            dp_main_loop_generic::<true>(
                q_ptr,
                t_ptr,
                &mut self.grid,
                q_len,
                t_len,
                &self.dsm_adjusted,
                &self.gap_gap_profile,
                &mut best,
            );
        } else {
            dp_main_loop_generic::<false>(
                q_ptr,
                t_ptr,
                &mut self.grid,
                q_len,
                t_len,
                &self.dsm_adjusted,
                &self.gap_gap_profile,
                &mut best,
            );
        }

        trace!(
            "{} result: score={} q_len={} t_len={}",
            view.dir,
            best.score,
            best.i,
            best.j,
        );

        ExtendResult {
            grid: &self.grid,
            dsm_adjusted: &self.dsm_adjusted,
            score: best.score,
            q_len: best.i,
            t_len: best.j,
            _model: std::marker::PhantomData,
        }
    }
}

impl<M: DsmModel> Default for DpExtender<M> {
    fn default() -> Self {
        Self::new()
    }
}
