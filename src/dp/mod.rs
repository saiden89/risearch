use crate::alignment::Pairing;
use crate::dsm::{stack_with_penalty, DsmModel};
use crate::seq::Sequence;
use crate::types::Base;
use log::trace;
use smallvec::SmallVec;

mod core;
mod init;
mod traceback;

use core::dp_main_loop_generic;
use init::init_frontier;
use traceback::traceback;

/// Maximum extension length for precomputed index arrays.
/// Matches the typical max_ext parameter (100-200 bases).
const MAX_EXT: usize = 256;

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
    query: &'a Sequence,
    target: &'a Sequence,
    q_anchor: usize,
    t_anchor: usize,
    pub dir: ExtendDir,
    pub q_len: usize,
    pub t_len: usize,
    penalty: i32,
    _model: std::marker::PhantomData<M>,
}

/// Gap base index constant
const GAP: usize = Base::Gap as usize;

impl<'a, M: DsmModel> DpView<'a, M> {
    #[inline(always)]
    fn base_or_gap(seq: &Sequence, pos: usize) -> Base {
        if pos >= seq.len() {
            Base::Gap
        } else {
            // SAFETY: bounds checked above
            unsafe { seq.get_unchecked(pos) }
        }
    }

    #[inline(always)]
    fn left_base(seq: &Sequence, anchor: usize, offset: usize) -> Base {
        if offset > anchor {
            Base::Gap
        } else {
            // SAFETY: offset <= anchor, and anchor < seq.len() by construction
            unsafe { seq.get_unchecked(anchor - offset) }
        }
    }

    #[inline(always)]
    fn right_base(seq: &Sequence, anchor: usize, offset: usize) -> Base {
        Self::base_or_gap(seq, anchor + offset)
    }

    /// Create a left extension view (query toward 5', target toward 3')
    pub fn left(
        query: &'a Sequence,
        target: &'a Sequence,
        q_start: usize,
        t_start: usize,
        max_ext: usize,
        penalty: i32,
    ) -> Self {
        Self {
            query,
            target,
            q_anchor: q_start,
            t_anchor: t_start,
            dir: ExtendDir::Left,
            q_len: (q_start + 1).min(max_ext),
            t_len: (target.len() - t_start).min(max_ext),
            penalty,
            _model: std::marker::PhantomData,
        }
    }

    /// Create a right extension view (query toward 3', target toward 5')
    pub fn right(
        query: &'a Sequence,
        target: &'a Sequence,
        q_end: usize,
        t_end: usize,
        max_ext: usize,
        penalty: i32,
    ) -> Self {
        Self {
            query,
            target,
            q_anchor: q_end,
            t_anchor: t_end,
            dir: ExtendDir::Right,
            q_len: (query.len() - q_end).min(max_ext),
            t_len: (t_end + 1).min(max_ext),
            penalty,
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

    /// Get target base index at DP position j (0 = anchor)
    #[inline(always)]
    pub fn t(&self, j: usize) -> usize {
        match self.dir {
            ExtendDir::Left => Self::right_base(self.target, self.t_anchor, j).idx(),
            ExtendDir::Right => Self::left_base(self.target, self.t_anchor, j).idx(),
        }
    }

    /// Stacking energy with direction-aware ordering.
    ///
    /// Always called with arguments in (prev, curr, prev, curr) order.
    /// Left extension internally swaps to (curr, prev, curr, prev).
    #[inline(always)]
    pub fn e(&self, q1: usize, q2: usize, t1: usize, t2: usize) -> i32 {
        match self.dir {
            // Left: DSM[curr, prev, curr, prev] - stacking toward 5'
            ExtendDir::Left => stack_with_penalty::<M>(q2, q1, t2, t1, self.penalty),
            // Right: DSM[prev, curr, prev, curr] - stacking toward 3'
            ExtendDir::Right => stack_with_penalty::<M>(q1, q2, t1, t2, self.penalty),
        }
    }

    /// Terminal penalty at position (i, j) - stacking with gap boundary.
    ///
    /// Terminal energy is direction-specific (cannot use view.e() swap):
    /// - LEFT: GAP at 5' side, base at 3' → DSM[GAP, q, GAP, t]
    /// - RIGHT: base at 5' side, GAP at 3' → DSM[q, GAP, t, GAP]
    #[inline(always)]
    pub fn terminal(&self, i: usize, j: usize) -> i32 {
        match self.dir {
            ExtendDir::Left => {
                stack_with_penalty::<M>(GAP, self.q(i), GAP, self.t(j), self.penalty)
            }
            ExtendDir::Right => {
                stack_with_penalty::<M>(self.q(i), GAP, self.t(j), GAP, self.penalty)
            }
        }
    }

    /// Match/mismatch energy at (i, j) from diagonal (i-1, j-1)
    #[inline(always)]
    pub fn match_e(&self, i: usize, j: usize) -> i32 {
        self.e(self.q(i - 1), self.q(i), self.t(j - 1), self.t(j))
    }

    // =========================================================================
    // DP TRANSITION HELPERS - exactly match the DP recurrence semantics
    // =========================================================================
    // All methods take arguments in (prev, curr) order; view.e() swaps for LEFT.
    // These match the actual energy lookups in the DP main loop.

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
// SCORE-ONLY GRID - Optimized for DP hot loop (matches C implementation)
// =============================================================================

/// Score-only grid for DP. No traceback storage - reconstructed post-hoc.
/// This matches the C implementation which stores only M, Bq, Bt scores.
pub struct ScoreGrid {
    data: Vec<i32>,
    width: usize,
}

impl ScoreGrid {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            data: vec![MIN_SCORE; width * height],
            width,
        }
    }

    /// Resize without clearing. Caller must ensure all accessed cells are initialized.
    /// This matches C behavior where matrices are allocated once and reused.
    #[inline]
    pub fn resize(&mut self, width: usize, height: usize) {
        let new_len = width * height;
        if self.data.len() < new_len {
            self.data.resize(new_len, MIN_SCORE);
        }
        self.width = width;
        // NOTE: No fill() - cells are initialized on demand during DP.
        // This is critical for performance (avoids O(n²) overhead per extension).
    }

    #[inline(always)]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        i * self.width + j
    }

    #[inline(always)]
    pub fn get(&self, i: usize, j: usize) -> i32 {
        let idx = i * self.width + j;
        debug_assert!(
            idx < self.data.len(),
            "ScoreOnlyGrid::get out of bounds: ({}, {}) idx={} len={}",
            i,
            j,
            idx,
            self.data.len()
        );
        unsafe { *self.data.get_unchecked(idx) }
    }

    #[inline(always)]
    pub fn set(&mut self, i: usize, j: usize, val: i32) {
        let idx = i * self.width + j;
        debug_assert!(
            idx < self.data.len(),
            "ScoreOnlyGrid::set out of bounds: ({}, {}) idx={} len={}",
            i,
            j,
            idx,
            self.data.len()
        );
        unsafe { *self.data.get_unchecked_mut(idx) = val }
    }

    #[inline(always)]
    pub fn ptr(&mut self) -> *mut i32 {
        self.data.as_mut_ptr()
    }

    pub fn width(&self) -> usize {
        self.width
    }
}

// DP EXTENDER - Stateful extension with reusable matrices (score-only)
// =============================================================================

/// Score-only DP matrices for extension (matches C implementation).
/// Traceback is reconstructed post-hoc by comparing scores.
pub(super) struct DpMatrices {
    pub(super) m: ScoreGrid,  // Match/mismatch state
    pub(super) bq: ScoreGrid, // Query bulge (gap in target)
    pub(super) bt: ScoreGrid, // Target bulge (gap in query)
}

impl DpMatrices {
    fn new(width: usize, height: usize) -> Self {
        Self {
            m: ScoreGrid::new(width, height),
            bq: ScoreGrid::new(width, height),
            bt: ScoreGrid::new(width, height),
        }
    }

    pub(super) fn resize(&mut self, width: usize, height: usize) {
        self.m.resize(width, height);
        self.bq.resize(width, height);
        self.bt.resize(width, height);
    }
}

/// Stateful DP extender with reusable matrices
pub struct DpExtender<M: DsmModel> {
    matrices: DpMatrices,
    penalty: i32,
    _model: std::marker::PhantomData<M>,
}

/// Guard holding a reference to DP matrices after forward pass.
/// While this exists, the extender cannot be used for another extension
/// (borrow checker enforces this via the lifetime on `matrices`).
pub struct ExtendResult<'a, M: DsmModel> {
    matrices: &'a DpMatrices,
    pub score: i32,
    pub q_len: usize,
    pub t_len: usize,
    _model: std::marker::PhantomData<M>,
}

impl<M: DsmModel> ExtendResult<'_, M> {
    /// Run traceback to reconstruct alignment as Pairings.
    /// Only call when alignment output is needed (skip for Minimal format).
    pub fn traceback(&self, view: &DpView<'_, M>) -> SmallVec<[Pairing; 64]> {
        let mut out = SmallVec::new();
        traceback::<M>(
            view,
            &self.matrices.m,
            &self.matrices.bq,
            &self.matrices.bt,
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

    pub fn with_penalty(penalty: i32) -> Self {
        Self {
            matrices: DpMatrices::new(200, 200),
            penalty,
            _model: std::marker::PhantomData,
        }
    }

    // =========================================================================
    // UNIFIED EXTEND - Direction-agnostic DP using DpView abstraction
    // =========================================================================

    #[cfg_attr(feature = "prof", inline(never))]
    /// Run DP forward pass and return an ExtendResult guard.
    /// Call `.traceback()` on the result if alignment is needed.
    pub fn extend(&mut self, view: &DpView<'_, M>) -> ExtendResult<'_, M> {
        let (q_len, t_len) = (view.q_len, view.t_len);

        trace!("{} q_len={} t_len={}", view.dir, q_len, t_len);

        // Initial score: terminal penalty for seed boundary
        let mut best = BestScore::new(view.terminal(0, 0));

        // Early return
        if q_len <= 1 || t_len <= 1 {
            return ExtendResult {
                matrices: &self.matrices,
                score: best.score,
                q_len: 0,
                t_len: 0,
                _model: std::marker::PhantomData,
            };
        }

        // Resize matrices
        self.matrices.resize(t_len + 1, q_len + 1);

        // =====================================================================
        // PRECOMPUTE Q/T BASE INDICES (used by init AND main loop)
        // =====================================================================
        // This eliminates repeated `view.q(i)` and `view.t(j)` calls
        // which have a match on direction each time.

        // Avoid per-call zeroing of these stack arrays (keeps memset out of hot path).
        let mut q_idx = std::mem::MaybeUninit::<[usize; MAX_EXT]>::uninit();
        let mut t_idx = std::mem::MaybeUninit::<[usize; MAX_EXT]>::uninit();

        // SAFETY invariants:
        // - We only read indices we explicitly write below.
        // - q_len and t_len are <= MAX_EXT.
        let q_ptr = q_idx.as_mut_ptr() as *mut usize;
        let t_ptr = t_idx.as_mut_ptr() as *mut usize;

        debug_assert!(q_len <= MAX_EXT, "q_len exceeds precomputed index capacity");
        debug_assert!(t_len <= MAX_EXT, "t_len exceeds precomputed index capacity");

        for i in 0..q_len.min(MAX_EXT) {
            unsafe { *q_ptr.add(i) = view.q(i) };
        }
        for j in 0..t_len.min(MAX_EXT) {
            unsafe { *t_ptr.add(j) = view.t(j) };
        }

        // =====================================================================
        // INITIALIZATION - Unconditional writes to avoid stale data
        // =====================================================================
        // Unlike C which uses calloc (zeroed memory), we reuse buffers.
        // We must explicitly set ALL cells that might be read, including NA cells.

        let has_main_region = init_frontier::<M>(
            view,
            q_ptr,
            t_ptr,
            &mut self.matrices,
            q_len,
            t_len,
            &mut best,
        );
        if !has_main_region {
            return ExtendResult {
                matrices: &self.matrices,
                score: best.score,
                q_len: best.i,
                t_len: best.j,
                _model: std::marker::PhantomData,
            };
        }

        // =======================================================================
        // MAIN DP LOOP (i >= 3, j >= 3) - OPTIMIZED
        // =======================================================================
        //
        // Optimizations applied:
        // 1. Precompute Q/T base indices into stack arrays (eliminates branch per access)
        // 2. Cache row offsets (eliminates multiplication per cell)
        // 3. Direct slice access (eliminates method call overhead)
        // 4. Raw DSM lookup (bypasses abstraction layers)

        if view.dir == ExtendDir::Left {
            dp_main_loop_generic::<true, M>(
                q_ptr,
                t_ptr,
                &mut self.matrices,
                q_len,
                t_len,
                self.penalty,
                &mut best,
            );
        } else {
            dp_main_loop_generic::<false, M>(
                q_ptr,
                t_ptr,
                &mut self.matrices,
                q_len,
                t_len,
                self.penalty,
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
            matrices: &self.matrices,
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
