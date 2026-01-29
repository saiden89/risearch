use crate::dsm::DsmModel;
use crate::seq::Sequence;
use crate::types::Base;
use log::trace;
use smallvec::SmallVec;

mod core;
mod init;
mod traceback;

use core::dp_main_loop_generic;
use init::{add_e, init_limited_cols, init_limited_rows, update_best_with_term};
use traceback::traceback;

/// Stack-allocated trace buffer. 64 ops covers most extensions without heap allocation.
/// DpOp is 1 byte, so 64 * 1 = 64 bytes on stack.
pub type TracebackPath = SmallVec<[DpOp; 64]>;

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
    _model: std::marker::PhantomData<M>,
}

/// Gap base index constant
const GAP: usize = Base::Gap as usize;
const DSM_DIM: usize = 6;

#[inline]
fn max_dsm_stack<M: DsmModel>() -> i32 {
    let mut max = i32::MIN;
    for q1 in 0..DSM_DIM {
        for q2 in 0..DSM_DIM {
            for t1 in 0..DSM_DIM {
                for t2 in 0..DSM_DIM {
                    let val = M::stack_idx(q1, q2, t1, t2);
                    if val > max {
                        max = val;
                    }
                }
            }
        }
    }
    max
}

#[inline]
fn max_dsm_terminal<M: DsmModel>() -> i32 {
    let mut max = i32::MIN;
    for q in 0..DSM_DIM {
        for t in 0..DSM_DIM {
            let left = M::stack_idx(GAP, q, GAP, t);
            let right = M::stack_idx(q, GAP, t, GAP);
            let val = if left > right { left } else { right };
            if val > max {
                max = val;
            }
        }
    }
    max
}

impl<'a, M: DsmModel> DpView<'a, M> {
    #[inline(always)]
    fn base_or_gap(seq: &Sequence, pos: usize) -> Base {
        if pos >= seq.len() {
            Base::Gap
        } else {
            seq[pos]
        }
    }

    #[inline(always)]
    fn left_base(seq: &Sequence, anchor: usize, offset: usize) -> Base {
        if offset > anchor {
            Base::Gap
        } else {
            seq[anchor - offset]
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
    ) -> Self {
        Self {
            query,
            target,
            q_anchor: q_start,
            t_anchor: t_start,
            dir: ExtendDir::Left,
            q_len: (q_start + 1).min(max_ext),
            t_len: (target.len() - t_start).min(max_ext),
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
    ) -> Self {
        Self {
            query,
            target,
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
            ExtendDir::Left => M::stack_idx(q2, q1, t2, t1),
            // Right: DSM[prev, curr, prev, curr] - stacking toward 3'
            ExtendDir::Right => M::stack_idx(q1, q2, t1, t2),
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
            ExtendDir::Left => M::stack_idx(GAP, self.q(i), GAP, self.t(j)),
            ExtendDir::Right => M::stack_idx(self.q(i), GAP, self.t(j), GAP),
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

/// Result of DP extension
#[derive(Debug, Clone)]
pub struct DpExtension {
    pub score: i32,
    pub q_len: usize,
    pub t_len: usize,
    pub trace: TracebackPath,
}

/// Alignment operation for traceback
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DpOp {
    #[default]
    Stop, // Traceback termination
    Match, // Diagonal move - paired bases
    GapQ,  // Gap in query, target base unpaired
    GapT,  // Gap in target, query base unpaired
}

/// Sentinel value for invalid/uninitialized score (~ -1 billion).
/// Safe for arithmetic: MIN_SCORE + energy (±2000) will not overflow/underflow i32.
const MIN_SCORE: i32 = -1_000_000_000;

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
pub struct DpMatrices {
    pub m: ScoreGrid,  // Match/mismatch state
    pub bq: ScoreGrid, // Query bulge (gap in target)
    pub bt: ScoreGrid, // Target bulge (gap in query)
}

impl DpMatrices {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            m: ScoreGrid::new(width, height),
            bq: ScoreGrid::new(width, height),
            bt: ScoreGrid::new(width, height),
        }
    }

    pub fn resize(&mut self, width: usize, height: usize) {
        self.m.resize(width, height);
        self.bq.resize(width, height);
        self.bt.resize(width, height);
    }
}

/// Stateful DP extender with reusable matrices
pub struct DpExtender<M: DsmModel> {
    matrices: DpMatrices,
    /// Reusable traceback buffer - cleared and reused on each extend() call.
    trace_buf: TracebackPath,
    max_stack: i32,
    max_terminal: i32,
    _model: std::marker::PhantomData<M>,
}

impl<M: DsmModel> DpExtender<M> {
    pub fn new() -> Self {
        let max_stack = max_dsm_stack::<M>();
        let max_terminal = max_dsm_terminal::<M>();
        debug_assert!(max_stack >= MIN_SCORE, "invalid DSM max stack");
        debug_assert!(max_terminal >= MIN_SCORE, "invalid DSM max terminal");

        Self {
            matrices: DpMatrices::new(200, 200),
            // SmallVec doesn't need with_capacity for inline storage
            trace_buf: TracebackPath::new(),
            max_stack,
            max_terminal,
            _model: std::marker::PhantomData,
        }
    }
    /// Extend to the left (query 5', target 3')
    pub fn extend_left(
        &mut self,
        query: &Sequence,
        target: &Sequence,
        q_start: usize,
        t_start: usize,
        max_ext: usize,
    ) -> DpExtension {
        let view = DpView::<M>::left(query, target, q_start, t_start, max_ext);
        self.extend(&view)
    }

    /// Extend to the right (query 3', target 5')
    pub fn extend_right(
        &mut self,
        query: &Sequence,
        target: &Sequence,
        q_end: usize,
        t_end: usize,
        max_ext: usize,
    ) -> DpExtension {
        let view = DpView::<M>::right(query, target, q_end, t_end, max_ext);
        self.extend(&view)
    }

    // =========================================================================
    // UNIFIED EXTEND - Direction-agnostic DP using DpView abstraction
    // =========================================================================

    #[cfg_attr(feature = "prof", inline(never))]
    pub fn extend(&mut self, view: &DpView<'_, M>) -> DpExtension {
        let (q_len, t_len) = (view.q_len, view.t_len);
        let max_stack = self.max_stack;
        let max_terminal = self.max_terminal;

        trace!("{} q_len={} t_len={}", view.dir, q_len, t_len);

        // Initial score: terminal penalty for seed boundary
        let mut best_e = view.terminal(0, 0);
        let mut best_i = 0usize;
        let mut best_j = 0usize;

        // Early return
        if q_len <= 1 || t_len <= 1 {
            return DpExtension {
                score: best_e,
                q_len: 0,
                t_len: 0,
                trace: TracebackPath::new(),
            };
        }

        // Resize matrices
        self.matrices.resize(t_len + 1, q_len + 1);
        let DpMatrices { m, bq, bt } = &mut self.matrices;

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

        // Avoid per-iteration branching in view.q/view.t by specializing on direction.
        if view.dir == ExtendDir::Left {
            for i in 0..q_len.min(MAX_EXT) {
                let qi = DpView::<M>::left_base(view.query, view.q_anchor, i).idx();
                unsafe { *q_ptr.add(i) = qi };
            }
            for j in 0..t_len.min(MAX_EXT) {
                let tj = DpView::<M>::right_base(view.target, view.t_anchor, j).idx();
                unsafe { *t_ptr.add(j) = tj };
            }
        } else {
            for i in 0..q_len.min(MAX_EXT) {
                let qi = DpView::<M>::right_base(view.query, view.q_anchor, i).idx();
                unsafe { *q_ptr.add(i) = qi };
            }
            for j in 0..t_len.min(MAX_EXT) {
                let tj = DpView::<M>::left_base(view.target, view.t_anchor, j).idx();
                unsafe { *t_ptr.add(j) = tj };
            }
        }

        // =====================================================================
        // INITIALIZATION - Unconditional writes to avoid stale data
        // =====================================================================
        // Unlike C which uses calloc (zeroed memory), we reuse buffers.
        // We must explicitly set ALL cells that might be read, including NA cells.
        let width = m.width();
        let m_ptr = m.ptr();
        let bq_ptr = bq.ptr();
        let bt_ptr = bt.ptr();

        // SAFETY invariants:
        // - matrices are sized to at least (q_len+1) x (t_len+1)
        // - indices used below are within those bounds
        debug_assert!(width > t_len, "matrix width too small for t_len");
        unsafe {
            let idx = |i: usize, j: usize| -> usize { i * width + j };

            // Corner cells: explicit NA values (like C code)
            *m_ptr.add(idx(0, 0)) = 0;
            *bq_ptr.add(idx(0, 0)) = MIN_SCORE;
            *bt_ptr.add(idx(0, 0)) = MIN_SCORE;
            *m_ptr.add(idx(0, 1)) = MIN_SCORE;
            *bq_ptr.add(idx(0, 1)) = MIN_SCORE;
            *m_ptr.add(idx(1, 0)) = MIN_SCORE;
            *bt_ptr.add(idx(1, 0)) = MIN_SCORE;
            *bq_ptr.add(idx(1, 1)) = MIN_SCORE;
            *bt_ptr.add(idx(1, 1)) = MIN_SCORE;

            // Valid corner values
            *bt_ptr.add(idx(0, 1)) = view.bt_open(0, 1);
            *bq_ptr.add(idx(1, 0)) = view.bq_open(1, 0);
            let m11 = view.match_e(1, 1);
            *m_ptr.add(idx(1, 1)) = m11;
            update_best_with_term(
                &mut best_e,
                &mut best_i,
                &mut best_j,
                m11,
                view.terminal(1, 1),
                1,
                1,
            );

            // Row 0 (Bt only) and Row 1 (M) - unconditional writes
            for k in 2..t_len {
                let prev = *bt_ptr.add(idx(0, k - 1));
                // Always write - use add_e to propagate MIN_SCORE
                let bt_val = add_e(prev, view.bt_ext(k));
                let m_val = add_e(prev, view.m_from_bt(1, k));
                *bt_ptr.add(idx(0, k)) = bt_val;
                *m_ptr.add(idx(0, k)) = MIN_SCORE; // M[0,k] is NA
                *bq_ptr.add(idx(0, k)) = MIN_SCORE; // Bq[0,k] is NA
                *bq_ptr.add(idx(1, k)) = MIN_SCORE; // Bq[1,k] is NA
                *m_ptr.add(idx(1, k)) = m_val;
                update_best_with_term(
                    &mut best_e,
                    &mut best_i,
                    &mut best_j,
                    m_val,
                    view.terminal(1, k),
                    1,
                    k,
                );
            }

            // Col 0 (Bq only) and Col 1 (M) - unconditional writes
            for k in 2..q_len {
                let prev = *bq_ptr.add(idx(k - 1, 0));
                let bq_val = add_e(prev, view.bq_ext(k));
                let m_val = add_e(prev, view.m_from_bq(k, 1));
                *bq_ptr.add(idx(k, 0)) = bq_val;
                *m_ptr.add(idx(k, 0)) = MIN_SCORE; // M[k,0] is NA
                *bt_ptr.add(idx(k, 0)) = MIN_SCORE; // Bt[k,0] is NA
                *bt_ptr.add(idx(k, 1)) = MIN_SCORE; // Bt[k,1] is NA
                *m_ptr.add(idx(k, 1)) = m_val;
                update_best_with_term(
                    &mut best_e,
                    &mut best_i,
                    &mut best_j,
                    m_val,
                    view.terminal(k, 1),
                    k,
                    1,
                );
            }
        }

        // Early return when either axis is too small for row/col 2 cells
        if q_len <= 2 || t_len <= 2 {
            return DpExtension {
                score: best_e,
                q_len: best_i,
                t_len: best_j,
                trace: TracebackPath::new(),
            };
        }

        // Cell (2,2) init: bridge corner to limited rows/cols
        let m11_val = m.get(1, 1);
        let bt12 = add_e(m11_val, view.bt_open(1, 2));
        let bq21 = add_e(m11_val, view.bq_open(2, 1));
        let m22 = add_e(m11_val, view.match_e(2, 2));
        bt.set(1, 2, bt12);
        bq.set(2, 1, bq21);
        m.set(2, 2, m22);
        update_best_with_term(
            &mut best_e,
            &mut best_i,
            &mut best_j,
            m22,
            view.terminal(2, 2),
            2,
            2,
        );

        let m12 = m.get(1, 2);
        let m21 = m.get(2, 1);
        bq.set(2, 2, add_e(m12, view.bq_open(2, 2)));
        bt.set(2, 2, add_e(m21, view.bt_open(2, 2)));

        // =======================================================================
        // LIMITED ROWS/COLUMNS - unified via macro (score-only)
        // =======================================================================

        init_limited_rows::<M>(
            view,
            q_ptr,
            t_ptr,
            m_ptr,
            bt_ptr,
            bq_ptr,
            width,
            q_len,
            t_len,
            &mut best_e,
            &mut best_i,
            &mut best_j,
        );

        init_limited_cols::<M>(
            view,
            q_ptr,
            t_ptr,
            m_ptr,
            bq_ptr,
            bt_ptr,
            width,
            q_len,
            &mut best_e,
            &mut best_i,
            &mut best_j,
        );

        // =======================================================================
        // MAIN DP LOOP (i >= 3, j >= 3) - OPTIMIZED
        // =======================================================================
        //
        // Optimizations applied:
        // 1. Precompute Q/T base indices into stack arrays (eliminates branch per access)
        // 2. Cache row offsets (eliminates multiplication per cell)
        // 3. Direct slice access (eliminates method call overhead)
        // 4. Raw DSM lookup (bypasses abstraction layers)

        if q_len >= 3 && t_len >= 3 {
            if view.dir == ExtendDir::Left {
                dp_main_loop_generic::<true, M>(
                    q_ptr,
                    t_ptr,
                    m,
                    bq,
                    bt,
                    width,
                    q_len,
                    t_len,
                    max_stack,
                    max_terminal,
                    &mut best_e,
                    &mut best_i,
                    &mut best_j,
                );
            } else {
                dp_main_loop_generic::<false, M>(
                    q_ptr,
                    t_ptr,
                    m,
                    bq,
                    bt,
                    width,
                    q_len,
                    t_len,
                    max_stack,
                    max_terminal,
                    &mut best_e,
                    &mut best_i,
                    &mut best_j,
                );
            }
        }

        // =======================================================================
        // TRACEBACK - Reconstructed from scores (like C implementation)
        // =======================================================================
        // C doesn't store traceback during DP. It reconstructs the path by
        // comparing scores to determine which transition was taken.

        trace!(
            "{} TB start: best=({},{}) score={}",
            view.dir, best_i, best_j, best_e
        );

        // Reuse traceback buffer (cleared each call, capacity preserved)
        self.trace_buf.clear();
        traceback::<M>(view, m, bq, bt, best_i, best_j, &mut self.trace_buf);

        trace!(
            "{} result: score={} q_len={} t_len={} trace={:?}",
            view.dir, best_e, best_i, best_j, self.trace_buf
        );
        // Move trace out using mem::take (zero-copy) - SmallVec::default() is inline-empty
        DpExtension {
            score: best_e,
            q_len: best_i,
            t_len: best_j,
            trace: std::mem::take(&mut self.trace_buf),
        }
    }
}

impl<M: DsmModel> Default for DpExtender<M> {
    fn default() -> Self {
        Self::new()
    }
}
