use crate::dsm::{EnergyModel, dsm_lookup_raw};
use crate::seq::Seq;
use crate::types::Base;
use log::trace;
use smallvec::SmallVec;

/// Stack-allocated trace buffer. 64 ops covers most extensions without heap allocation.
/// DpOp is 1 byte, so 64 * 1 = 64 bytes on stack.
pub type TraceVec = SmallVec<[DpOp; 64]>;

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
pub struct DpView<'a> {
    query: &'a Seq<'a>,
    target: &'a Seq<'a>,
    q_anchor: usize,
    t_anchor: usize,
    pub dir: ExtendDir,
    pub q_len: usize,
    pub t_len: usize,
}

/// Gap base index constant
const GAP: usize = Base::Gap as usize;

impl<'a> DpView<'a> {
    /// Create a left extension view (query toward 5', target toward 3')
    pub fn left(
        query: &'a Seq<'a>,
        target: &'a Seq<'a>,
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
        }
    }

    /// Create a right extension view (query toward 3', target toward 5')
    pub fn right(
        query: &'a Seq<'a>,
        target: &'a Seq<'a>,
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
        }
    }

    /// Get query base index at DP position i (0 = anchor)
    #[inline(always)]
    pub fn q(&self, i: usize) -> usize {
        match self.dir {
            ExtendDir::Left => self.query.left(self.q_anchor, i).idx(),
            ExtendDir::Right => self.query.right(self.q_anchor, i).idx(),
        }
    }

    /// Get target base index at DP position j (0 = anchor)
    #[inline(always)]
    pub fn t(&self, j: usize) -> usize {
        match self.dir {
            ExtendDir::Left => self.target.right(self.t_anchor, j).idx(),
            ExtendDir::Right => self.target.left(self.t_anchor, j).idx(),
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
            ExtendDir::Left => EnergyModel::T04.stack_idx(q2, q1, t2, t1),
            // Right: DSM[prev, curr, prev, curr] - stacking toward 3'
            ExtendDir::Right => EnergyModel::T04.stack_idx(q1, q2, t1, t2),
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
            ExtendDir::Left => EnergyModel::T04.stack_idx(GAP, self.q(i), GAP, self.t(j)),
            ExtendDir::Right => EnergyModel::T04.stack_idx(self.q(i), GAP, self.t(j), GAP),
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
    pub trace: TraceVec,
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
pub struct ScoreOnlyGrid {
    data: Vec<i32>,
    width: usize,
}

impl ScoreOnlyGrid {
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

    #[inline(always)]
    pub fn width(&self) -> usize {
        self.width
    }
}

// =============================================================================
// DP EXTENDER - Stateful extension with reusable matrices (score-only)
// =============================================================================

/// Score-only DP matrices for extension (matches C implementation).
/// Traceback is reconstructed post-hoc by comparing scores.
pub struct DpMatrices {
    pub m: ScoreOnlyGrid,  // Match/mismatch state
    pub bq: ScoreOnlyGrid, // Query bulge (gap in target)
    pub bt: ScoreOnlyGrid, // Target bulge (gap in query)
}

impl DpMatrices {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            m: ScoreOnlyGrid::new(width, height),
            bq: ScoreOnlyGrid::new(width, height),
            bt: ScoreOnlyGrid::new(width, height),
        }
    }

    pub fn resize(&mut self, width: usize, height: usize) {
        self.m.resize(width, height);
        self.bq.resize(width, height);
        self.bt.resize(width, height);
    }
}

/// Stateful DP extender with reusable matrices
pub struct DpExtender {
    matrices: DpMatrices,
    /// Reusable traceback buffer - cleared and reused on each extend() call.
    trace_buf: TraceVec,
}

// =============================================================================
// INIT HELPERS - Reduce code duplication in DP initialization
// =============================================================================

/// Pick best value from two sources.
/// Uses std::cmp::max which LLVM compiles to branchless CMOV.
#[inline(always)]
fn pick_best(val_a: i32, val_b: i32) -> i32 {
    std::cmp::max(val_a, val_b)
}

#[inline(always)]
fn update_best_with_term(
    best_e: &mut i32,
    best_i: &mut usize,
    best_j: &mut usize,
    val: i32,
    term: i32,
    i: usize,
    j: usize,
) {
    if val > MIN_SCORE {
        let curr = val + term;
        if curr > *best_e {
            *best_e = curr;
            *best_i = i;
            *best_j = j;
        }
    }
}

/// Helper: add energy if base is valid (not MIN_SCORE).
/// Simple branch - LLVM optimizes to CMOV when beneficial.
#[inline(always)]
fn add_e(base: i32, energy: i32) -> i32 {
    if base > MIN_SCORE {
        base + energy
    } else {
        MIN_SCORE
    }
}

/// max of 2 values - branchless via std::cmp::max (compiles to CMOV)
#[inline(always)]
fn max2(a: i32, b: i32) -> i32 {
    std::cmp::max(a, b)
}

/// max of 3 values - branchless
#[inline(always)]
fn max3(a: i32, b: i32, c: i32) -> i32 {
    std::cmp::max(std::cmp::max(a, b), c)
}

/// Initialize limited rows (t_len axis) - score-only version.
/// All writes are unconditional to avoid reading stale data.
fn init_limited_rows(
    view: &DpView<'_>,
    q_ptr: *const usize,
    t_ptr: *const usize,
    m_ptr: *mut i32,
    bt_ptr: *mut i32,
    bq_ptr: *mut i32,
    width: usize,
    q_len: usize,
    t_len: usize,
    best_e: &mut i32,
    best_i: &mut usize,
    best_j: &mut usize,
) {
    // SAFETY invariants:
    // - q_ptr and t_ptr are valid for indices [0, q_len) and [0, t_len)
    // - q_len and t_len are >= 3 (we index 1 and 2)
    // - m_ptr/bt_ptr/bq_ptr point to matrices sized at least (q_len+1) x (t_len+1)
    let left = view.dir == ExtendDir::Left;
    let stack = |q1: usize, q2: usize, t1: usize, t2: usize| -> i32 {
        if left {
            dsm_lookup_raw(q2, q1, t2, t1)
        } else {
            dsm_lookup_raw(q1, q2, t1, t2)
        }
    };

    debug_assert!(q_len > 2 && t_len > 2, "limited rows require q_len/t_len >= 3");
    debug_assert!(q_len <= MAX_EXT, "q_len exceeds precomputed index capacity");
    debug_assert!(t_len <= MAX_EXT, "t_len exceeds precomputed index capacity");
    unsafe {
        let idx = |i: usize, j: usize| -> usize { i * width + j };
        let qi1 = *q_ptr.add(1);
        let qi2 = *q_ptr.add(2);

        // Rolling values for row 1 and 2
        let mut m1_prev = *m_ptr.add(idx(1, 2));  // m(1, k-1)
        let mut bt1_prev = *bt_ptr.add(idx(1, 2)); // bt(1, k-1)
        let mut m2_prev = *m_ptr.add(idx(2, 2));  // m(2, k-1)
        let mut bt2_prev = *bt_ptr.add(idx(2, 2)); // bt(2, k-1)

        for k in 3..t_len {
            let tj = *t_ptr.add(k);
            let tj_prev = *t_ptr.add(k - 1);

            // Primary[1,k] = Bt
            let from_m = add_e(m1_prev, stack(qi1, GAP, tj_prev, tj));
            let from_b = add_e(bt1_prev, stack(GAP, GAP, tj_prev, tj));
            let bt1 = pick_best(from_m, from_b);
            *bt_ptr.add(idx(1, k)) = bt1;

            // M[2,k]
            let from_m = add_e(m1_prev, stack(qi1, qi2, tj_prev, tj));
            let from_b = add_e(bt1_prev, stack(GAP, qi2, tj_prev, tj));
            let m2 = pick_best(from_m, from_b);
            *m_ptr.add(idx(2, k)) = m2;
            let term = if left {
                dsm_lookup_raw(GAP, qi2, GAP, tj)
            } else {
                dsm_lookup_raw(qi2, GAP, tj, GAP)
            };
            update_best_with_term(best_e, best_i, best_j, m2, term, 2, k);

            // Secondary[2,k] = Bq
            let m1k = *m_ptr.add(idx(1, k));
            *bq_ptr.add(idx(2, k)) = add_e(m1k, stack(qi1, qi2, tj, GAP));

            // Primary[2,k] = Bt
            let from_m = add_e(m2_prev, stack(qi2, GAP, tj_prev, tj));
            let from_b = add_e(bt2_prev, stack(GAP, GAP, tj_prev, tj));
            let bt2 = pick_best(from_m, from_b);
            *bt_ptr.add(idx(2, k)) = bt2;

            m1_prev = *m_ptr.add(idx(1, k));
            bt1_prev = bt1;
            m2_prev = m2;
            bt2_prev = bt2;
        }
    }
}

/// Initialize limited columns (q_len axis) - score-only version.
/// All writes are unconditional to avoid reading stale data.
fn init_limited_cols(
    view: &DpView<'_>,
    q_ptr: *const usize,
    t_ptr: *const usize,
    m_ptr: *mut i32,
    bq_ptr: *mut i32,
    bt_ptr: *mut i32,
    width: usize,
    q_len: usize,
    best_e: &mut i32,
    best_i: &mut usize,
    best_j: &mut usize,
) {
    // SAFETY invariants:
    // - q_ptr and t_ptr are valid for indices [0, q_len) and [0, t_len)
    // - q_len and t_len are >= 3 (we index 1 and 2)
    // - m_ptr/bq_ptr/bt_ptr point to matrices sized at least (q_len+1) x (t_len+1)
    let left = view.dir == ExtendDir::Left;
    let stack = |q1: usize, q2: usize, t1: usize, t2: usize| -> i32 {
        if left {
            dsm_lookup_raw(q2, q1, t2, t1)
        } else {
            dsm_lookup_raw(q1, q2, t1, t2)
        }
    };

    debug_assert!(q_len > 2, "limited cols require q_len >= 3");
    debug_assert!(q_len <= MAX_EXT, "q_len exceeds precomputed index capacity");
    unsafe {
        let idx = |i: usize, j: usize| -> usize { i * width + j };
        let tj1 = *t_ptr.add(1);
        let tj2 = *t_ptr.add(2);

        for k in 3..q_len {
            let qi = *q_ptr.add(k);
            let qi_prev = *q_ptr.add(k - 1);

            // Primary[ k,1 ] = Bq
            let from_m = add_e(*m_ptr.add(idx(k - 1, 1)), stack(qi_prev, qi, tj1, GAP));
            let from_b = add_e(*bq_ptr.add(idx(k - 1, 1)), stack(qi_prev, qi, GAP, GAP));
            *bq_ptr.add(idx(k, 1)) = pick_best(from_m, from_b);

            // M[ k,2 ]
            let from_m = add_e(*m_ptr.add(idx(k - 1, 1)), stack(qi_prev, qi, tj1, tj2));
            let from_b = add_e(*bq_ptr.add(idx(k - 1, 1)), stack(qi_prev, qi, GAP, tj2));
            let val = pick_best(from_m, from_b);
            *m_ptr.add(idx(k, 2)) = val;
            let term = if left {
                dsm_lookup_raw(GAP, qi, GAP, tj2)
            } else {
                dsm_lookup_raw(qi, GAP, tj2, GAP)
            };
            update_best_with_term(best_e, best_i, best_j, val, term, k, 2);

            // Secondary[ k,2 ] = Bt
            let m_k1 = *m_ptr.add(idx(k, 1));
            *bt_ptr.add(idx(k, 2)) = add_e(m_k1, stack(qi, GAP, tj1, tj2));

            // Primary[ k,2 ] = Bq
            let from_m = add_e(*m_ptr.add(idx(k - 1, 2)), stack(qi_prev, qi, tj2, GAP));
            let from_b = add_e(*bq_ptr.add(idx(k - 1, 2)), stack(qi_prev, qi, GAP, GAP));
            *bq_ptr.add(idx(k, 2)) = pick_best(from_m, from_b);
        }
    }
}

impl DpExtender {
    pub fn new() -> Self {
        Self {
            matrices: DpMatrices::new(200, 200),
            // SmallVec doesn't need with_capacity for inline storage
            trace_buf: TraceVec::new(),
        }
    }
    /// Extend to the left (query 5', target 3')
    pub fn extend_left(
        &mut self,
        query: &Seq,
        target: &Seq,
        q_start: usize,
        t_start: usize,
        max_ext: usize,
    ) -> DpExtension {
        let view = DpView::left(query, target, q_start, t_start, max_ext);
        self.extend(&view)
    }

    /// Extend to the right (query 3', target 5')
    pub fn extend_right(
        &mut self,
        query: &Seq,
        target: &Seq,
        q_end: usize,
        t_end: usize,
        max_ext: usize,
    ) -> DpExtension {
        let view = DpView::right(query, target, q_end, t_end, max_ext);
        self.extend(&view)
    }

    // =========================================================================
    // UNIFIED EXTEND - Direction-agnostic DP using DpView abstraction
    // =========================================================================

    #[cfg_attr(feature = "prof", inline(never))]
    pub fn extend(&mut self, view: &DpView) -> DpExtension {
        let (q_len, t_len) = (view.q_len, view.t_len);

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
                trace: TraceVec::new(),
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
                let qi = view.query.left(view.q_anchor, i).idx();
                unsafe { *q_ptr.add(i) = qi };
            }
            for j in 0..t_len.min(MAX_EXT) {
                let tj = view.target.right(view.t_anchor, j).idx();
                unsafe { *t_ptr.add(j) = tj };
            }
        } else {
            for i in 0..q_len.min(MAX_EXT) {
                let qi = view.query.right(view.q_anchor, i).idx();
                unsafe { *q_ptr.add(i) = qi };
            }
            for j in 0..t_len.min(MAX_EXT) {
                let tj = view.target.left(view.t_anchor, j).idx();
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
        debug_assert!(width >= t_len + 1, "matrix width too small for t_len");
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
                trace: TraceVec::new(),
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

        init_limited_rows(
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

        init_limited_cols(
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

        // Skip if too short for main loop
        if q_len >= 3 && t_len >= 3 {
            // Get direct pointer access to score data - eliminates struct indirection

            // =================================================================
            // DIRECTION-SPECIFIC MAIN LOOP (SCORES ONLY - like C)
            // =================================================================
            // C stores only scores, reconstructs traceback from scores post-hoc.
            // This reduces memory bandwidth: 3 i32 per cell instead of 6.

            macro_rules! dp_main_loop {
                ($match_e:expr, $m_from_bq:expr, $m_from_bt:expr, $term:expr,
                 $bq_open:expr, $bq_ext:expr, $bt_open:expr, $bt_ext:expr) => {
                    // SAFETY: All indices are bounded by grid dimensions.
                    unsafe {
                        let m_ptr = m.ptr();
                        let bq_ptr = bq.ptr();
                        let bt_ptr = bt.ptr();
                        for i in 3..q_len {
                            let row_i = i * width;
                            let row_prev = (i - 1) * width;
                            let qi = *q_ptr.add(i);
                            let qi_prev = *q_ptr.add(i - 1);

                            for j in 3..t_len {
                                let diag_idx = row_prev + j - 1;
                                let up_idx = row_prev + j;
                                let left_idx = row_i + j - 1;
                                let curr_idx = row_i + j;
                                let tj = *t_ptr.add(j);
                                let tj_prev = *t_ptr.add(j - 1);

                                // M[i,j] = max3(M[i-1,j-1]+match, Bq[i-1,j-1]+m_from_bq, Bt[i-1,j-1]+m_from_bt)
                                let m_diag = *m_ptr.add(diag_idx);
                                let bq_diag = *bq_ptr.add(diag_idx);
                                let bt_diag = *bt_ptr.add(diag_idx);

                                let s_mm = add_e(m_diag, $match_e(qi, qi_prev, tj, tj_prev));
                                let s_mq = add_e(bq_diag, $m_from_bq(qi, qi_prev, tj));
                                let s_mt = add_e(bt_diag, $m_from_bt(qi, tj, tj_prev));

                                // max3 like C - no traceback storage
                                let val_m = max3(s_mm, s_mq, s_mt);

                                // Update best (with terminal penalty)
                                update_best_with_term(
                                    &mut best_e,
                                    &mut best_i,
                                    &mut best_j,
                                    val_m,
                                    $term(qi, tj),
                                    i,
                                    j,
                                );

                                *m_ptr.add(curr_idx) = val_m;

                                // Bq[i,j] = max(M[i-1,j]+bq_open, Bq[i-1,j]+bq_ext)
                                let m_up = *m_ptr.add(up_idx);
                                let bq_up = *bq_ptr.add(up_idx);
                                let s_qm = add_e(m_up, $bq_open(qi, qi_prev, tj));
                                let s_qq = add_e(bq_up, $bq_ext(qi, qi_prev));
                                *bq_ptr.add(curr_idx) = max2(s_qm, s_qq);

                                // Bt[i,j] = max(M[i,j-1]+bt_open, Bt[i,j-1]+bt_ext)
                                let m_left = *m_ptr.add(left_idx);
                                let bt_left = *bt_ptr.add(left_idx);
                                let s_tm = add_e(m_left, $bt_open(qi, tj, tj_prev));
                                let s_tt = add_e(bt_left, $bt_ext(tj, tj_prev));
                                *bt_ptr.add(curr_idx) = max2(s_tm, s_tt);
                            }
                        }
                    }
                };
            }

            if view.dir == ExtendDir::Left {
                // LEFT: DSM[curr, prev, curr, prev]
                dp_main_loop!(
                    |qi, qi_prev, tj, tj_prev| dsm_lookup_raw(qi, qi_prev, tj, tj_prev),
                    |qi, qi_prev, tj| dsm_lookup_raw(qi, qi_prev, tj, GAP),
                    |qi, tj, tj_prev| dsm_lookup_raw(qi, GAP, tj, tj_prev),
                    |qi, tj| dsm_lookup_raw(GAP, qi, GAP, tj),
                    |qi, qi_prev, tj| dsm_lookup_raw(qi, qi_prev, GAP, tj),
                    |qi, qi_prev| dsm_lookup_raw(qi, qi_prev, GAP, GAP),
                    |qi, tj, tj_prev| dsm_lookup_raw(GAP, qi, tj, tj_prev),
                    |tj, tj_prev| dsm_lookup_raw(GAP, GAP, tj, tj_prev)
                );
            } else {
                // RIGHT: DSM[prev, curr, prev, curr]
                dp_main_loop!(
                    |qi, qi_prev, tj, tj_prev| dsm_lookup_raw(qi_prev, qi, tj_prev, tj),
                    |qi, qi_prev, tj| dsm_lookup_raw(qi_prev, qi, GAP, tj),
                    |qi, tj, tj_prev| dsm_lookup_raw(GAP, qi, tj_prev, tj),
                    |qi, tj| dsm_lookup_raw(qi, GAP, tj, GAP),
                    |qi, qi_prev, tj| dsm_lookup_raw(qi_prev, qi, tj, GAP),
                    |qi, qi_prev| dsm_lookup_raw(qi_prev, qi, GAP, GAP),
                    |qi, tj, tj_prev| dsm_lookup_raw(qi, GAP, tj_prev, tj),
                    |tj, tj_prev| dsm_lookup_raw(GAP, GAP, tj_prev, tj)
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
        traceback(
            view,
            m,
            bq,
            bt,
            best_i,
            best_j,
            &mut self.trace_buf,
        );

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

#[cfg_attr(feature = "prof", inline(never))]
fn traceback(
    view: &DpView<'_>,
    m: &ScoreOnlyGrid,
    bq: &ScoreOnlyGrid,
    bt: &ScoreOnlyGrid,
    best_i: usize,
    best_j: usize,
    trace_buf: &mut TraceVec,
) {
    let (mut i, mut j) = (best_i, best_j);
    let mut state = DpOp::Match;

    while i > 0 || j > 0 {
        match state {
            DpOp::Stop => break,
            DpOp::Match if i > 0 && j > 0 => {
                trace_buf.push(DpOp::Match);
                let m_val = m.get(i, j);
                let m_diag = m.get(i - 1, j - 1);
                let bq_diag = bq.get(i - 1, j - 1);
                let bt_diag = bt.get(i - 1, j - 1);

                // Check which transition produced this M value
                let match_e = view.match_e(i, j);
                let m_from_bq = view.m_from_bq(i, j);
                let m_from_bt = view.m_from_bt(i, j);

                let next = if m_diag > MIN_SCORE && m_val == m_diag + match_e {
                    DpOp::Match
                } else if bq_diag > MIN_SCORE && m_val == bq_diag + m_from_bq {
                    DpOp::GapQ
                } else if bt_diag > MIN_SCORE && m_val == bt_diag + m_from_bt {
                    DpOp::GapT
                } else {
                    DpOp::Stop
                };

                trace!("{} TB Match({},{}): next={:?}", view.dir, i, j, next);
                i -= 1;
                j -= 1;
                state = next;
            }
            DpOp::GapQ if i > 0 => {
                trace_buf.push(DpOp::GapQ);
                let bq_val = bq.get(i, j);
                let m_up = m.get(i - 1, j);
                let bq_up = bq.get(i - 1, j);

                let bq_open = view.bq_open(i, j);
                let bq_ext = view.bq_ext(i);

                let next = if m_up > MIN_SCORE && bq_val == m_up + bq_open {
                    DpOp::Match
                } else if bq_up > MIN_SCORE && bq_val == bq_up + bq_ext {
                    DpOp::GapQ
                } else {
                    DpOp::Stop
                };

                trace!("{} TB GapQ({},{}): next={:?}", view.dir, i, j, next);
                i -= 1;
                state = next;
            }
            DpOp::GapT if j > 0 => {
                trace_buf.push(DpOp::GapT);
                let bt_val = bt.get(i, j);
                let m_left = m.get(i, j - 1);
                let bt_left = bt.get(i, j - 1);

                let bt_open = view.bt_open(i, j);
                let bt_ext = view.bt_ext(j);

                let next = if m_left > MIN_SCORE && bt_val == m_left + bt_open {
                    DpOp::Match
                } else if bt_left > MIN_SCORE && bt_val == bt_left + bt_ext {
                    DpOp::GapT
                } else {
                    DpOp::Stop
                };

                trace!("{} TB GapT({},{}): next={:?}", view.dir, i, j, next);
                j -= 1;
                state = next;
            }
            _ => break,
        }
    }
}

impl Default for DpExtender {
    fn default() -> Self {
        Self::new()
    }
}
