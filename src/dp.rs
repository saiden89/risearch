//! Clean DP extension implementation from first principles.
//!
//! Uses nearest-neighbor stacking energy (DSM tables) with affine gap penalties.
//! Based on rust-bio patterns but adapted for RNA duplex alignment.
//!
//! # Architecture
//!
//! The key abstraction is `DpView` which encapsulates direction-specific index
//! mapping and DSM stacking order. This allows a single DP implementation to
//! handle both left and right extensions.

use crate::dsm::EnergyModel;
use crate::seq::Seq;
use crate::types::Base;
use log::trace;

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

    /// Query gap open: entering Bq state from M at (i-1, j)
    #[inline(always)]
    pub fn gap_q_open(&self, i: usize, j: usize) -> i32 {
        self.e(self.q(i - 1), self.q(i), GAP, self.t(j))
    }

    /// Query gap extend: staying in Bq state
    #[inline(always)]
    pub fn gap_q_ext(&self, i: usize) -> i32 {
        self.e(self.q(i - 1), self.q(i), GAP, GAP)
    }

    /// Target gap open: entering Bt state from M at (i, j-1)
    #[inline(always)]
    pub fn gap_t_open(&self, i: usize, j: usize) -> i32 {
        self.e(GAP, self.q(i), self.t(j - 1), self.t(j))
    }

    /// Target gap extend: staying in Bt state
    #[inline(always)]
    pub fn gap_t_ext(&self, j: usize) -> i32 {
        self.e(GAP, GAP, self.t(j - 1), self.t(j))
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
    pub trace: Vec<DpOp>,
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

// =============================================================================
// GRID TYPES - Reusable DP matrix storage
// =============================================================================

/// 2D grid for DP matrix storage
#[derive(Clone, Debug)]
pub struct Grid<T> {
    data: Vec<T>,
    width: usize,
    height: usize,
}

impl<T: Clone + Copy + Default> Grid<T> {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            data: vec![T::default(); width * height],
            width,
            height,
        }
    }

    #[inline(always)]
    pub fn get(&self, r: usize, c: usize) -> T {
        self.data[r * self.width + c]
    }

    #[inline(always)]
    pub fn set(&mut self, r: usize, c: usize, val: T) {
        self.data[r * self.width + c] = val;
    }

    /// Access diagonal neighbor (i-1, j-1)
    #[inline(always)]
    pub fn diag(&self, i: usize, j: usize) -> T {
        self.get(i - 1, j - 1)
    }

    /// Access upper neighbor (i-1, j) - Gap in Target
    #[inline(always)]
    pub fn up(&self, i: usize, j: usize) -> T {
        self.get(i - 1, j)
    }

    /// Access left neighbor (i, j-1) - Gap in Query
    #[inline(always)]
    pub fn left(&self, i: usize, j: usize) -> T {
        self.get(i, j - 1)
    }

    pub fn resize(&mut self, width: usize, height: usize) {
        let new_len = width * height;
        self.data.clear();
        self.data.resize(new_len, T::default());
        self.width = width;
        self.height = height;
    }
}

/// Score grid using Option<i32> instead of sentinel values
pub type ScoreGrid = Grid<Option<i32>>;

// =============================================================================
// DP STATE GRID - Paired score + traceback
// =============================================================================

/// Paired score and traceback grid for a single DP state (M, Bq, or Bt).
pub struct DpStateGrid {
    score: ScoreGrid,
    tb: Grid<DpOp>,
}

impl DpStateGrid {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            score: ScoreGrid::new(width, height),
            tb: Grid::new(width, height),
        }
    }

    pub fn resize(&mut self, width: usize, height: usize) {
        self.score.resize(width, height);
        self.tb.resize(width, height);
    }

    #[inline(always)]
    pub fn get(&self, i: usize, j: usize) -> Option<i32> {
        self.score.get(i, j)
    }

    #[inline(always)]
    pub fn set(&mut self, i: usize, j: usize, val: Option<i32>, step: DpOp) {
        self.score.set(i, j, val);
        self.tb.set(i, j, step);
    }

    #[inline(always)]
    pub fn diag(&self, i: usize, j: usize) -> Option<i32> {
        self.score.diag(i, j)
    }

    #[inline(always)]
    pub fn up(&self, i: usize, j: usize) -> Option<i32> {
        self.score.up(i, j)
    }

    #[inline(always)]
    pub fn left(&self, i: usize, j: usize) -> Option<i32> {
        self.score.left(i, j)
    }

    #[inline(always)]
    pub fn tb(&self, i: usize, j: usize) -> DpOp {
        self.tb.get(i, j)
    }
}

// =============================================================================
// DP EXTENDER - Stateful extension with reusable matrices
// =============================================================================

/// DP matrices for extension (reusable to avoid allocations)
pub struct DpMatrices {
    pub m: DpStateGrid,  // Match/mismatch state
    pub bq: DpStateGrid, // Query bulge (gap in target)
    pub bt: DpStateGrid, // Target bulge (gap in query)
}

impl DpMatrices {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            m: DpStateGrid::new(width, height),
            bq: DpStateGrid::new(width, height),
            bt: DpStateGrid::new(width, height),
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
}

// =============================================================================
// INIT HELPERS - Reduce code duplication in DP initialization
// =============================================================================

/// Pick best value from two optional sources with tie-breaking.
/// Returns (best_value, which_source).
/// When both have equal value, prefers `a` (first argument).
#[inline(always)]
fn pick_best(a: Option<i32>, b: Option<i32>, step_a: DpOp, step_b: DpOp) -> (Option<i32>, DpOp) {
    match (a, b) {
        (Some(va), Some(vb)) if va >= vb => (Some(va), step_a),
        (Some(_), Some(vb)) => (Some(vb), step_b),
        (Some(va), None) => (Some(va), step_a),
        (None, Some(vb)) => (Some(vb), step_b),
        (None, None) => (None, DpOp::Stop),
    }
}

/// Initialize limited rows or columns - unified via macro.
/// Row axis: Bt primary, Bq secondary, indices (fixed, k)
/// Col axis: Bq primary, Bt secondary, indices (k, fixed)
macro_rules! init_limited_axis {
    ($len:expr, $view:expr, $m:expr, $update_best:expr,
     $primary:expr, $secondary:expr, $gap_op:expr,
     $idx:expr, $prev_idx:expr,
     $b_open:expr, $b_ext:expr, $match_e:expr, $m_from_b:expr, $s_open:expr) => {
        for k in 3..$len {
            let (pi1, pj1) = $prev_idx(1, k);
            let (i1, j1) = $idx(1, k);
            let (i2, j2) = $idx(2, k);
            let (pi2, pj2) = $prev_idx(2, k);
            let (i2p, j2p) = $idx(2, k - 1);

            // Primary[1,k]
            let from_m = $m.get(pi1, pj1).map(|v| v + $b_open($view, 1, k));
            let from_b = $primary.get(pi1, pj1).map(|v| v + $b_ext($view, k));
            let (val, step) = pick_best(from_m, from_b, DpOp::Match, $gap_op);
            $primary.set(i1, j1, val, step);

            // M[2,k]
            let from_m = $m.get(pi1, pj1).map(|v| v + $match_e($view, 2, k));
            let from_b = $primary.get(pi1, pj1).map(|v| v + $m_from_b($view, 2, k));
            let (val, step) = pick_best(from_m, from_b, DpOp::Match, $gap_op);
            $m.set(i2, j2, val, step);
            $update_best(val, i2, j2);

            // Secondary[2,k]
            if let Some(m1k) = $m.get(i1, j1) {
                $secondary.set(i2, j2, Some(m1k + $s_open($view, 2, k)), DpOp::Match);
            }

            // Primary[2,k]
            let from_m = $m.get(pi2, pj2).map(|v| v + $b_open($view, 2, k));
            let from_b = $primary.get(i2p, j2p).map(|v| v + $b_ext($view, k));
            let (val, step) = pick_best(from_m, from_b, DpOp::Match, $gap_op);
            $primary.set(i2, j2, val, step);
        }
    };
}

impl DpExtender {
    pub fn new() -> Self {
        Self {
            matrices: DpMatrices::new(200, 200),
        }
    }
    /// Extend to the left (query 5', target 3')
    ///
    /// Thin wrapper around `extend()` - creates a LEFT view and delegates.
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
    ///
    /// Thin wrapper around `extend()` - creates a RIGHT view and delegates.
    ///
    /// NOTE: Both `extend_left` and `extend_right` exist for API compatibility.
    /// Consider using `extend(&DpView)` directly for new code.
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

    /// Unified DP extension using DpView for direction-aware energy calculations.
    ///
    /// This method handles both left and right extensions through the DpView
    /// abstraction, which encapsulates:
    /// - Index mapping (view.q(i), view.t(j))
    /// - DSM stacking order (view.e() swaps for LEFT extension)
    /// - Terminal energy (view.terminal() handles direction-specific gap position)
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
                trace: Vec::new(),
            };
        }

        // Resize matrices
        self.matrices.resize(t_len + 1, q_len + 1);
        let DpMatrices { m, bq, bt } = &mut self.matrices;

        // =====================================================================
        // HELPER: Update best score if value + terminal is better
        // =====================================================================
        let mut update_best = |val: Option<i32>, i: usize, j: usize| {
            if let Some(v) = val {
                let curr = v + view.terminal(i, j);
                if curr > best_e {
                    best_e = curr;
                    best_i = i;
                    best_j = j;
                }
            }
        };

        // =====================================================================
        // INITIALIZATION using view helper methods
        // =====================================================================

        // Corner: M[0,0], Bt[0,1], Bq[1,0], M[1,1]
        m.set(0, 0, Some(0), DpOp::Stop);
        bt.set(0, 1, Some(view.bt_open(0, 1)), DpOp::Stop);
        bq.set(1, 0, Some(view.bq_open(1, 0)), DpOp::Stop);
        let m11 = Some(view.match_e(1, 1));
        m.set(1, 1, m11, DpOp::Stop);
        update_best(m11, 1, 1);

        // Boundary init: Row 0 (Bt) and Col 0 (Bq) are symmetric
        for k in 2..t_len {
            if let Some(prev) = bt.left(0, k) {
                bt.set(0, k, Some(prev + view.bt_ext(k)), DpOp::GapT);
                let new_m = Some(prev + view.m_from_bt(1, k));
                m.set(1, k, new_m, DpOp::GapT);
                update_best(new_m, 1, k);
            }
        }
        for k in 2..q_len {
            if let Some(prev) = bq.up(k, 0) {
                bq.set(k, 0, Some(prev + view.bq_ext(k)), DpOp::GapQ);
                let new_m = Some(prev + view.m_from_bq(k, 1));
                m.set(k, 1, new_m, DpOp::GapQ);
                update_best(new_m, k, 1);
            }
        }

        // Cell (2,2) init: bridge corner to limited rows/cols
        if q_len >= 3 && t_len >= 3 {
            if let Some(m11_val) = m.diag(2, 2) {
                bt.set(1, 2, Some(m11_val + view.bt_open(1, 2)), DpOp::Match);
                bq.set(2, 1, Some(m11_val + view.bq_open(2, 1)), DpOp::Match);
                let m22 = Some(m11_val + view.match_e(2, 2));
                m.set(2, 2, m22, DpOp::Match);
                update_best(m22, 2, 2);
            }
            if let Some(m12) = m.up(2, 2) {
                bq.set(2, 2, Some(m12 + view.bq_open(2, 2)), DpOp::Match);
            }
            if let Some(m21) = m.left(2, 2) {
                bt.set(2, 2, Some(m21 + view.bt_open(2, 2)), DpOp::Match);
            }
        }

        // =======================================================================
        // LIMITED ROWS/COLUMNS - unified via macro
        // =======================================================================

        init_limited_axis!(
            t_len,
            view,
            m,
            update_best,
            bt,
            bq,
            DpOp::GapT,
            |f, k| (f, k),
            |f, k| (f, k - 1),
            |v: &DpView, f, k| v.bt_open(f, k),
            |v: &DpView, k| v.bt_ext(k),
            |v: &DpView, f, k| v.match_e(f, k),
            |v: &DpView, f, k| v.m_from_bt(f, k),
            |v: &DpView, f, k| v.bq_open(f, k)
        );

        init_limited_axis!(
            q_len,
            view,
            m,
            update_best,
            bq,
            bt,
            DpOp::GapQ,
            |f, k| (k, f),
            |f, k| (k - 1, f),
            |v: &DpView, f, k| v.bq_open(k, f),
            |v: &DpView, k| v.bq_ext(k),
            |v: &DpView, f, k| v.match_e(k, f),
            |v: &DpView, f, k| v.m_from_bq(k, f),
            |v: &DpView, f, k| v.bt_open(k, f)
        );

        // =======================================================================
        // MAIN DP LOOP (i >= 3, j >= 3)
        // =======================================================================

        for i in 3..q_len {
            for j in 3..t_len {
                // M[i,j] - pick best of three sources
                let s_mm = m.diag(i, j).map(|v| v + view.match_e(i, j));
                let s_mq = bq.diag(i, j).map(|v| v + view.m_from_bq(i, j));
                let s_mt = bt.diag(i, j).map(|v| v + view.m_from_bt(i, j));

                // Tie-breaking: Match > GapT > GapQ (C's max3 priority)
                let (val_m, step_m) = [(s_mq, DpOp::GapQ), (s_mt, DpOp::GapT), (s_mm, DpOp::Match)]
                    .into_iter()
                    .filter_map(|(opt, step)| opt.map(|v| (v, step)))
                    .max_by_key(|(v, _)| *v)
                    .unwrap_or((i32::MIN, DpOp::Stop));

                let val_m = if val_m == i32::MIN { None } else { Some(val_m) };
                m.set(i, j, val_m, step_m);
                update_best(val_m, i, j);

                // Bq[i,j] - query bulge state (gap in target)
                let s_qm = m.up(i, j).map(|v| v + view.bq_open(i, j));
                let s_qq = bq.up(i, j).map(|v| v + view.bq_ext(i));
                let (bq_val, step) = pick_best(s_qm, s_qq, DpOp::Match, DpOp::GapQ);
                bq.set(i, j, bq_val, step);

                // Bt[i,j] - target bulge state (gap in query)
                let s_tm = m.left(i, j).map(|v| v + view.bt_open(i, j));
                let s_tt = bt.left(i, j).map(|v| v + view.bt_ext(j));
                let (bt_val, step) = pick_best(s_tm, s_tt, DpOp::Match, DpOp::GapT);
                bt.set(i, j, bt_val, step);
            }
        }

        // =======================================================================
        // TRACEBACK
        // =======================================================================

        trace!(
            "{} TB start: best=({},{}) score={}",
            view.dir, best_i, best_j, best_e
        );
        let mut trace_vec = Vec::new();
        let (mut i, mut j) = (best_i, best_j);
        let mut state = DpOp::Match;

        while i > 0 || j > 0 {
            let (next, di, dj) = match state {
                DpOp::Stop => break,
                DpOp::Match if i > 0 && j > 0 => (m.tb(i, j), 1, 1),
                DpOp::GapQ if i > 0 => (bq.tb(i, j), 1, 0),
                DpOp::GapT if j > 0 => (bt.tb(i, j), 0, 1),
                _ => break,
            };
            trace_vec.push(state);
            trace!("{} TB {:?}({},{}): next={:?}", view.dir, state, i, j, next);
            i -= di;
            j -= dj;
            if next == DpOp::Stop {
                break;
            }
            state = next;
        }

        trace!(
            "{} result: score={} q_len={} t_len={} trace={:?}",
            view.dir, best_e, best_i, best_j, trace_vec
        );
        DpExtension {
            score: best_e,
            q_len: best_i,
            t_len: best_j,
            trace: trace_vec,
        }
    }
}

impl Default for DpExtender {
    fn default() -> Self {
        Self::new()
    }
}
