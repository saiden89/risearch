//! Clean DP extension implementation from first principles.
//!
//! Uses nearest-neighbor stacking energy (DSM tables) with affine gap penalties.
//! Based on rust-bio patterns but adapted for RNA duplex alignment.

use crate::dsm::{EnergyModel, StackPair};
use crate::seq::Seq;
use crate::types::Base;

/// Extension direction - determines terminal stacking order
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ExtendDir {
    Left,  // Terminal: Gap→Q, Gap→T (extending into sequence from gap)
    Right, // Terminal: Q→Gap, T→Gap (extending out of sequence into gap)
}

/// Extend alignment to the left (query 5', target 3')
///
/// Query extends toward 5' (decreasing index), Target extends toward 3' (increasing index).
pub fn extend_left(
    query: &Seq,
    target: &Seq,
    q_start: usize,
    t_start: usize,
    max_ext: usize,
) -> DpExtension {
    let q_len = (q_start + 1).min(max_ext);
    let t_len = (target.len() - t_start - 1).min(max_ext);

    extend(
        |i| query.left(q_start, i),
        |j| target.right(t_start, j),
        q_len,
        t_len,
        ExtendDir::Left,
    )
}

/// Extend alignment to the right (query 3', target 5')
///
/// Query extends toward 3' (increasing index), Target extends toward 5' (decreasing index).
pub fn extend_right(
    query: &Seq,
    target: &Seq,
    q_end: usize,
    t_end: usize,
    max_ext: usize,
) -> DpExtension {
    let q_len = (query.len() - q_end).min(max_ext);
    let t_len = (t_end + 1).min(max_ext);

    extend(
        |i| query.right(q_end, i),
        |j| target.left(t_end, j),
        q_len,
        t_len,
        ExtendDir::Right,
    )
}

/// Result of DP extension
#[derive(Debug, Clone)]
pub struct DpExtension {
    pub score: i32,
    pub q_len: usize,
    pub t_len: usize,
    pub trace: Vec<DpOp>,
}

// =============================================================================
// EXTENDER TRAIT - Allows swappable DP implementations
// =============================================================================

/// Trait for DP extension implementations.
///
/// Allows different algorithms to be used interchangeably:
/// - `LegacyExtender`: Current algorithm with C parity
/// - `OptimizedExtender`: Future 2-row optimized version
pub trait Extender {
    /// Extend alignment to the left (query 5', target 3')
    fn extend_left(
        &mut self,
        query: &Seq,
        target: &Seq,
        q_start: usize,
        t_start: usize,
        max_ext: usize,
    ) -> DpExtension;

    /// Extend alignment to the right (query 3', target 5')
    fn extend_right(
        &mut self,
        query: &Seq,
        target: &Seq,
        q_end: usize,
        t_end: usize,
        max_ext: usize,
    ) -> DpExtension;
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

/// DP state for traceback navigation
#[derive(Debug, Clone, Copy)]
enum DpState {
    Match,
    GapQ,
    GapT,
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
// DP EXTENDER - Stateful extension with reusable matrices
// =============================================================================

/// DP matrices for extension (reusable to avoid allocations)
pub struct DpMatrices {
    // Score matrices (None = not reachable)
    pub m: ScoreGrid,  // Match/mismatch state
    pub bq: ScoreGrid, // Query bulge (gap in target)
    pub bt: ScoreGrid, // Target bulge (gap in query)

    // Traceback matrices
    pub tb_m: Grid<DpOp>,
    pub tb_bq: Grid<DpOp>,
    pub tb_bt: Grid<DpOp>,
}

impl DpMatrices {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            m: ScoreGrid::new(width, height),
            bq: ScoreGrid::new(width, height),
            bt: ScoreGrid::new(width, height),
            tb_m: Grid::new(width, height),
            tb_bq: Grid::new(width, height),
            tb_bt: Grid::new(width, height),
        }
    }

    pub fn resize(&mut self, width: usize, height: usize) {
        self.m.resize(width, height);
        self.bq.resize(width, height);
        self.bt.resize(width, height);
        self.tb_m.resize(width, height);
        self.tb_bq.resize(width, height);
        self.tb_bt.resize(width, height);
    }
}

/// Stateful DP extender with reusable matrices
pub struct DpExtender {
    matrices: DpMatrices,
}

impl DpExtender {
    pub fn new() -> Self {
        Self {
            matrices: DpMatrices::new(200, 200),
        }
    }

    /// Extend to the left (query 5', target 3') - Legacy algorithm with C parity
    pub fn extend_left(
        &mut self,
        query: &Seq,
        target: &Seq,
        q_start: usize,
        t_start: usize,
        max_ext: usize,
    ) -> DpExtension {
        // DP Left: Extend Query toward 5' (decreasing index), Target toward 3' (increasing index)
        let q_len = (q_start + 1).min(max_ext);
        let t_len = (target.len() - t_start - 1).min(max_ext);

        // Energy lookup using centralized EnergyModel
        #[inline(always)]
        fn e(q1: usize, q2: usize, t1: usize, t2: usize) -> i32 {
            EnergyModel::T04.stack_idx(q1, q2, t1, t2)
        }
        const GAP: usize = Base::Gap as usize;

        // Base index accessors
        let q_idx = |i: usize| -> usize { query.left(q_start, i).idx() };
        let t_idx = |j: usize| -> usize { target.right(t_start, j).idx() };

        // Initial score: terminal penalty for seed boundary
        let mut best_e = e(GAP, q_idx(0), GAP, t_idx(0));
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
        let DpMatrices {
            m,
            bq,
            bt,
            tb_m,
            tb_bq,
            tb_bt,
        } = &mut self.matrices;

        // Init M[0,0]
        m.set(0, 0, Some(0));

        // Init (0,1), (1,0), (1,1)
        bt.set(0, 1, Some(e(GAP, q_idx(0), t_idx(1), t_idx(0))));
        bq.set(1, 0, Some(e(q_idx(1), q_idx(0), GAP, t_idx(0))));
        m.set(1, 1, Some(e(q_idx(1), q_idx(0), t_idx(1), t_idx(0))));

        if let Some(m11) = m.get(1, 1) {
            let val = m11 + e(GAP, q_idx(1), GAP, t_idx(1));
            if val > best_e {
                best_e = val;
                best_i = 1;
                best_j = 1;
            }
        }

        // Row 0 (j from 2 to t_len-1)
        for j in 2..t_len {
            if let Some(prev_bt) = bt.left(0, j) {
                bt.set(0, j, Some(prev_bt + e(GAP, GAP, t_idx(j), t_idx(j - 1))));
                tb_bt.set(0, j, DpOp::GapT);

                if q_len >= 1 {
                    let new_m = prev_bt + e(q_idx(1), GAP, t_idx(j), t_idx(j - 1));
                    m.set(1, j, Some(new_m));
                    tb_m.set(1, j, DpOp::GapT);
                    let val = new_m + e(GAP, q_idx(1), GAP, t_idx(j));
                    if val > best_e {
                        best_e = val;
                        best_i = 1;
                        best_j = j;
                    }
                }
            }
        }

        // Col 0 (i from 2 to q_len-1)
        for i in 2..q_len {
            if let Some(prev_bq) = bq.up(i, 0) {
                bq.set(i, 0, Some(prev_bq + e(q_idx(i), q_idx(i - 1), GAP, GAP)));
                tb_bq.set(i, 0, DpOp::GapQ);

                if t_len >= 1 {
                    let new_m = prev_bq + e(q_idx(i), q_idx(i - 1), t_idx(1), GAP);
                    m.set(i, 1, Some(new_m));
                    tb_m.set(i, 1, DpOp::GapQ);
                    let val = new_m + e(GAP, q_idx(i), GAP, t_idx(1));
                    if val > best_e {
                        best_e = val;
                        best_i = i;
                        best_j = 1;
                    }
                }
            }
        }

        // 2,2 Init
        if q_len >= 3 && t_len >= 3 {
            if let Some(m11) = m.diag(2, 2) {
                bt.set(1, 2, Some(m11 + e(GAP, q_idx(1), t_idx(2), t_idx(1))));
                tb_bt.set(1, 2, DpOp::Match);

                bq.set(2, 1, Some(m11 + e(q_idx(2), q_idx(1), GAP, t_idx(1))));
                tb_bq.set(2, 1, DpOp::Match);

                let m22 = m11 + e(q_idx(2), q_idx(1), t_idx(2), t_idx(1));
                m.set(2, 2, Some(m22));
                tb_m.set(2, 2, DpOp::Match);

                let val = m22 + e(GAP, q_idx(2), GAP, t_idx(2));
                if val > best_e {
                    best_e = val;
                    best_i = 2;
                    best_j = 2;
                }
            }
            if let Some(m12) = m.up(2, 2) {
                bq.set(2, 2, Some(m12 + e(q_idx(2), q_idx(1), GAP, t_idx(2))));
                tb_bq.set(2, 2, DpOp::Match);
            }
            if let Some(m21) = m.left(2, 2) {
                bt.set(2, 2, Some(m21 + e(GAP, q_idx(2), t_idx(2), t_idx(1))));
                tb_bt.set(2, 2, DpOp::Match);
            }
        }

        // Main DP loop
        for i in 2..q_len {
            for j in 2..t_len {
                if i == 2 && j == 2 {
                    continue;
                }

                // M[i,j] - pick best of three sources
                let s_mm = m
                    .diag(i, j)
                    .map(|v| v + e(q_idx(i), q_idx(i - 1), t_idx(j), t_idx(j - 1)));
                let s_mq = bq
                    .diag(i, j)
                    .map(|v| v + e(q_idx(i), q_idx(i - 1), t_idx(j), GAP));
                let s_mt = bt
                    .diag(i, j)
                    .map(|v| v + e(q_idx(i), GAP, t_idx(j), t_idx(j - 1)));

                // Tie-breaking: GapQ wins
                let (val_m, step_m) = [(s_mt, DpOp::GapT), (s_mm, DpOp::Match), (s_mq, DpOp::GapQ)]
                    .into_iter()
                    .filter_map(|(opt, step)| opt.map(|v| (v, step)))
                    .max_by_key(|(v, _)| *v)
                    .unwrap_or((i32::MIN, DpOp::Stop));

                let val_m = if val_m == i32::MIN { None } else { Some(val_m) };
                m.set(i, j, val_m);
                tb_m.set(i, j, step_m);

                if let Some(v) = val_m {
                    let curr_e = v + e(GAP, q_idx(i), GAP, t_idx(j));
                    if curr_e > best_e {
                        best_e = curr_e;
                        best_i = i;
                        best_j = j;
                    }
                }

                // Bq[i,j]
                if i > 2 || (i == 2 && j > 2) {
                    let s_qm = m
                        .up(i, j)
                        .map(|v| v + e(q_idx(i), q_idx(i - 1), GAP, t_idx(j)));
                    let s_qq = bq.up(i, j).map(|v| v + e(q_idx(i), q_idx(i - 1), GAP, GAP));
                    match (s_qq, s_qm) {
                        (Some(qq), Some(qm)) if qq > qm => {
                            bq.set(i, j, Some(qq));
                            tb_bq.set(i, j, DpOp::GapQ);
                        }
                        (_, Some(qm)) => {
                            bq.set(i, j, Some(qm));
                            tb_bq.set(i, j, DpOp::Match);
                        }
                        (Some(qq), None) => {
                            bq.set(i, j, Some(qq));
                            tb_bq.set(i, j, DpOp::GapQ);
                        }
                        _ => {}
                    }
                }

                // Bt[i,j]
                if j > 2 || (j == 2 && i > 2) {
                    let s_tm = m
                        .left(i, j)
                        .map(|v| v + e(GAP, q_idx(i), t_idx(j), t_idx(j - 1)));
                    let s_tt = bt
                        .left(i, j)
                        .map(|v| v + e(GAP, GAP, t_idx(j), t_idx(j - 1)));
                    match (s_tt, s_tm) {
                        (Some(tt), Some(tm)) if tt > tm => {
                            bt.set(i, j, Some(tt));
                            tb_bt.set(i, j, DpOp::GapT);
                        }
                        (_, Some(tm)) => {
                            bt.set(i, j, Some(tm));
                            tb_bt.set(i, j, DpOp::Match);
                        }
                        (Some(tt), None) => {
                            bt.set(i, j, Some(tt));
                            tb_bt.set(i, j, DpOp::GapT);
                        }
                        _ => {}
                    }
                }
            }
        }

        // Traceback
        let mut trace_vec = Vec::new();
        let (mut i, mut j) = (best_i, best_j);
        let mut state = DpState::Match;

        while i > 0 || j > 0 {
            match state {
                DpState::Match => {
                    if i == 0 || j == 0 {
                        break;
                    }
                    let step = tb_m.get(i, j);
                    i -= 1;
                    j -= 1;
                    trace_vec.push(step);
                    state = match step {
                        DpOp::Stop => break,
                        DpOp::Match => DpState::Match,
                        DpOp::GapQ => DpState::GapQ,
                        DpOp::GapT => DpState::GapT,
                    };
                }
                DpState::GapQ => {
                    let step = tb_bq.get(i, j);
                    trace_vec.push(step);
                    if i > 0 {
                        i -= 1;
                    } else {
                        break;
                    }
                    state = match step {
                        DpOp::Stop => break,
                        DpOp::Match => DpState::Match,
                        DpOp::GapQ => DpState::GapQ,
                        DpOp::GapT => DpState::GapT,
                    };
                }
                DpState::GapT => {
                    let step = tb_bt.get(i, j);
                    trace_vec.push(step);
                    if j > 0 {
                        j -= 1;
                    } else {
                        break;
                    }
                    state = match step {
                        DpOp::Stop => break,
                        DpOp::Match => DpState::Match,
                        DpOp::GapQ => DpState::GapQ,
                        DpOp::GapT => DpState::GapT,
                    };
                }
            }
        }

        DpExtension {
            score: best_e,
            q_len: best_i,
            t_len: best_j,
            trace: trace_vec,
        }
    }

    /// Extend to the right (query 3', target 5') - Legacy algorithm with C parity
    pub fn extend_right(
        &mut self,
        query: &Seq,
        target: &Seq,
        q_end: usize,
        t_end: usize,
        max_ext: usize,
    ) -> DpExtension {
        // DP Right: Extend Query toward 3' (increasing), Target toward 5' (decreasing)
        let q_len = (query.len() - q_end).min(max_ext);
        let t_len = (t_end + 1).min(max_ext);

        // Energy lookup
        #[inline(always)]
        fn e(q1: usize, q2: usize, t1: usize, t2: usize) -> i32 {
            EnergyModel::T04.stack_idx(q1, q2, t1, t2)
        }
        const GAP: usize = Base::Gap as usize;

        // Base index accessors (reversed direction)
        let q_idx = |i: usize| -> usize { query.right(q_end, i).idx() };
        let t_idx = |j: usize| -> usize { target.left(t_end, j).idx() };

        // Initial score: terminal penalty (stacking order for right extension)
        let mut best_e = e(q_idx(0), GAP, t_idx(0), GAP);
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
        let DpMatrices {
            m,
            bq,
            bt,
            tb_m,
            tb_bq,
            tb_bt,
        } = &mut self.matrices;

        // Init M[0,0]
        m.set(0, 0, Some(0));

        // Init (0,1), (1,0), (1,1) - note reversed stacking order for right extension
        bt.set(0, 1, Some(e(q_idx(0), GAP, t_idx(0), t_idx(1))));
        bq.set(1, 0, Some(e(q_idx(0), q_idx(1), t_idx(0), GAP)));
        m.set(1, 1, Some(e(q_idx(0), q_idx(1), t_idx(0), t_idx(1))));

        if let Some(m11) = m.get(1, 1) {
            let val = m11 + e(q_idx(1), GAP, t_idx(1), GAP);
            if val > best_e {
                best_e = val;
                best_i = 1;
                best_j = 1;
            }
        }

        // Row 0 (j from 2 to t_len-1)
        for j in 2..t_len {
            if let Some(prev_bt) = bt.left(0, j) {
                bt.set(0, j, Some(prev_bt + e(GAP, GAP, t_idx(j - 1), t_idx(j))));
                tb_bt.set(0, j, DpOp::GapT);

                if q_len >= 1 {
                    let new_m = prev_bt + e(GAP, q_idx(1), t_idx(j - 1), t_idx(j));
                    m.set(1, j, Some(new_m));
                    tb_m.set(1, j, DpOp::GapT);
                    let val = new_m + e(q_idx(1), GAP, t_idx(j), GAP);
                    if val > best_e {
                        best_e = val;
                        best_i = 1;
                        best_j = j;
                    }
                }
            }
        }

        // Col 0 (i from 2 to q_len-1)
        for i in 2..q_len {
            if let Some(prev_bq) = bq.up(i, 0) {
                bq.set(i, 0, Some(prev_bq + e(q_idx(i - 1), q_idx(i), GAP, GAP)));
                tb_bq.set(i, 0, DpOp::GapQ);

                if t_len >= 1 {
                    let new_m = prev_bq + e(q_idx(i - 1), q_idx(i), GAP, t_idx(1));
                    m.set(i, 1, Some(new_m));
                    tb_m.set(i, 1, DpOp::GapQ);
                    let val = new_m + e(q_idx(i), GAP, t_idx(1), GAP);
                    if val > best_e {
                        best_e = val;
                        best_i = i;
                        best_j = 1;
                    }
                }
            }
        }

        // 2,2 Init
        if q_len >= 3 && t_len >= 3 {
            if let Some(m11) = m.diag(2, 2) {
                bt.set(1, 2, Some(m11 + e(q_idx(1), GAP, t_idx(1), t_idx(2))));
                tb_bt.set(1, 2, DpOp::Match);

                bq.set(2, 1, Some(m11 + e(q_idx(1), q_idx(2), t_idx(1), GAP)));
                tb_bq.set(2, 1, DpOp::Match);

                let m22 = m11 + e(q_idx(1), q_idx(2), t_idx(1), t_idx(2));
                m.set(2, 2, Some(m22));
                tb_m.set(2, 2, DpOp::Match);

                let val = m22 + e(q_idx(2), GAP, t_idx(2), GAP);
                if val > best_e {
                    best_e = val;
                    best_i = 2;
                    best_j = 2;
                }
            }
            if let Some(m12) = m.up(2, 2) {
                bq.set(2, 2, Some(m12 + e(q_idx(1), q_idx(2), t_idx(2), GAP)));
                tb_bq.set(2, 2, DpOp::Match);
            }
            if let Some(m21) = m.left(2, 2) {
                bt.set(2, 2, Some(m21 + e(q_idx(2), GAP, t_idx(1), t_idx(2))));
                tb_bt.set(2, 2, DpOp::Match);
            }
        }

        // Main DP loop
        for i in 2..q_len {
            for j in 2..t_len {
                if i == 2 && j == 2 {
                    continue;
                }

                // M[i,j] - reversed stacking order for right extension
                let s_mm = m
                    .diag(i, j)
                    .map(|v| v + e(q_idx(i - 1), q_idx(i), t_idx(j - 1), t_idx(j)));
                let s_mq = bq
                    .diag(i, j)
                    .map(|v| v + e(q_idx(i - 1), q_idx(i), GAP, t_idx(j)));
                let s_mt = bt
                    .diag(i, j)
                    .map(|v| v + e(GAP, q_idx(i), t_idx(j - 1), t_idx(j)));

                // Tie-breaking: GapQ wins
                let (val_m, step_m) = [(s_mt, DpOp::GapT), (s_mm, DpOp::Match), (s_mq, DpOp::GapQ)]
                    .into_iter()
                    .filter_map(|(opt, step)| opt.map(|v| (v, step)))
                    .max_by_key(|(v, _)| *v)
                    .unwrap_or((i32::MIN, DpOp::Stop));

                let val_m = if val_m == i32::MIN { None } else { Some(val_m) };
                m.set(i, j, val_m);
                tb_m.set(i, j, step_m);

                if let Some(v) = val_m {
                    let curr_e = v + e(q_idx(i), GAP, t_idx(j), GAP);
                    if curr_e > best_e {
                        best_e = curr_e;
                        best_i = i;
                        best_j = j;
                    }
                }

                // Bq[i,j]
                if i > 2 || (i == 2 && j > 2) {
                    let s_qm = m
                        .up(i, j)
                        .map(|v| v + e(q_idx(i - 1), q_idx(i), t_idx(j), GAP));
                    let s_qq = bq.up(i, j).map(|v| v + e(q_idx(i - 1), q_idx(i), GAP, GAP));
                    match (s_qq, s_qm) {
                        (Some(qq), Some(qm)) if qq > qm => {
                            bq.set(i, j, Some(qq));
                            tb_bq.set(i, j, DpOp::GapQ);
                        }
                        (_, Some(qm)) => {
                            bq.set(i, j, Some(qm));
                            tb_bq.set(i, j, DpOp::Match);
                        }
                        (Some(qq), None) => {
                            bq.set(i, j, Some(qq));
                            tb_bq.set(i, j, DpOp::GapQ);
                        }
                        _ => {}
                    }
                }

                // Bt[i,j]
                if j > 2 || (j == 2 && i > 2) {
                    let s_tm = m
                        .left(i, j)
                        .map(|v| v + e(q_idx(i), GAP, t_idx(j - 1), t_idx(j)));
                    let s_tt = bt
                        .left(i, j)
                        .map(|v| v + e(GAP, GAP, t_idx(j - 1), t_idx(j)));
                    match (s_tt, s_tm) {
                        (Some(tt), Some(tm)) if tt > tm => {
                            bt.set(i, j, Some(tt));
                            tb_bt.set(i, j, DpOp::GapT);
                        }
                        (_, Some(tm)) => {
                            bt.set(i, j, Some(tm));
                            tb_bt.set(i, j, DpOp::Match);
                        }
                        (Some(tt), None) => {
                            bt.set(i, j, Some(tt));
                            tb_bt.set(i, j, DpOp::GapT);
                        }
                        _ => {}
                    }
                }
            }
        }

        // Traceback
        let mut trace_vec = Vec::new();
        let (mut i, mut j) = (best_i, best_j);
        let mut state = DpState::Match;

        while i > 0 || j > 0 {
            match state {
                DpState::Match => {
                    if i == 0 || j == 0 {
                        break;
                    }
                    let step = tb_m.get(i, j);
                    i -= 1;
                    j -= 1;
                    trace_vec.push(step);
                    state = match step {
                        DpOp::Stop => break,
                        DpOp::Match => DpState::Match,
                        DpOp::GapQ => DpState::GapQ,
                        DpOp::GapT => DpState::GapT,
                    };
                }
                DpState::GapQ => {
                    let step = tb_bq.get(i, j);
                    trace_vec.push(step);
                    if i > 0 {
                        i -= 1;
                    } else {
                        break;
                    }
                    state = match step {
                        DpOp::Stop => break,
                        DpOp::Match => DpState::Match,
                        DpOp::GapQ => DpState::GapQ,
                        DpOp::GapT => DpState::GapT,
                    };
                }
                DpState::GapT => {
                    let step = tb_bt.get(i, j);
                    trace_vec.push(step);
                    if j > 0 {
                        j -= 1;
                    } else {
                        break;
                    }
                    state = match step {
                        DpOp::Stop => break,
                        DpOp::Match => DpState::Match,
                        DpOp::GapQ => DpState::GapQ,
                        DpOp::GapT => DpState::GapT,
                    };
                }
            }
        }

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

impl Extender for DpExtender {
    fn extend_left(
        &mut self,
        query: &Seq,
        target: &Seq,
        q_start: usize,
        t_start: usize,
        max_ext: usize,
    ) -> DpExtension {
        // Delegate to inherent method
        DpExtender::extend_left(self, query, target, q_start, t_start, max_ext)
    }

    fn extend_right(
        &mut self,
        query: &Seq,
        target: &Seq,
        q_end: usize,
        t_end: usize,
        max_ext: usize,
    ) -> DpExtension {
        // Delegate to inherent method
        DpExtender::extend_right(self, query, target, q_end, t_end, max_ext)
    }
}

/// Core DP extension function
#[allow(clippy::needless_range_loop)] // Index i is used for matrix access across multiple arrays
pub fn extend<Q, T>(q: Q, t: T, q_len: usize, t_len: usize, dir: ExtendDir) -> DpExtension
where
    Q: Fn(usize) -> Base,
    T: Fn(usize) -> Base,
{
    // Terminal stacking order: Gap→Q, Gap→T for both directions (matches C behavior)
    let terminal_fn = |q_base: Base, t_base: Base| -> i32 {
        StackPair::new(Base::Gap, q_base, Base::Gap, t_base).energy() as i32
    };

    // Early return if nothing to extend - but still return initial terminal
    if q_len == 0 || t_len == 0 {
        let initial_terminal = terminal_fn(q(0), t(0));
        return DpExtension {
            score: initial_terminal,
            q_len: 0,
            t_len: 0,
            trace: vec![],
        };
    }

    let stack = |q1: Base, q2: Base, t1: Base, t2: Base| -> i32 {
        StackPair::new(q1, q2, t1, t2).energy() as i32
    };

    // DP matrices (2-row optimization)
    let mut m_prev = vec![i32::MIN / 2; t_len + 1];
    let mut m_curr = vec![i32::MIN / 2; t_len + 1];
    let mut bq_prev = vec![i32::MIN / 2; t_len + 1];
    let mut bq_curr = vec![i32::MIN / 2; t_len + 1];
    let mut bt_prev = vec![i32::MIN / 2; t_len + 1];
    let mut bt_curr = vec![i32::MIN / 2; t_len + 1];

    // Traceback storage
    let mut tb: Vec<Vec<DpOp>> = vec![vec![DpOp::Match; t_len + 1]; q_len + 1];

    // Best score tracking - initialize with terminal penalty for zero extension
    let initial_terminal = terminal_fn(q(0), t(0));
    let mut best_score = initial_terminal;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    // Initialize M[0,0] = 0 (set in m_curr because loop swaps first)
    m_curr[0] = 0;

    for i in 1..=q_len {
        std::mem::swap(&mut m_prev, &mut m_curr);
        std::mem::swap(&mut bq_prev, &mut bq_curr);
        std::mem::swap(&mut bt_prev, &mut bt_curr);

        for val in m_curr.iter_mut() {
            *val = i32::MIN / 2;
        }
        for val in bq_curr.iter_mut() {
            *val = i32::MIN / 2;
        }
        for val in bt_curr.iter_mut() {
            *val = i32::MIN / 2;
        }

        for j in 1..=t_len {
            let s = stack(q(i - 1), q(i), t(j - 1), t(j));

            // Match: transition from M, Bq, or Bt diagonal
            let m_from_m = m_prev[j - 1] + s;
            let m_from_bq = bq_prev[j - 1] + s;
            let m_from_bt = bt_prev[j - 1] + s;

            // Find best score and traceback for M[i,j]
            // Tie-breaking matches C's max3(M, Bq, Bt): GapQ wins ties
            // Process in order: GapT, Match, GapQ - last one wins on >= so GapQ wins ties
            let (m_best, m_tb) = [
                (m_from_bt, DpOp::GapT), // lowest priority
                (m_from_m, DpOp::Match), // middle priority
                (m_from_bq, DpOp::GapQ), // highest priority (wins ties)
            ]
            .into_iter()
            .max_by_key(|(score, _)| *score)
            .unwrap_or((i32::MIN / 2, DpOp::Stop));

            m_curr[j] = m_best;
            tb[i][j] = m_tb;

            // Bq: gap in query (from above)
            let bq_open = m_prev[j] + stack(q(i - 1), q(i), Base::Gap, t(j));
            let bq_ext = bq_prev[j] + stack(q(i - 1), q(i), Base::Gap, Base::Gap);
            bq_curr[j] = bq_open.max(bq_ext);

            // Bt: gap in target (from left)
            let bt_open = m_curr[j - 1] + stack(q(i), Base::Gap, t(j - 1), t(j));
            let bt_ext = bt_curr[j - 1] + stack(Base::Gap, Base::Gap, t(j - 1), t(j));
            bt_curr[j] = bt_open.max(bt_ext);

            // Update best with terminal penalty (direction-dependent)
            let terminal = terminal_fn(q(i), t(j));
            let score_with_term = m_curr[j] + terminal;
            if score_with_term > best_score {
                best_score = score_with_term;
                best_i = i;
                best_j = j;
            }
        }
    }

    // Traceback
    let mut trace = Vec::new();
    let (mut i, mut j) = (best_i, best_j);
    while i > 0 && j > 0 {
        let op = tb[i][j];
        trace.push(op);
        match op {
            DpOp::Stop => break,
            DpOp::Match => {
                i -= 1;
                j -= 1;
            }
            DpOp::GapQ => {
                i -= 1;
            }
            DpOp::GapT => {
                j -= 1;
            }
        }
    }
    trace.reverse();

    DpExtension {
        score: best_score,
        q_len: best_i,
        t_len: best_j,
        trace,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Strand;

    #[test]
    fn test_extend_empty() {
        let result = extend(|_| Base::A, |_| Base::U, 0, 0, ExtendDir::Right);
        // Empty extension returns terminal penalty: DSM[Gap][A][Gap][U] (Gap→Q, Gap→T order)
        let expected = StackPair::new(Base::Gap, Base::A, Base::Gap, Base::U).energy() as i32;
        assert_eq!(result.score, expected);
        assert_eq!(result.q_len, 0);
    }

    #[test]
    fn test_single_match() {
        // Simplest case: extend by exactly 1 position in each direction
        // Q: A-C  (positions 0,1)
        // T: U-G  (positions 0,1, antiparallel)
        //
        // At (1,1): stack(q0=A, q1=C, t0=U, t1=G)
        // This should give us a stacking energy value
        let q = [Base::A, Base::C];
        let t = [Base::U, Base::G];

        let result = extend(|i| q[i.min(1)], |j| t[j.min(1)], 1, 1, ExtendDir::Right);

        // Calculate expected:
        // M[1,1] = M[0,0] + stack(A,C,U,G) = 0 + stack
        // Terminal = stack(Gap, C, Gap, G) (Gap→Q, Gap→T order per implementation)
        // Score = M[1,1] + terminal
        let stack_val = StackPair::new(Base::A, Base::C, Base::U, Base::G).energy() as i32;
        let terminal = StackPair::new(Base::Gap, Base::C, Base::Gap, Base::G).energy() as i32;
        let initial_term = StackPair::new(Base::Gap, Base::A, Base::Gap, Base::U).energy() as i32;

        println!("Single match test:");
        println!("  stack(A,C,U,G) = {}", stack_val);
        println!("  terminal(Gap,C,Gap,G) = {}", terminal);
        println!("  initial_terminal(Gap,A,Gap,U) = {}", initial_term);
        println!(
            "  Expected M[1,1] + term = {} + {} = {}",
            stack_val,
            terminal,
            stack_val + terminal
        );
        println!(
            "  Actual result: score={}, q_len={}, t_len={}",
            result.score, result.q_len, result.t_len
        );

        // The result should either be:
        // - initial_terminal (if no extension is better)
        // - stack_val + terminal (if extending is better)
        let expected = (stack_val + terminal).max(initial_term);
        assert_eq!(
            result.score, expected,
            "Score should match hand calculation"
        );
    }

    #[test]
    fn test_no_extension_better() {
        // Case where NOT extending gives better score than extending
        // Use bases that give unfavorable stacking
        let q = [Base::A, Base::A]; // AA
        let t = [Base::A, Base::A]; // AA (not complementary, should be unfavorable)

        let result = extend(|i| q[i.min(1)], |j| t[j.min(1)], 1, 1, ExtendDir::Right);

        let stack_val = StackPair::new(Base::A, Base::A, Base::A, Base::A).energy() as i32;
        let terminal = StackPair::new(Base::A, Base::Gap, Base::A, Base::Gap).energy() as i32;
        let initial_term = StackPair::new(Base::Gap, Base::A, Base::Gap, Base::A).energy() as i32;

        println!("No extension test:");
        println!("  stack(A,A,A,A) = {}", stack_val);
        println!("  terminal(A,Gap,A,Gap) = {}", terminal);
        println!("  initial_terminal(Gap,A,Gap,A) = {}", initial_term);
        println!(
            "  Extend score = {} + {} = {}",
            stack_val,
            terminal,
            stack_val + terminal
        );
        println!("  No extend score = {}", initial_term);
        println!("  Actual: score={}", result.score);
    }

    #[test]
    fn test_extend_right_with_seq() {
        let query = Seq::new(b"ACGUACGU", Strand::Forward);
        let target = Seq::new(b"UGCAUGCA", Strand::Forward);
        let result = extend_right(&query, &target, 2, 5, 10);
        assert!(result.score != i32::MIN / 2, "Should compute valid score");
    }

    #[test]
    fn test_extend_left_with_seq() {
        let query = Seq::new(b"ACGUACGU", Strand::Forward);
        let target = Seq::new(b"UGCAUGCA", Strand::Forward);
        let result = extend_left(&query, &target, 5, 2, 10);
        assert!(result.score != i32::MIN / 2, "Should compute valid score");
    }
}
