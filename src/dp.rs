//! Clean DP extension implementation from first principles.
//!
//! Uses nearest-neighbor stacking energy (DSM tables) with affine gap penalties.
//! Based on rust-bio patterns but adapted for RNA duplex alignment.

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
        // q_len: positions available from q_start down to 0 (q_start + 1 values)
        // t_len: positions available from t_start up to target end (target.len() - t_start values)
        let q_len = (q_start + 1).min(max_ext);
        let t_len = (target.len() - t_start).min(max_ext);

        const DIR: ExtendDir = ExtendDir::Left;
        trace!(
            "{} q_start={} t_start={} q_len={} t_len={}",
            DIR, q_start, t_start, q_len, t_len
        );

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
                trace!("{} best@(1,1): {} -> {}", DIR, best_e, val);
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

        // =======================================================================
        // LIMITED ROWS/COLUMNS INITIALIZATION (C parity: lines 314-351)
        // =======================================================================
        // WHY THIS EXISTS:
        // The DP has three states: M (match), Bq (query bulge), Bt (target bulge).
        // Transitions between states have structural constraints:
        //
        // 1. At row i=1: Bq is invalid (can't have query bulge with only 1 query base)
        //    Therefore M[2,j] can only come from M or Bt, never from Bq
        //
        // 2. At col j=1: Bt is invalid (can't have target bulge with only 1 target base)
        //    Therefore M[i,2] can only come from M or Bq, never from Bt
        //
        // The main DP loop (i≥3, j≥3) uses the general recurrence that considers
        // all three states. But rows i=1,2 and cols j=1,2 need special handling
        // because some transitions are impossible.
        //
        // Without this initialization, cells like Bt[1,3], M[2,3], etc. remain
        // unset, blocking paths through the DP matrix and producing suboptimal
        // alignments.
        // =======================================================================

        // Limited rows: Initialize Bt[1,j], M[2,j], Bq[2,j], Bt[2,j] for j >= 3
        for j in 3..t_len {
            // Bt[1, j] - gap extension in row 1
            let bt_1j = match (m.get(1, j - 1), bt.get(1, j - 1)) {
                (Some(m_val), Some(bt_val)) => {
                    let from_m = m_val + e(GAP, q_idx(1), t_idx(j), t_idx(j - 1));
                    let from_bt = bt_val + e(GAP, GAP, t_idx(j), t_idx(j - 1));
                    if from_m >= from_bt {
                        tb_bt.set(1, j, DpOp::Match);
                        Some(from_m)
                    } else {
                        tb_bt.set(1, j, DpOp::GapT);
                        Some(from_bt)
                    }
                }
                (Some(m_val), None) => {
                    tb_bt.set(1, j, DpOp::Match);
                    Some(m_val + e(GAP, q_idx(1), t_idx(j), t_idx(j - 1)))
                }
                (None, Some(bt_val)) => {
                    tb_bt.set(1, j, DpOp::GapT);
                    Some(bt_val + e(GAP, GAP, t_idx(j), t_idx(j - 1)))
                }
                (None, None) => None,
            };
            bt.set(1, j, bt_1j);

            // M[2, j] - can only come from M or Bt (not Bq, per C code comment)
            let m_2j = match (m.get(1, j - 1), bt.get(1, j - 1)) {
                (Some(m_val), Some(bt_val)) => {
                    let from_m = m_val + e(q_idx(2), q_idx(1), t_idx(j), t_idx(j - 1));
                    let from_bt = bt_val + e(q_idx(2), GAP, t_idx(j), t_idx(j - 1));
                    if from_m >= from_bt {
                        tb_m.set(2, j, DpOp::Match);
                        Some(from_m)
                    } else {
                        tb_m.set(2, j, DpOp::GapT);
                        Some(from_bt)
                    }
                }
                (Some(m_val), None) => {
                    tb_m.set(2, j, DpOp::Match);
                    Some(m_val + e(q_idx(2), q_idx(1), t_idx(j), t_idx(j - 1)))
                }
                (None, Some(bt_val)) => {
                    tb_m.set(2, j, DpOp::GapT);
                    Some(bt_val + e(q_idx(2), GAP, t_idx(j), t_idx(j - 1)))
                }
                (None, None) => None,
            };
            m.set(2, j, m_2j);

            // Check best for M[2, j]
            if let Some(v) = m_2j {
                let val = v + e(GAP, q_idx(2), GAP, t_idx(j));
                if val > best_e {
                    best_e = val;
                    best_i = 2;
                    best_j = j;
                }
            }

            // Bq[2, j]
            if let Some(m1j) = m.get(1, j) {
                bq.set(2, j, Some(m1j + e(q_idx(2), q_idx(1), GAP, t_idx(j))));
                tb_bq.set(2, j, DpOp::Match);
            }

            // Bt[2, j]
            let bt_2j = match (m.get(2, j - 1), bt.get(2, j - 1)) {
                (Some(m_val), Some(bt_val)) => {
                    let from_m = m_val + e(GAP, q_idx(2), t_idx(j), t_idx(j - 1));
                    let from_bt = bt_val + e(GAP, GAP, t_idx(j), t_idx(j - 1));
                    if from_m >= from_bt {
                        tb_bt.set(2, j, DpOp::Match);
                        Some(from_m)
                    } else {
                        tb_bt.set(2, j, DpOp::GapT);
                        Some(from_bt)
                    }
                }
                (Some(m_val), None) => {
                    tb_bt.set(2, j, DpOp::Match);
                    Some(m_val + e(GAP, q_idx(2), t_idx(j), t_idx(j - 1)))
                }
                (None, Some(bt_val)) => {
                    tb_bt.set(2, j, DpOp::GapT);
                    Some(bt_val + e(GAP, GAP, t_idx(j), t_idx(j - 1)))
                }
                (None, None) => None,
            };
            bt.set(2, j, bt_2j);
        }

        // Initialize limited columns (C parity: lines 334-351)
        // Col 1 Bq and col 2 M/Bt/Bq for i >= 3
        for i in 3..q_len {
            // Bq[i, 1] - gap extension in col 1
            let bq_i1 = match (m.get(i - 1, 1), bq.get(i - 1, 1)) {
                (Some(m_val), Some(bq_val)) => {
                    let from_m = m_val + e(q_idx(i), q_idx(i - 1), GAP, t_idx(1));
                    let from_bq = bq_val + e(q_idx(i), q_idx(i - 1), GAP, GAP);
                    if from_m >= from_bq {
                        tb_bq.set(i, 1, DpOp::Match);
                        Some(from_m)
                    } else {
                        tb_bq.set(i, 1, DpOp::GapQ);
                        Some(from_bq)
                    }
                }
                (Some(m_val), None) => {
                    tb_bq.set(i, 1, DpOp::Match);
                    Some(m_val + e(q_idx(i), q_idx(i - 1), GAP, t_idx(1)))
                }
                (None, Some(bq_val)) => {
                    tb_bq.set(i, 1, DpOp::GapQ);
                    Some(bq_val + e(q_idx(i), q_idx(i - 1), GAP, GAP))
                }
                (None, None) => None,
            };
            bq.set(i, 1, bq_i1);

            // M[i, 2] - can only come from M or Bq (not Bt, per C code comment)
            let m_i2 = match (m.get(i - 1, 1), bq.get(i - 1, 1)) {
                (Some(m_val), Some(bq_val)) => {
                    let from_m = m_val + e(q_idx(i), q_idx(i - 1), t_idx(2), t_idx(1));
                    let from_bq = bq_val + e(q_idx(i), q_idx(i - 1), t_idx(2), GAP);
                    if from_m >= from_bq {
                        tb_m.set(i, 2, DpOp::Match);
                        Some(from_m)
                    } else {
                        tb_m.set(i, 2, DpOp::GapQ);
                        Some(from_bq)
                    }
                }
                (Some(m_val), None) => {
                    tb_m.set(i, 2, DpOp::Match);
                    Some(m_val + e(q_idx(i), q_idx(i - 1), t_idx(2), t_idx(1)))
                }
                (None, Some(bq_val)) => {
                    tb_m.set(i, 2, DpOp::GapQ);
                    Some(bq_val + e(q_idx(i), q_idx(i - 1), t_idx(2), GAP))
                }
                (None, None) => None,
            };
            m.set(i, 2, m_i2);

            // Check best for M[i, 2]
            if let Some(v) = m_i2 {
                let val = v + e(GAP, q_idx(i), GAP, t_idx(2));
                if val > best_e {
                    best_e = val;
                    best_i = i;
                    best_j = 2;
                }
            }

            // Bt[i, 2]
            if let Some(mi1) = m.get(i, 1) {
                bt.set(i, 2, Some(mi1 + e(GAP, q_idx(i), t_idx(2), t_idx(1))));
                tb_bt.set(i, 2, DpOp::Match);
            }

            // Bq[i, 2]
            let bq_i2 = match (m.get(i - 1, 2), bq.get(i - 1, 2)) {
                (Some(m_val), Some(bq_val)) => {
                    let from_m = m_val + e(q_idx(i), q_idx(i - 1), GAP, t_idx(2));
                    let from_bq = bq_val + e(q_idx(i), q_idx(i - 1), GAP, GAP);
                    if from_m >= from_bq {
                        tb_bq.set(i, 2, DpOp::Match);
                        Some(from_m)
                    } else {
                        tb_bq.set(i, 2, DpOp::GapQ);
                        Some(from_bq)
                    }
                }
                (Some(m_val), None) => {
                    tb_bq.set(i, 2, DpOp::Match);
                    Some(m_val + e(q_idx(i), q_idx(i - 1), GAP, t_idx(2)))
                }
                (None, Some(bq_val)) => {
                    tb_bq.set(i, 2, DpOp::GapQ);
                    Some(bq_val + e(q_idx(i), q_idx(i - 1), GAP, GAP))
                }
                (None, None) => None,
            };
            bq.set(i, 2, bq_i2);
        }

        // Main DP loop (starts at i=3, j=3 since rows/cols 0-2 are initialized above)
        for i in 3..q_len {
            for j in 3..t_len {
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

                // Tie-breaking: Match > GapT > GapQ (C's max3 priority)
                let (val_m, step_m) = [(s_mq, DpOp::GapQ), (s_mt, DpOp::GapT), (s_mm, DpOp::Match)]
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
                        trace!(
                            "{} best@({},{}): {} -> {} step={:?}",
                            DIR, i, j, best_e, curr_e, step_m
                        );
                        best_e = curr_e;
                        best_i = i;
                        best_j = j;
                    }
                }

                // Bq[i,j] - query bulge state (gap in target)
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

                // Bt[i,j] - target bulge state (gap in query)
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

        // Traceback
        trace!(
            "{} TB start: best=({},{}) score={}",
            DIR, best_i, best_j, best_e
        );
        let mut trace_vec = Vec::new();
        let (mut i, mut j) = (best_i, best_j);
        let mut state = DpOp::Match;

        while i > 0 || j > 0 {
            match state {
                DpOp::Stop => break,
                DpOp::Match => {
                    if i == 0 || j == 0 {
                        break;
                    }
                    // Current operation is Match (diagonal move)
                    trace_vec.push(DpOp::Match);
                    let next_state = tb_m.get(i, j);
                    trace!("{} TB M({},{}): next={:?}", DIR, i, j, next_state);
                    i -= 1;
                    j -= 1;
                    state = match next_state {
                        DpOp::Stop => break,
                        DpOp::Match => DpOp::Match,
                        DpOp::GapQ => DpOp::GapQ,
                        DpOp::GapT => DpOp::GapT,
                    };
                }
                DpOp::GapQ => {
                    // Current operation is GapQ (query bulge, gap in target)
                    trace_vec.push(DpOp::GapQ);
                    let next_state = tb_bq.get(i, j);
                    trace!("{} TB Bq({},{}): next={:?}", DIR, i, j, next_state);
                    if i > 0 {
                        i -= 1;
                    } else {
                        break;
                    }
                    state = match next_state {
                        DpOp::Stop => break,
                        DpOp::Match => DpOp::Match,
                        DpOp::GapQ => DpOp::GapQ,
                        DpOp::GapT => DpOp::GapT,
                    };
                }
                DpOp::GapT => {
                    // Current operation is GapT (target bulge, gap in query)
                    trace_vec.push(DpOp::GapT);
                    let next_state = tb_bt.get(i, j);
                    trace!("{} TB Bt({},{}): next={:?}", DIR, i, j, next_state);
                    if j > 0 {
                        j -= 1;
                    } else {
                        break;
                    }
                    state = match next_state {
                        DpOp::Stop => break,
                        DpOp::Match => DpOp::Match,
                        DpOp::GapQ => DpOp::GapQ,
                        DpOp::GapT => DpOp::GapT,
                    };
                }
            }
        }

        trace!(
            "{} result: score={} q_len={} t_len={} trace={:?}",
            DIR, best_e, best_i, best_j, trace_vec
        );
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
        const DIR: ExtendDir = ExtendDir::Right;
        let mut best_e = e(q_idx(0), GAP, t_idx(0), GAP);
        let mut best_i = 0usize;
        let mut best_j = 0usize;
        trace!(
            "{} q_end={} t_end={} q_len={} t_len={} Q(0)={} T(0)={} init_e=DSM[{}][0][{}][0]={}",
            DIR,
            q_end,
            t_end,
            q_len,
            t_len,
            q_idx(0),
            t_idx(0),
            q_idx(0),
            t_idx(0),
            best_e
        );

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

        // =======================================================================
        // LIMITED ROWS/COLUMNS INITIALIZATION (C parity)
        // =======================================================================
        // Same rationale as extend_left - boundary rows/cols need special handling
        // because some DP state transitions are impossible at the edges.
        // Stacking order for right extension: [q(i-1)][q(i)][t(j-1)][t(j)]
        // =======================================================================

        // Limited rows: Initialize Bt[1,j], M[2,j], Bq[2,j], Bt[2,j] for j >= 3
        for j in 3..t_len {
            // Bt[1, j] - gap extension in row 1
            let bt_1j = match (m.get(1, j - 1), bt.get(1, j - 1)) {
                (Some(m_val), Some(bt_val)) => {
                    let from_m = m_val + e(q_idx(1), GAP, t_idx(j - 1), t_idx(j));
                    let from_bt = bt_val + e(GAP, GAP, t_idx(j - 1), t_idx(j));
                    if from_m >= from_bt {
                        tb_bt.set(1, j, DpOp::Match);
                        Some(from_m)
                    } else {
                        tb_bt.set(1, j, DpOp::GapT);
                        Some(from_bt)
                    }
                }
                (Some(m_val), None) => {
                    tb_bt.set(1, j, DpOp::Match);
                    Some(m_val + e(q_idx(1), GAP, t_idx(j - 1), t_idx(j)))
                }
                (None, Some(bt_val)) => {
                    tb_bt.set(1, j, DpOp::GapT);
                    Some(bt_val + e(GAP, GAP, t_idx(j - 1), t_idx(j)))
                }
                (None, None) => None,
            };
            bt.set(1, j, bt_1j);

            // M[2, j] - can only come from M or Bt (not Bq)
            let m_2j = match (m.get(1, j - 1), bt.get(1, j - 1)) {
                (Some(m_val), Some(bt_val)) => {
                    let from_m = m_val + e(q_idx(1), q_idx(2), t_idx(j - 1), t_idx(j));
                    let from_bt = bt_val + e(GAP, q_idx(2), t_idx(j - 1), t_idx(j));
                    if from_m >= from_bt {
                        tb_m.set(2, j, DpOp::Match);
                        Some(from_m)
                    } else {
                        tb_m.set(2, j, DpOp::GapT);
                        Some(from_bt)
                    }
                }
                (Some(m_val), None) => {
                    tb_m.set(2, j, DpOp::Match);
                    Some(m_val + e(q_idx(1), q_idx(2), t_idx(j - 1), t_idx(j)))
                }
                (None, Some(bt_val)) => {
                    tb_m.set(2, j, DpOp::GapT);
                    Some(bt_val + e(GAP, q_idx(2), t_idx(j - 1), t_idx(j)))
                }
                (None, None) => None,
            };
            m.set(2, j, m_2j);

            // Check best for M[2, j]
            if let Some(v) = m_2j {
                let val = v + e(q_idx(2), GAP, t_idx(j), GAP);
                if val > best_e {
                    best_e = val;
                    best_i = 2;
                    best_j = j;
                }
            }

            // Bq[2, j]
            if let Some(m1j) = m.get(1, j) {
                bq.set(2, j, Some(m1j + e(q_idx(1), q_idx(2), t_idx(j), GAP)));
                tb_bq.set(2, j, DpOp::Match);
            }

            // Bt[2, j]
            let bt_2j = match (m.get(2, j - 1), bt.get(2, j - 1)) {
                (Some(m_val), Some(bt_val)) => {
                    let from_m = m_val + e(q_idx(2), GAP, t_idx(j - 1), t_idx(j));
                    let from_bt = bt_val + e(GAP, GAP, t_idx(j - 1), t_idx(j));
                    if from_m >= from_bt {
                        tb_bt.set(2, j, DpOp::Match);
                        Some(from_m)
                    } else {
                        tb_bt.set(2, j, DpOp::GapT);
                        Some(from_bt)
                    }
                }
                (Some(m_val), None) => {
                    tb_bt.set(2, j, DpOp::Match);
                    Some(m_val + e(q_idx(2), GAP, t_idx(j - 1), t_idx(j)))
                }
                (None, Some(bt_val)) => {
                    tb_bt.set(2, j, DpOp::GapT);
                    Some(bt_val + e(GAP, GAP, t_idx(j - 1), t_idx(j)))
                }
                (None, None) => None,
            };
            bt.set(2, j, bt_2j);
        }

        // Limited columns: Initialize Bq[i,1], M[i,2], Bt[i,2], Bq[i,2] for i >= 3
        for i in 3..q_len {
            // Bq[i, 1] - gap extension in col 1
            let bq_i1 = match (m.get(i - 1, 1), bq.get(i - 1, 1)) {
                (Some(m_val), Some(bq_val)) => {
                    let from_m = m_val + e(q_idx(i - 1), q_idx(i), t_idx(1), GAP);
                    let from_bq = bq_val + e(q_idx(i - 1), q_idx(i), GAP, GAP);
                    if from_m >= from_bq {
                        tb_bq.set(i, 1, DpOp::Match);
                        Some(from_m)
                    } else {
                        tb_bq.set(i, 1, DpOp::GapQ);
                        Some(from_bq)
                    }
                }
                (Some(m_val), None) => {
                    tb_bq.set(i, 1, DpOp::Match);
                    Some(m_val + e(q_idx(i - 1), q_idx(i), t_idx(1), GAP))
                }
                (None, Some(bq_val)) => {
                    tb_bq.set(i, 1, DpOp::GapQ);
                    Some(bq_val + e(q_idx(i - 1), q_idx(i), GAP, GAP))
                }
                (None, None) => None,
            };
            bq.set(i, 1, bq_i1);

            // M[i, 2] - can only come from M or Bq (not Bt)
            let m_i2 = match (m.get(i - 1, 1), bq.get(i - 1, 1)) {
                (Some(m_val), Some(bq_val)) => {
                    let from_m = m_val + e(q_idx(i - 1), q_idx(i), t_idx(1), t_idx(2));
                    let from_bq = bq_val + e(q_idx(i - 1), q_idx(i), GAP, t_idx(2));
                    if from_m >= from_bq {
                        tb_m.set(i, 2, DpOp::Match);
                        Some(from_m)
                    } else {
                        tb_m.set(i, 2, DpOp::GapQ);
                        Some(from_bq)
                    }
                }
                (Some(m_val), None) => {
                    tb_m.set(i, 2, DpOp::Match);
                    Some(m_val + e(q_idx(i - 1), q_idx(i), t_idx(1), t_idx(2)))
                }
                (None, Some(bq_val)) => {
                    tb_m.set(i, 2, DpOp::GapQ);
                    Some(bq_val + e(q_idx(i - 1), q_idx(i), GAP, t_idx(2)))
                }
                (None, None) => None,
            };
            m.set(i, 2, m_i2);

            // Check best for M[i, 2]
            if let Some(v) = m_i2 {
                let val = v + e(q_idx(i), GAP, t_idx(2), GAP);
                if val > best_e {
                    best_e = val;
                    best_i = i;
                    best_j = 2;
                }
            }

            // Bt[i, 2]
            if let Some(mi1) = m.get(i, 1) {
                bt.set(i, 2, Some(mi1 + e(q_idx(i), GAP, t_idx(1), t_idx(2))));
                tb_bt.set(i, 2, DpOp::Match);
            }

            // Bq[i, 2]
            let bq_i2 = match (m.get(i - 1, 2), bq.get(i - 1, 2)) {
                (Some(m_val), Some(bq_val)) => {
                    let from_m = m_val + e(q_idx(i - 1), q_idx(i), t_idx(2), GAP);
                    let from_bq = bq_val + e(q_idx(i - 1), q_idx(i), GAP, GAP);
                    if from_m >= from_bq {
                        tb_bq.set(i, 2, DpOp::Match);
                        Some(from_m)
                    } else {
                        tb_bq.set(i, 2, DpOp::GapQ);
                        Some(from_bq)
                    }
                }
                (Some(m_val), None) => {
                    tb_bq.set(i, 2, DpOp::Match);
                    Some(m_val + e(q_idx(i - 1), q_idx(i), t_idx(2), GAP))
                }
                (None, Some(bq_val)) => {
                    tb_bq.set(i, 2, DpOp::GapQ);
                    Some(bq_val + e(q_idx(i - 1), q_idx(i), GAP, GAP))
                }
                (None, None) => None,
            };
            bq.set(i, 2, bq_i2);
        }

        // Main DP loop (starts at i=3, j=3 since rows/cols 0-2 are initialized above)
        for i in 3..q_len {
            for j in 3..t_len {

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

                // Debug trace for specific cells
                if (i == 5 && j == 7)
                    || (i == 4 && j == 6)
                    || (i == 3 && j == 5)
                    || (i == 2 && j == 4)
                {
                    let stack_e = e(q_idx(i - 1), q_idx(i), t_idx(j - 1), t_idx(j));
                    trace!(
                        "{} M[{},{}] stack=DSM[{}][{}][{}][{}]={} diag_m={:?} s_mm={:?}",
                        DIR,
                        i,
                        j,
                        q_idx(i - 1),
                        q_idx(i),
                        t_idx(j - 1),
                        t_idx(j),
                        stack_e,
                        m.diag(i, j),
                        s_mm
                    );
                }

                // Tie-breaking: Match > GapT > GapQ (C's max3 priority)
                let (val_m, step_m) = [(s_mq, DpOp::GapQ), (s_mt, DpOp::GapT), (s_mm, DpOp::Match)]
                    .into_iter()
                    .filter_map(|(opt, step)| opt.map(|v| (v, step)))
                    .max_by_key(|(v, _)| *v)
                    .unwrap_or((i32::MIN, DpOp::Stop));

                let val_m = if val_m == i32::MIN { None } else { Some(val_m) };
                m.set(i, j, val_m);
                tb_m.set(i, j, step_m);

                if let Some(v) = val_m {
                    let qi = q_idx(i);
                    let tj = t_idx(j);
                    let term_e = e(qi, GAP, tj, GAP);
                    let curr_e = v + term_e;
                    if curr_e > best_e {
                        trace!(
                            "{} best@({},{}): {} -> {} M={} term=DSM[{}][0][{}][0]={}",
                            DIR, i, j, best_e, curr_e, v, qi, tj, term_e
                        );
                        best_e = curr_e;
                        best_i = i;
                        best_j = j;
                    }
                }

                // Bq[i,j] - query bulge state (gap in target)
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

                // Bt[i,j] - target bulge state (gap in query)
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

        // Traceback
        let mut trace_vec = Vec::new();
        let (mut i, mut j) = (best_i, best_j);
        let mut state = DpOp::Match;

        while i > 0 || j > 0 {
            match state {
                DpOp::Stop => break,
                DpOp::Match => {
                    if i == 0 || j == 0 {
                        break;
                    }
                    // Current operation is Match (diagonal move)
                    trace_vec.push(DpOp::Match);
                    let next_state = tb_m.get(i, j);
                    i -= 1;
                    j -= 1;
                    state = match next_state {
                        DpOp::Stop => break,
                        DpOp::Match => DpOp::Match,
                        DpOp::GapQ => DpOp::GapQ,
                        DpOp::GapT => DpOp::GapT,
                    };
                }
                DpOp::GapQ => {
                    // Current operation is GapQ (query bulge, gap in target)
                    trace_vec.push(DpOp::GapQ);
                    let next_state = tb_bq.get(i, j);
                    if i > 0 {
                        i -= 1;
                    } else {
                        break;
                    }
                    state = match next_state {
                        DpOp::Stop => break,
                        DpOp::Match => DpOp::Match,
                        DpOp::GapQ => DpOp::GapQ,
                        DpOp::GapT => DpOp::GapT,
                    };
                }
                DpOp::GapT => {
                    // Current operation is GapT (target bulge, gap in query)
                    trace_vec.push(DpOp::GapT);
                    let next_state = tb_bt.get(i, j);
                    if j > 0 {
                        j -= 1;
                    } else {
                        break;
                    }
                    state = match next_state {
                        DpOp::Stop => break,
                        DpOp::Match => DpOp::Match,
                        DpOp::GapQ => DpOp::GapQ,
                        DpOp::GapT => DpOp::GapT,
                    };
                }
            }
        }

        trace!(
            "{} result: score={} q_len={} t_len={} trace={:?}",
            DIR, best_e, best_i, best_j, trace_vec
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
