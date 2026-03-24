use crate::alignment::PairClass;
use crate::config::{ExtendConfig, ScoreConfig};
use crate::dp::gotoh::Gotoh;
use crate::types::Base;
use smallvec::SmallVec;

mod core;
pub mod gotoh;
mod init;
mod traceback;
use std::cmp::max;

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

/// Extension direction — determines sequence indexing polarity.
///
/// After scoring tables are built for each direction, extension direction
/// only affects sequence coordinate arithmetic. Stacking order is resolved
/// by the direction-canonical `Gotoh` tables.
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
// DP VIEW - Direction-aware sequence access (no scoring)
// =============================================================================

/// View into sequences for DP extension.
///
/// Pure sequence access — knows nothing about scoring or energy.
/// Direction-dependent base indexing is resolved here;
/// scoring is done via `Gotoh` transition tables.
pub struct DpView<'a> {
    query: &'a [Base],
    target: &'a [Base],
    q_anchor: usize,
    t_anchor: usize,
    pub(super) dir: ExtendDir,
    pub q_len: usize,
    pub t_len: usize,
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

    /// Create a directional extension view anchored at a seed boundary.
    pub fn new(
        query: &'a [Base],
        target: &'a [Base],
        q_anchor: usize,
        t_anchor: usize,
        dir: ExtendDir,
        max_ext: usize,
    ) -> Self {
        Self {
            query,
            target,
            q_anchor,
            t_anchor,
            dir,
            q_len: match dir {
                ExtendDir::Left => (q_anchor + 1).min(max_ext),
                ExtendDir::Right => (query.len() - q_anchor).min(max_ext),
            },
            t_len: match dir {
                ExtendDir::Left => (target.len() - t_anchor).min(max_ext),
                ExtendDir::Right => (t_anchor + 1).min(max_ext),
            },
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
    /// Bases are passed straight from the transformed index.
    #[inline(always)]
    pub fn t(&self, j: usize) -> usize {
        match self.dir {
            ExtendDir::Left => Self::right_base(self.target, self.t_anchor, j).idx(),
            ExtendDir::Right => Self::left_base(self.target, self.t_anchor, j).idx(),
        }
    }
}

/// Negative infinity for the (max, +) semiring over DP scores.
///
/// Must satisfy two invariants (enforced by compile-time assert below):
/// 1. Invalid scores can never drift into valid range through accumulated adds
/// 2. No i32 underflow from accumulated negative energy
pub(super) const NEG_INF: i32 = -1_000_000_000;

/// Conservative upper bound on |energy| from a single scoring table lookup.
/// Source tables are i16 (max 32767); penalty adds modest overhead.
/// Real values are ~300-400 (0.01 kcal/mol units), but we bound generously.
/// Enforced at runtime in ScoringModel::new.
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
pub struct DpGrid {
    data: Vec<DpCell>,
    width: usize,
}

impl DpGrid {
    pub fn new(max_extension: usize) -> Self {
        let side = max_extension.min(MAX_EXT).saturating_add(1).max(1);
        Self {
            data: vec![DpCell::EMPTY; side * side],
            width: side,
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
        // SAFETY: Bounds verified by debug_assert above. In release, callers
        // (init, core, traceback) only access indices within the grid allocation
        // established by grid.resize(t_len+1, q_len+1).
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

/// Guard holding a reference to the DP grid after forward pass.
/// While this exists, the grid cannot be reused for another extension
/// (borrow checker enforces this via the lifetime on `grid`).
/// Stores the `Gotoh` reference used during the forward pass so that
/// `traceback()` is guaranteed to use the same scoring tables.
pub struct ExtendResult<'a> {
    grid: &'a DpGrid,
    gotoh: &'a Gotoh,
    pub score: i32,
    pub q_len: usize,
    pub t_len: usize,
}

impl<'a> ExtendResult<'a> {
    /// Run traceback to reconstruct alignment as Pairings.
    /// Only call when alignment output is needed (skip for Minimal format).
    pub fn traceback(&self, view: &DpView<'_>) -> SmallVec<[PairClass; 64]> {
        let mut out = SmallVec::new();
        traceback(
            view, self.gotoh, self.grid, self.q_len, self.t_len, &mut out,
        );
        out
    }
}
