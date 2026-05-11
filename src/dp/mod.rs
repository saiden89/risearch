use crate::config::ExtendConfig;
use crate::types::Base;

mod core;
pub mod gotoh;
mod init;
use std::cmp::max;

/// Maximum extension length for precomputed DP lookup arrays.
/// Matches the typical max_ext parameter (100-200 bases).
pub(crate) const MAX_EXT: usize = 256;

/// DP runtime configuration derived from high-level search configs.
#[derive(Clone, Copy, Debug)]
pub(crate) struct DpConfig {
    max_extension: usize,
}

impl DpConfig {
    #[inline(always)]
    pub(crate) const fn max_extension(self) -> usize {
        self.max_extension
    }
}

impl From<&ExtendConfig> for DpConfig {
    fn from(extend: &ExtendConfig) -> Self {
        Self {
            max_extension: usize::from(extend.max_extension).min(MAX_EXT),
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

/// One step in a Gotoh traceback: which DP state was active.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TraceOp {
    Paired,
    GapQ,
    GapT,
}

// =============================================================================
// DP VIEW - Direction-aware sequence access (no scoring)
// =============================================================================

/// View into sequences for DP extension.
///
/// Pure sequence/coordinate view — knows nothing about scoring or energy.
/// It defines the extension window and how DP offsets map back to semantic
/// query/target bases. Scoring-specific lookup indices are materialized once
/// in `Gotoh::extend`.
pub struct DpView<'a> {
    query: &'a [Base],
    target: &'a [Base],
    q_anchor: usize,
    t_anchor: usize,
    dir: ExtendDir,
    q_len: usize,
    t_len: usize,
}

impl<'a> DpView<'a> {
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

    #[inline(always)]
    pub(crate) fn is_empty(&self) -> bool {
        self.q_len == 0 || self.t_len == 0
    }

    #[inline(always)]
    pub(crate) fn dir(&self) -> ExtendDir {
        self.dir
    }

    #[inline(always)]
    pub(crate) fn q_anchor_base(&self) -> Base {
        unsafe { *self.query.get_unchecked(self.q_anchor) }
    }

    #[inline(always)]
    pub(crate) fn t_anchor_base(&self) -> Base {
        unsafe { *self.target.get_unchecked(self.t_anchor) }
    }

    /// Get query base at DP position i (0 = anchor).
    #[inline(always)]
    pub(crate) fn q_base(&self, i: usize) -> Base {
        debug_assert!(i < self.q_len);
        let pos = match self.dir {
            ExtendDir::Left => self.q_anchor - i,
            ExtendDir::Right => self.q_anchor + i,
        };
        // SAFETY: q_len is derived from q_anchor/query.len() for the chosen dir,
        // so any i < q_len maps to an in-bounds query position.
        unsafe { *self.query.get_unchecked(pos) }
    }

    /// Get target base at DP position j (0 = anchor).
    #[inline(always)]
    pub(crate) fn t_base(&self, j: usize) -> Base {
        debug_assert!(j < self.t_len);
        let pos = match self.dir {
            ExtendDir::Left => self.t_anchor + j,
            ExtendDir::Right => self.t_anchor - j,
        };
        // SAFETY: t_len is derived from t_anchor/target.len() for the chosen dir,
        // so any j < t_len maps to an in-bounds target position.
        unsafe { *self.target.get_unchecked(pos) }
    }
}

/// Negative infinity for the (max, +) semiring over DP scores.
///
/// Must satisfy two invariants (enforced by compile-time assert below):
/// 1. Invalid scores can never drift into valid range through accumulated adds
/// 2. No i32 underflow from accumulated negative energy
const NEG_INF: i32 = -1_500_000_000;

#[inline(always)]
pub(crate) fn is_valid_score(score: i32) -> bool {
    score > NEG_INF
}

/// Conservative upper bound on |score| from a single scoring table lookup.
/// Tables use RIsearch3 score units (score = -kcal/mol * 10000). Worst case:
/// TSV energy 20 kcal/mol (200k score) + penalty 50 kcal/mol (500k score) × 2 = 1.2M.
/// Enforced at runtime in ScoringModel::from_source_table.
pub(crate) const MAX_ENERGY: i64 = 1_200_000;

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
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BestScore {
    pub energy: i32,
    pub q_idx: usize,
    pub t_idx: usize,
}

impl BestScore {
    pub(super) fn new(energy: i32) -> Self {
        Self {
            energy,
            q_idx: 0,
            t_idx: 0,
        }
    }

    /// Update if `val + term` exceeds current best.
    #[inline(always)]
    pub(super) fn update_if_better(&mut self, val: i32, term: i32, q_idx: usize, t_idx: usize) {
        let curr = val + term;
        if curr > self.energy {
            self.energy = curr;
            self.q_idx = q_idx;
            self.t_idx = t_idx;
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
    pub(super) fn as_mut_ptr(&mut self) -> *mut DpCell {
        self.data.as_mut_ptr()
    }

    #[inline(always)]
    pub(super) fn width(&self) -> usize {
        self.width
    }
}
