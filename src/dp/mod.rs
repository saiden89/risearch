//! Gotoh 3-state DP over ranked symbol streams.
//!
//! A window arrives as dense `u8` ranks in DP order, every transition score
//! comes from a [`GotohScoring`] implementation, and results are plain `i32`.
//! The caller owns the symbols' meaning and any interpretation of a trace.

mod core;
pub mod gotoh;
mod init;
pub mod scoring;

pub use scoring::{GotohRowProfile, GotohScoring};
use std::cmp::max;

/// Maximum window length, per side, that [`gotoh::Gotoh::extend`] accepts.
///
/// Matches the typical max_ext parameter (100-200 symbols). Public because it is
/// a precondition of `extend`: the NEG_INF drift proof below is only valid
/// within this bound.
pub const MAX_EXT: usize = 256;

/// One step in a Gotoh traceback, named for the active DP state.
///
/// The variants specify only coordinate consumption:
/// - `Match`: one symbol from each input
/// - `GapQ`: one symbol from `q`
/// - `GapT`: one symbol from `t`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum TraceOp {
    Match,
    GapQ,
    GapT,
}

/// Negative infinity for the (max, +) semiring over DP scores.
///
/// Must satisfy two invariants (enforced by compile-time assert below):
/// 1. Invalid scores can never drift into valid range through accumulated adds
/// 2. No i32 underflow from accumulated negative scores
const NEG_INF: i32 = -1_500_000_000;

#[inline(always)]
pub(crate) fn is_valid_score(score: i32) -> bool {
    score > NEG_INF
}

/// Conservative upper bound on |score| from a single scoring transition.
/// This is a DP arithmetic bound; implementations of [`GotohScoring`] must keep
/// individual transition scores within it to prevent overflow during DP
/// accumulation.
pub(crate) const MAX_TRANSITION_SCORE: i64 = 1_200_000;

// Compile-time proof that NEG_INF arithmetic is safe for MAX_EXT.
const _: () = {
    // Longest path through MAX_EXT × MAX_EXT grid
    let max_drift = 2 * MAX_EXT as i64 * MAX_TRANSITION_SCORE;
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
    pub score: i32,
    /// Query-side offset of the best cell, in DP order from the anchor column.
    pub q_idx: usize,
    /// Target-side offset of the best cell, in DP order from the anchor column.
    pub t_idx: usize,
}

impl BestScore {
    pub(super) fn new(score: i32) -> Self {
        Self {
            score,
            q_idx: 0,
            t_idx: 0,
        }
    }

    /// Update if `val + term` exceeds current best.
    #[inline(always)]
    pub(super) fn update_if_better(&mut self, val: i32, term: i32, q_idx: usize, t_idx: usize) {
        let curr = val + term;
        if curr > self.score {
            self.score = curr;
            self.q_idx = q_idx;
            self.t_idx = t_idx;
        }
    }
}

/// Compares and returns the maximum of three values.
#[inline(always)]
pub(super) fn max3(a: i32, b: i32, c: i32) -> i32 {
    max(max(a, b), c)
}

/// Single DP cell: all three state scores packed together.
///
/// `repr(C)` guarantees field order and no padding (3 × i32 = 12 bytes).
/// Every access in the DP touches multiple states at the same (i,j),
/// so interleaving them maximizes cache line utilization.
#[repr(C)]
#[derive(Clone, Copy)]
pub(super) struct DpCell {
    pub(super) m: i32,     // Diagonal state: consume q and t
    pub(super) gap_q: i32, // Q-only state: consume q
    pub(super) gap_t: i32, // T-only state: consume t
}

impl DpCell {
    pub(super) const EMPTY: Self = Self {
        m: NEG_INF,
        gap_q: NEG_INF,
        gap_t: NEG_INF,
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
