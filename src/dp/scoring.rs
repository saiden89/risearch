//! DP scoring interface.
//!
//! This module defines the [`GotohScoring`] trait which decouples Gotoh
//! state-transition semantics from specific scoring implementations.
//!
//! Two access patterns are provided:
//!
//! - [`GotohRowProfile`]: row-local lookup lanes
//!   for the hot inner loop (fixed `(qp, qc)` per row, varying `(tp, tc)` per column).
//! - [`GotohScoring`]: Trait for point-lookup used in initialization and traceback.

/// Point-lookup interface for DP scoring.
///
/// Implementations must provide scores for all Gotoh 3-state transitions
/// given rank-indexed symbols.
pub trait GotohScoring {
    /// Associated row-profile type for hot-path optimization.
    type RowProfile: GotohRowProfile;

    /// Build a row-local lookup for the given (qp, qc) pair.
    fn row_profile(&self, qp: u8, qc: u8) -> Self::RowProfile;

    /// M ← M transition score: continue the match run.
    fn r#match(&self, qp: u8, qc: u8, tp: u8, tc: u8) -> i32;

    /// M ← GapQ transition score: close a query gap and return to matching.
    fn close_query_gap(&self, qp: u8, qc: u8, tc: u8) -> i32;

    /// M ← GapT transition score: close a target gap and return to matching.
    fn close_target_gap(&self, qc: u8, tp: u8, tc: u8) -> i32;

    /// GapQ ← M transition score: open a query gap.
    fn open_query_gap(&self, qp: u8, qc: u8, tc: u8) -> i32;

    /// GapQ ← GapQ transition score: extend a query gap.
    fn extend_query_gap(&self, qp: u8, qc: u8) -> i32;

    /// GapT ← M transition score: open a target gap.
    fn open_target_gap(&self, qc: u8, tp: u8, tc: u8) -> i32;

    /// GapT ← GapT transition score: extend a target gap.
    fn extend_target_gap(&self, tp: u8, tc: u8) -> i32;

    /// Boundary transition score.
    fn boundary(&self, qc: u8, tc: u8) -> i32;
}

/// Row-local precomputed lookup for the hot inner loop.
pub trait GotohRowProfile {
    /// Number of rank-indexed symbols in each row-local lookup lane.
    const SYMBOL_COUNT: usize;

    fn match_ptr(&self) -> *const i32;
    fn close_query_gap_ptr(&self) -> *const i32;
    fn open_query_gap_ptr(&self) -> *const i32;
    fn close_target_gap_ptr(&self) -> *const i32;
    fn open_target_gap_ptr(&self) -> *const i32;
    fn extend_target_gap_ptr(&self) -> *const i32;
    fn boundary_ptr(&self) -> *const i32;
    fn ext_qgap(&self) -> i32;
}
