use crate::dp::gotoh::Gotoh;
use crate::dp::scoring::{GotohRowProfile, GotohScoring};
use crate::dsm::{flat_idx, ScoringModel};
use crate::types::{BASE_COUNT, GAP};

/// Type alias for a Gotoh engine powered by the standard scoring model.
pub type GotohModel = Gotoh<ScoringModel>;

/// Row-local symbol lookup into the standard flat scoring table.
///
/// Constructed once per outer-loop iteration (fixed qp, qc pair). Inner-loop
/// lanes take only target-side symbol indices and return transition scores without
/// exposing the flat-table layout or the GAP=0 offset trick.
///
///  SAFETY: All symbol indices (`tp`, `tc`) must be in `0..SYMBOL_COUNT`.
///
/// The backing table pointer must remain valid for the lifetime of this struct.
/// Since `Gotoh` owns its scoring source, and `RowProfile` is
/// created and consumed within a single method call on `Gotoh`, the raw pointers are guaranteed
/// to be valid as long as the `Gotoh` instance and its owned scoring source
/// are alive.
#[derive(Clone, Copy)]
pub struct RowProfile {
    qp_qc_base: *const i32,
    close_qgap: *const i32,
    open_qgap: *const i32,
    close_tgap: *const i32,
    open_tgap: *const i32,
    extend_tgap: *const i32,
    boundary: *const i32,
    ext_qgap: i32,
}

impl RowProfile {
    /// Build a row-local lookup for the given (qp, qc) pair.
    ///
    /// `qp` and `qc` must be valid rank-indexed symbols for this scoring model.
    /// The scoring model must outlive the returned `RowProfile`.
    #[inline(always)]
    pub(crate) fn new(model: &ScoringModel, qp: u8, qc: u8) -> Self {
        // SAFETY: ScoringModel's table is a fixed-size array on the heap.
        // RowProfile is used strictly during the lifetime of the ScoringModel
        // reference passed to the Gotoh engine.
        let table = model.table_ptr();
        unsafe {
            RowProfile {
                qp_qc_base: table.add(flat_idx(qp, qc, 0, 0)),
                close_qgap: table.add(flat_idx(qp, qc, GAP, 0)),
                open_qgap: table.add(flat_idx(qp, qc, 0, GAP)),
                close_tgap: table.add(flat_idx(GAP, qc, 0, 0)),
                open_tgap: table.add(flat_idx(qc, GAP, 0, 0)),
                extend_tgap: table.add(flat_idx(GAP, GAP, 0, 0)),
                boundary: table.add(flat_idx(qc, GAP, 0, GAP)),
                ext_qgap: *table.add(flat_idx(qp, qc, GAP, GAP)),
            }
        }
    }
}

impl GotohRowProfile for RowProfile {
    /// DP sees rank-indexed symbols; this adapter stores one lane per concrete base.
    const SYMBOL_COUNT: usize = BASE_COUNT;

    #[inline(always)]
    fn match_ptr(&self) -> *const i32 {
        self.qp_qc_base
    }

    #[inline(always)]
    fn close_query_gap_ptr(&self) -> *const i32 {
        self.close_qgap
    }

    #[inline(always)]
    fn open_query_gap_ptr(&self) -> *const i32 {
        self.open_qgap
    }

    #[inline(always)]
    fn close_target_gap_ptr(&self) -> *const i32 {
        self.close_tgap
    }

    #[inline(always)]
    fn open_target_gap_ptr(&self) -> *const i32 {
        self.open_tgap
    }

    #[inline(always)]
    fn extend_target_gap_ptr(&self) -> *const i32 {
        self.extend_tgap
    }

    #[inline(always)]
    fn boundary_ptr(&self) -> *const i32 {
        self.boundary
    }

    #[inline(always)]
    fn ext_qgap(&self) -> i32 {
        self.ext_qgap
    }
}

impl GotohScoring for ScoringModel {
    type RowProfile = RowProfile;

    #[inline(always)]
    fn row_profile(&self, qp: u8, qc: u8) -> Self::RowProfile {
        RowProfile::new(self, qp, qc)
    }

    #[inline(always)]
    fn r#match(&self, qp: u8, qc: u8, tp: u8, tc: u8) -> i32 {
        self.score(qp, qc, tp, tc)
    }

    #[inline(always)]
    fn close_query_gap(&self, qp: u8, qc: u8, tc: u8) -> i32 {
        self.score(qp, qc, GAP, tc)
    }

    #[inline(always)]
    fn close_target_gap(&self, qc: u8, tp: u8, tc: u8) -> i32 {
        self.score(GAP, qc, tp, tc)
    }

    #[inline(always)]
    fn open_query_gap(&self, qp: u8, qc: u8, tc: u8) -> i32 {
        self.score(qp, qc, tc, GAP)
    }

    #[inline(always)]
    fn extend_query_gap(&self, qp: u8, qc: u8) -> i32 {
        self.score(qp, qc, GAP, GAP)
    }

    #[inline(always)]
    fn open_target_gap(&self, qc: u8, tp: u8, tc: u8) -> i32 {
        self.score(qc, GAP, tp, tc)
    }

    #[inline(always)]
    fn extend_target_gap(&self, tp: u8, tc: u8) -> i32 {
        self.score(GAP, GAP, tp, tc)
    }

    #[inline(always)]
    fn boundary(&self, qc: u8, tc: u8) -> i32 {
        self.score(qc, GAP, tc, GAP)
    }
}
