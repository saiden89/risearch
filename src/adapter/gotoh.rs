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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Base, Energy};

    // Decimal positions make every tensor coordinate distinguishable. The
    // expected index order never goes through flat_idx or ScoringModel::score.
    fn source(a: usize, b: usize, c: usize, d: usize) -> i32 {
        (1000 * a + 100 * b + 10 * c + d) as i32
    }

    fn expected(a: usize, b: usize, c: usize, d: usize) -> i32 {
        let consumed = match (b, d) {
            (0, 0) => 0,
            (0, _) | (_, 0) => 1,
            _ if a == 0
                && c == 0
                && !Base::from_u8(b as u8)
                    .pair_type(Base::from_u8(d as u8))
                    .is_match(true) =>
            {
                0
            }
            _ => 2,
        };
        source(a, b, c, d) - 17 * consumed
    }

    #[test]
    fn scoring_tensor_and_every_adapter_lane_have_independent_expectations() {
        let mut table = [[[[0; 6]; 6]; 6]; 6];
        for (a, slab) in table.iter_mut().enumerate() {
            for (b, plane) in slab.iter_mut().enumerate() {
                for (c, row) in plane.iter_mut().enumerate() {
                    for (d, value) in row.iter_mut().enumerate() {
                        *value = source(a, b, c, d);
                    }
                }
            }
        }
        let model = ScoringModel::new(&table, Energy(12345), Energy(17));
        let left = model.transpose();
        for a in 0..6 {
            for b in 0..6 {
                let row = model.row_profile(a as u8, b as u8);
                assert_eq!(
                    model.extend_query_gap(a as u8, b as u8),
                    expected(a, b, 0, 0)
                );
                assert_eq!(row.ext_qgap(), expected(a, b, 0, 0));
                for d in 0..6 {
                    assert_eq!(
                        model.close_query_gap(a as u8, b as u8, d as u8),
                        expected(a, b, 0, d)
                    );
                    assert_eq!(
                        model.open_query_gap(a as u8, b as u8, d as u8),
                        expected(a, b, d, 0)
                    );
                    assert_eq!(model.boundary(b as u8, d as u8), expected(b, 0, d, 0));
                    // SAFETY: row borrows model and all offsets are within six-symbol lanes.
                    unsafe {
                        assert_eq!(*row.close_query_gap_ptr().add(d), expected(a, b, 0, d));
                        assert_eq!(*row.open_query_gap_ptr().add(d * 6), expected(a, b, d, 0));
                        assert_eq!(*row.boundary_ptr().add(d * 6), expected(b, 0, d, 0));
                    }
                }
                for c in 0..6 {
                    for d in 0..6 {
                        let (a8, b8, c8, d8) = (a as u8, b as u8, c as u8, d as u8);
                        assert_eq!(model.score(a8, b8, c8, d8), expected(a, b, c, d));
                        assert_eq!(left.score(a8, b8, c8, d8), expected(b, a, d, c));
                        assert_eq!(model.r#match(a8, b8, c8, d8), expected(a, b, c, d));
                        assert_eq!(model.close_target_gap(b8, c8, d8), expected(0, b, c, d));
                        assert_eq!(model.open_target_gap(b8, c8, d8), expected(b, 0, c, d));
                        assert_eq!(model.extend_target_gap(c8, d8), expected(0, 0, c, d));
                        // SAFETY: the model outlives row; each lane is indexed only
                        // within its documented six-symbol tensor dimensions.
                        unsafe {
                            assert_eq!(*row.match_ptr().add(c * 6 + d), expected(a, b, c, d));
                            assert_eq!(
                                *row.close_target_gap_ptr().add(c * 6 + d),
                                expected(0, b, c, d)
                            );
                            assert_eq!(
                                *row.open_target_gap_ptr().add(c * 6 + d),
                                expected(b, 0, c, d)
                            );
                            assert_eq!(
                                *row.extend_target_gap_ptr().add(c * 6 + d),
                                expected(0, 0, c, d)
                            );
                        }
                    }
                }
            }
        }
        assert_eq!(model.binding_energy(Energy(13579), 12), Energy(-1438));
    }
}
