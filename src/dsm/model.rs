use crate::dp::MAX_ENERGY;
use crate::types::{Base, Energy, BASE_COUNT, GAP};

use super::{flat_idx, DsmTable, DSM_FLAT_SIZE};

const BASES: [Base; BASE_COUNT] = [Base::Gap, Base::A, Base::C, Base::G, Base::N, Base::U];

/// Row-local precomputed lookup into the flat scoring table.
///
/// Constructed once per outer-loop iteration (fixed qp, qc pair). Inner-loop
/// methods take only target-side indices and return transition scores without
/// exposing the flat-table layout or the GAP=0 offset trick.
///
/// # Safety
///
/// All index arguments (`tp`, `tc`) must be in `0..BASE_COUNT`.
/// The backing table pointer must remain valid for the lifetime of this struct.
#[derive(Clone, Copy)]
pub(crate) struct RowLookup {
    qp_qc_base: *const i32,
    close_qgap: *const i32,
    open_qgap: *const i32,
    close_tgap: *const i32,
    open_tgap: *const i32,
    extend_tgap: *const i32,
    terminal_base: *const i32,
    pub(crate) ext_qgap: i32,
}

impl RowLookup {
    #[inline(always)]
    pub(crate) fn stack(&self, tp: usize, tc: usize) -> i32 {
        unsafe { *self.qp_qc_base.add(tp * BASE_COUNT + tc) }
    }

    #[inline(always)]
    pub(crate) fn close_query_gap(&self, tc: usize) -> i32 {
        unsafe { *self.close_qgap.add(tc) }
    }

    #[inline(always)]
    pub(crate) fn open_query_gap(&self, tc: usize) -> i32 {
        unsafe { *self.open_qgap.add(tc * BASE_COUNT) }
    }

    #[inline(always)]
    pub(crate) fn close_target_gap(&self, tp: usize, tc: usize) -> i32 {
        unsafe { *self.close_tgap.add(tp * BASE_COUNT + tc) }
    }

    #[inline(always)]
    pub(crate) fn open_target_gap(&self, tp: usize, tc: usize) -> i32 {
        unsafe { *self.open_tgap.add(tp * BASE_COUNT + tc) }
    }

    #[inline(always)]
    pub(crate) fn extend_target_gap(&self, tp: usize, tc: usize) -> i32 {
        unsafe { *self.extend_tgap.add(tp * BASE_COUNT + tc) }
    }

    #[inline(always)]
    pub(crate) fn terminal(&self, tc: usize) -> i32 {
        unsafe { *self.terminal_base.add(tc * BASE_COUNT) }
    }
}

/// Penalty-adjusted flat scoring table for one orientation.
///
/// External inputs and outputs are kcal/mol. Internally DP uses integer scores:
/// `score = -delta_g * 10000`.
#[derive(Clone, Debug)]
pub struct ScoringModel {
    table: [i32; DSM_FLAT_SIZE],
    initiation: Energy,
    penalty: Energy,
}

impl ScoringModel {
    #[inline]
    fn ext_penalty_mult(q1: Base, q2: Base, t1_orig: Base, t2_orig: Base) -> i32 {
        if q2 == Base::Gap && t2_orig == Base::Gap {
            0
        } else if q2 == Base::Gap || t2_orig == Base::Gap {
            1
        } else if q1 == Base::Gap && t1_orig == Base::Gap {
            if q2.pair_type(t2_orig.complement()).is_match(true) {
                2
            } else {
                0
            }
        } else {
            2
        }
    }

    /// Create a new scoring model from a source table and initiation energy.
    pub fn new(source_table: &DsmTable, initiation: Energy, penalty: Energy) -> Self {
        let mut table = [0i32; DSM_FLAT_SIZE];
        for q1 in BASES {
            for q2 in BASES {
                for t1 in BASES {
                    let t1_orig = t1.complement();
                    for t2 in BASES {
                        let t2_orig = t2.complement();

                        let idx = flat_idx(q1.as_u8(), q2.as_u8(), t1.as_u8(), t2.as_u8());
                        let ext_penalty_mult = Self::ext_penalty_mult(q1, q2, t1_orig, t2_orig);

                        table[idx] = source_table[q1.as_usize()][q2.as_usize()][t1_orig.as_usize()]
                            [t2_orig.as_usize()]
                            - penalty.0 * ext_penalty_mult;
                    }
                }
            }
        }

        for (i, &val) in table.iter().enumerate() {
            assert!(
                (val as i64).unsigned_abs() <= MAX_ENERGY as u64,
                "table entry {} = {} exceeds MAX_ENERGY bound {}",
                i,
                val,
                MAX_ENERGY,
            );
        }

        Self {
            table,
            initiation,
            penalty,
        }
    }

    /// Compute the total binding free energy for a state with the given total
    /// stacking stability score and physical length.
    pub fn binding_energy(&self, stacking_stability: Energy, length: usize) -> Energy {
        self.initiation - stacking_stability - (self.penalty * length)
    }

    #[inline(always)]
    pub fn table_ptr(&self) -> *const i32 {
        self.table.as_ptr()
    }

    /// Build a row-local lookup for the given (qp, qc) pair.
    ///
    /// # Safety
    ///
    /// `qp` and `qc` must be in `0..BASE_COUNT` (valid base rank indices).
    #[inline(always)]
    pub(crate) fn row_lookup(&self, qp: u8, qc: u8) -> RowLookup {
        let table = self.table_ptr();
        unsafe {
            RowLookup {
                qp_qc_base: table.add(flat_idx(qp, qc, 0, 0)),
                close_qgap: table.add(flat_idx(qp, qc, GAP, 0)),
                open_qgap: table.add(flat_idx(qp, qc, 0, GAP)),
                close_tgap: table.add(flat_idx(GAP, qc, 0, 0)),
                open_tgap: table.add(flat_idx(qc, GAP, 0, 0)),
                extend_tgap: table.add(flat_idx(GAP, GAP, 0, 0)),
                terminal_base: table.add(flat_idx(qc, GAP, 0, GAP)),
                ext_qgap: *table.add(flat_idx(qp, qc, GAP, GAP)),
            }
        }
    }

    /// M ← M transition score: continue the paired stack.
    #[inline(always)]
    pub(crate) fn stack_score(&self, qp: u8, qc: u8, tp: u8, tc: u8) -> i32 {
        self.transition_score(qp, qc, tp, tc)
    }

    /// M ← Bq transition score: close a query gap and return to the stack.
    #[inline(always)]
    pub(crate) fn close_query_gap_score(&self, qp: u8, qc: u8, tc: u8) -> i32 {
        self.transition_score(qp, qc, GAP, tc)
    }

    /// M ← Bt transition score: close a target gap and return to the stack.
    #[inline(always)]
    pub(crate) fn close_target_gap_score(&self, qc: u8, tp: u8, tc: u8) -> i32 {
        self.transition_score(GAP, qc, tp, tc)
    }

    /// Bq ← M transition score: open a query gap.
    #[inline(always)]
    pub(crate) fn open_query_gap_score(&self, qp: u8, qc: u8, tc: u8) -> i32 {
        self.transition_score(qp, qc, tc, GAP)
    }

    /// Bq ← Bq transition score: extend a query gap.
    #[inline(always)]
    pub(crate) fn extend_query_gap_score(&self, qp: u8, qc: u8) -> i32 {
        self.transition_score(qp, qc, GAP, GAP)
    }

    /// Bt ← M transition score: open a target gap.
    #[inline(always)]
    pub(crate) fn open_target_gap_score(&self, qc: u8, tp: u8, tc: u8) -> i32 {
        self.transition_score(qc, GAP, tp, tc)
    }

    /// Bt ← Bt transition score: extend a target gap.
    #[inline(always)]
    pub(crate) fn extend_target_gap_score(&self, tp: u8, tc: u8) -> i32 {
        self.transition_score(GAP, GAP, tp, tc)
    }

    /// Terminal (boundary) penalty.
    #[inline(always)]
    pub(crate) fn terminal_penalty(&self, qc: u8, tc: u8) -> i32 {
        self.transition_score(qc, GAP, tc, GAP)
    }

    /// Produce a left-canonical (transposed) copy.
    pub fn transpose(&self) -> ScoringModel {
        let mut table = [0i32; DSM_FLAT_SIZE];
        for q1 in 0..BASE_COUNT as u8 {
            for q2 in 0..BASE_COUNT as u8 {
                for t1 in 0..BASE_COUNT as u8 {
                    for t2 in 0..BASE_COUNT as u8 {
                        let dst = flat_idx(q1, q2, t1, t2);
                        let src = flat_idx(q2, q1, t2, t1);
                        table[dst] = self.table[src];
                    }
                }
            }
        }
        Self {
            table,
            initiation: self.initiation,
            penalty: self.penalty,
        }
    }

    /// Point query for one transition score in the canonical 4-base tensor.
    #[inline(always)]
    pub fn transition_score(&self, q1: u8, q2: u8, t1: u8, t2: u8) -> i32 {
        debug_assert!(q1 < 6 && q2 < 6 && t1 < 6 && t2 < 6);
        let idx = flat_idx(q1, q2, t1, t2);
        // SAFETY: All args are in 0..6. Max idx = 1295 < 1296.
        unsafe { *self.table.get_unchecked(idx) }
    }

    /// Point query for one transition score using semantic bases.
    #[inline(always)]
    pub fn transition_score_bases(&self, q1: Base, q2: Base, t1: Base, t2: Base) -> i32 {
        self.transition_score(q1.as_u8(), q2.as_u8(), t1.as_u8(), t2.as_u8())
    }

    /// Calculate the thermodynamic score of a continuous, ungapped anti-parallel duplex.
    pub fn ungapped_duplex_score(
        &self,
        query: &[Base],
        target: &[Base],
        q_pos: usize,
        t_pos: usize,
        len: usize,
    ) -> Energy {
        if len <= 1 {
            return Energy(0);
        }
        let mut score = 0;
        let t_match_end = t_pos + len - 1;
        for i in 0..(len - 1) {
            score += self.transition_score_bases(
                query[q_pos + i],
                query[q_pos + i + 1],
                target[t_match_end - i],
                target[t_match_end - i - 1],
            );
        }
        Energy(score)
    }
}
