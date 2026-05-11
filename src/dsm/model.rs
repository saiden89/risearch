use crate::dp::MAX_ENERGY;
use crate::types::{Base, Energy, BASE_COUNT};

use super::{flat_idx, DsmTable, DSM_FLAT_SIZE};

const BASES: [Base; BASE_COUNT] = [Base::Gap, Base::A, Base::C, Base::G, Base::N, Base::U];

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
