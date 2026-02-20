use log::trace;
use smallvec::SmallVec;

use super::{DpGrid, DpView, MIN_SCORE};
use crate::alignment::PairClass;
use crate::dsm::{DsmModel, DSM_FLAT_SIZE};
use crate::types::Base;

/// Traceback state — internal to this module, never stored.
#[derive(Debug, Clone, Copy)]
enum State {
    Match,
    GapQ,
    GapT,
}

/// Check if a DP transition is valid: predecessor score is valid and produces the expected value.
#[inline(always)]
fn is_transition(val: i32, pred: i32, energy: i32) -> bool {
    pred > MIN_SCORE && val == pred + energy
}

/// Reconstruct alignment from score-only DP matrices, emitting `PairClass` path directly.
///
/// Walks backward from (best_i, best_j) comparing scores to determine transitions.
/// Uses the `DpView` to access sequence bases and classify each step in-place,
/// eliminating the need for a separate alignment reconstruction pass.
#[cfg_attr(feature = "prof", inline(never))]
pub(super) fn traceback<M: DsmModel>(
    view: &DpView<'_, M>,
    grid: &DpGrid,
    dsm_adjusted: &[i32; DSM_FLAT_SIZE],
    best_i: usize,
    best_j: usize,
    out: &mut SmallVec<[PairClass; 64]>,
) {
    let (mut i, mut j) = (best_i, best_j);
    let mut state = State::Match;

    while i > 0 || j > 0 {
        match state {
            State::Match if i > 0 && j > 0 => {
                let q_base = Base::from_idx(view.q(i));
                let t_base = Base::from_idx(view.t(j));
                out.push(PairClass::from_bases(q_base, t_base));

                let c = grid.get(i, j);
                let m_val = c.m;
                let diag = grid.get(i - 1, j - 1);

                let match_e = view.match_e(i, j, dsm_adjusted);
                let m_from_bq = view.m_from_bq(i, j, dsm_adjusted);
                let m_from_bt = view.m_from_bt(i, j, dsm_adjusted);

                let next = if is_transition(m_val, diag.m, match_e) {
                    Some(State::Match)
                } else if is_transition(m_val, diag.bq, m_from_bq) {
                    Some(State::GapQ)
                } else if is_transition(m_val, diag.bt, m_from_bt) {
                    Some(State::GapT)
                } else {
                    None
                };

                trace!(
                    "{} TB {:?}({},{}): next={:?}",
                    view.dir,
                    State::Match,
                    i,
                    j,
                    next
                );
                i -= 1;
                j -= 1;
                if let Some(next_state) = next {
                    state = next_state;
                } else {
                    break;
                }
            }
            State::GapQ if i > 0 => {
                out.push(PairClass::QueryBulge);

                let c = grid.get(i, j);
                let bq_val = c.bq;
                let up = grid.get(i - 1, j);

                let bq_open = view.bq_open(i, j, dsm_adjusted);
                let bq_ext = view.bq_ext(i, dsm_adjusted);

                let next = if is_transition(bq_val, up.m, bq_open) {
                    Some(State::Match)
                } else if is_transition(bq_val, up.bq, bq_ext) {
                    Some(State::GapQ)
                } else {
                    None
                };

                trace!(
                    "{} TB {:?}({},{}): next={:?}",
                    view.dir,
                    State::GapQ,
                    i,
                    j,
                    next
                );
                i -= 1;
                if let Some(next_state) = next {
                    state = next_state;
                } else {
                    break;
                }
            }
            State::GapT if j > 0 => {
                out.push(PairClass::TargetBulge);

                let c = grid.get(i, j);
                let bt_val = c.bt;
                let left = grid.get(i, j - 1);

                let bt_open = view.bt_open(i, j, dsm_adjusted);
                let bt_ext = view.bt_ext(j, dsm_adjusted);

                let next = if is_transition(bt_val, left.m, bt_open) {
                    Some(State::Match)
                } else if is_transition(bt_val, left.bt, bt_ext) {
                    Some(State::GapT)
                } else {
                    None
                };

                trace!(
                    "{} TB {:?}({},{}): next={:?}",
                    view.dir,
                    State::GapT,
                    i,
                    j,
                    next
                );
                j -= 1;
                if let Some(next_state) = next {
                    state = next_state;
                } else {
                    break;
                }
            }
            _ => break,
        }
    }
}
