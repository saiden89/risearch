use log::trace;
use smallvec::SmallVec;

use super::{DpView, ScoreGrid, MIN_SCORE};
use crate::alignment::Pairing;
use crate::dsm::DsmModel;
use crate::types::Base;

/// Traceback state — internal to this module, never stored.
#[derive(Debug, Clone, Copy)]
enum State {
    Match,
    GapQ,
    GapT,
    Stop,
}

impl std::fmt::Display for State {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            State::Match => write!(f, "Match"),
            State::GapQ => write!(f, "GapQ"),
            State::GapT => write!(f, "GapT"),
            State::Stop => write!(f, "Stop"),
        }
    }
}

/// Reconstruct alignment from score-only DP matrices, emitting `Pairing` directly.
///
/// Walks backward from (best_i, best_j) comparing scores to determine transitions.
/// Uses the `DpView` to access sequence bases and convert to `Pairing` in-place,
/// eliminating the need for a separate alignment reconstruction pass.
#[cfg_attr(feature = "prof", inline(never))]
pub(super) fn traceback<M: DsmModel>(
    view: &DpView<'_, M>,
    m: &ScoreGrid,
    bq: &ScoreGrid,
    bt: &ScoreGrid,
    best_i: usize,
    best_j: usize,
    out: &mut SmallVec<[Pairing; 64]>,
) {
    let (mut i, mut j) = (best_i, best_j);
    let mut state = State::Match;

    while i > 0 || j > 0 {
        match state {
            State::Stop => break,
            State::Match if i > 0 && j > 0 => {
                let q_base = Base::from_idx(view.q(i));
                let t_base = Base::from_idx(view.t(j));
                out.push(Pairing::from_bases(q_base, t_base));

                let m_val = m.get(i, j);
                let m_diag = m.get(i - 1, j - 1);
                let bq_diag = bq.get(i - 1, j - 1);
                let bt_diag = bt.get(i - 1, j - 1);

                let match_e = view.match_e(i, j);
                let m_from_bq = view.m_from_bq(i, j);
                let m_from_bt = view.m_from_bt(i, j);

                let next = if m_diag > MIN_SCORE && m_val == m_diag + match_e {
                    State::Match
                } else if bq_diag > MIN_SCORE && m_val == bq_diag + m_from_bq {
                    State::GapQ
                } else if bt_diag > MIN_SCORE && m_val == bt_diag + m_from_bt {
                    State::GapT
                } else {
                    State::Stop
                };

                trace!("{} TB {}({},{}): next={}", view.dir, State::Match, i, j, next);
                i -= 1;
                j -= 1;
                state = next;
            }
            State::GapQ if i > 0 => {
                let q_base = Base::from_idx(view.q(i));
                out.push(Pairing::QueryBulge(q_base));

                let bq_val = bq.get(i, j);
                let m_up = m.get(i - 1, j);
                let bq_up = bq.get(i - 1, j);

                let bq_open = view.bq_open(i, j);
                let bq_ext = view.bq_ext(i);

                let next = if m_up > MIN_SCORE && bq_val == m_up + bq_open {
                    State::Match
                } else if bq_up > MIN_SCORE && bq_val == bq_up + bq_ext {
                    State::GapQ
                } else {
                    State::Stop
                };

                trace!("{} TB {}({},{}): next={}", view.dir, State::GapQ, i, j, next);
                i -= 1;
                state = next;
            }
            State::GapT if j > 0 => {
                let t_base = Base::from_idx(view.t(j));
                out.push(Pairing::TargetBulge(t_base));

                let bt_val = bt.get(i, j);
                let m_left = m.get(i, j - 1);
                let bt_left = bt.get(i, j - 1);

                let bt_open = view.bt_open(i, j);
                let bt_ext = view.bt_ext(j);

                let next = if m_left > MIN_SCORE && bt_val == m_left + bt_open {
                    State::Match
                } else if bt_left > MIN_SCORE && bt_val == bt_left + bt_ext {
                    State::GapT
                } else {
                    State::Stop
                };

                trace!("{} TB {}({},{}): next={}", view.dir, State::GapT, i, j, next);
                j -= 1;
                state = next;
            }
            _ => break,
        }
    }
}
