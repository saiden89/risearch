use log::trace;

use super::{DpOp, DpView, MIN_SCORE, ScoreGrid, TracebackPath};

#[cfg_attr(feature = "prof", inline(never))]
pub(super) fn traceback(
    view: &DpView<'_>,
    m: &ScoreGrid,
    bq: &ScoreGrid,
    bt: &ScoreGrid,
    best_i: usize,
    best_j: usize,
    trace: &mut TracebackPath,
) {
    let (mut i, mut j) = (best_i, best_j);
    let mut state = DpOp::Match;

    while i > 0 || j > 0 {
        match state {
            DpOp::Stop => break,
            DpOp::Match if i > 0 && j > 0 => {
                trace.push(DpOp::Match);
                let m_val = m.get(i, j);
                let m_diag = m.get(i - 1, j - 1);
                let bq_diag = bq.get(i - 1, j - 1);
                let bt_diag = bt.get(i - 1, j - 1);

                // Check which transition produced this M value
                let match_e = view.match_e(i, j);
                let m_from_bq = view.m_from_bq(i, j);
                let m_from_bt = view.m_from_bt(i, j);

                let next = if m_diag > MIN_SCORE && m_val == m_diag + match_e {
                    DpOp::Match
                } else if bq_diag > MIN_SCORE && m_val == bq_diag + m_from_bq {
                    DpOp::GapQ
                } else if bt_diag > MIN_SCORE && m_val == bt_diag + m_from_bt {
                    DpOp::GapT
                } else {
                    DpOp::Stop
                };

                trace!("{} TB Match({},{}): next={:?}", view.dir, i, j, next);
                i -= 1;
                j -= 1;
                state = next;
            }
            DpOp::GapQ if i > 0 => {
                trace.push(DpOp::GapQ);
                let bq_val = bq.get(i, j);
                let m_up = m.get(i - 1, j);
                let bq_up = bq.get(i - 1, j);

                let bq_open = view.bq_open(i, j);
                let bq_ext = view.bq_ext(i);

                let next = if m_up > MIN_SCORE && bq_val == m_up + bq_open {
                    DpOp::Match
                } else if bq_up > MIN_SCORE && bq_val == bq_up + bq_ext {
                    DpOp::GapQ
                } else {
                    DpOp::Stop
                };

                trace!("{} TB GapQ({},{}): next={:?}", view.dir, i, j, next);
                i -= 1;
                state = next;
            }
            DpOp::GapT if j > 0 => {
                trace.push(DpOp::GapT);
                let bt_val = bt.get(i, j);
                let m_left = m.get(i, j - 1);
                let bt_left = bt.get(i, j - 1);

                let bt_open = view.bt_open(i, j);
                let bt_ext = view.bt_ext(j);

                let next = if m_left > MIN_SCORE && bt_val == m_left + bt_open {
                    DpOp::Match
                } else if bt_left > MIN_SCORE && bt_val == bt_left + bt_ext {
                    DpOp::GapT
                } else {
                    DpOp::Stop
                };

                trace!("{} TB GapT({},{}): next={:?}", view.dir, i, j, next);
                j -= 1;
                state = next;
            }
            _ => break,
        }
    }
}
