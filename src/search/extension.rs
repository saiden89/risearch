use smallvec::SmallVec;

use crate::alignment::PairClass;
use crate::dp::gotoh::Gotoh;
use crate::dp::{DpGrid, DpView, ExtendDir, NEG_INF};
use crate::dsm::ScoringModel;
use crate::types::Base;

use super::ExtensionResult;

#[derive(Clone, Copy)]
enum State {
    Match,
    GapQ,
    GapT,
}

#[inline(always)]
fn is_transition(val: i32, pred: i32, energy: i32) -> bool {
    pred > NEG_INF && val == pred + energy
}

/// Per-worker extension engine. Each worker owns one; `&mut` is safe because
/// only `grid` is mutated — the rest is read-only extension state.
pub(super) struct ExtensionEngine {
    grid: DpGrid,
    max_extension: usize,
    gotoh_left: Gotoh,
    gotoh_right: Gotoh,
}

impl ExtensionEngine {
    pub(super) fn new(max_extension: usize, model: &ScoringModel) -> Self {
        let left_model = model.transpose();
        let gotoh_right = Gotoh::new(model);
        let gotoh_left = Gotoh::new(&left_model);
        Self {
            grid: DpGrid::new(max_extension),
            max_extension,
            gotoh_left,
            gotoh_right,
        }
    }

    pub(super) fn extend(
        &mut self,
        query: &[Base],
        target: &[Base],
        q_anchor: usize,
        t_anchor: usize,
        dir: ExtendDir,
        include_alignment: bool,
    ) -> ExtensionResult {
        let gotoh = match dir {
            ExtendDir::Left => &self.gotoh_left,
            ExtendDir::Right => &self.gotoh_right,
        };
        if self.max_extension == 0 {
            return ExtensionResult {
                score: gotoh.terminal_bases(query[q_anchor], target[t_anchor]),
                q_ext: 0,
                t_ext: 0,
                pairs: None,
            };
        }

        let view = DpView::new(query, target, q_anchor, t_anchor, dir, self.max_extension);
        let result = gotoh.extend(&view, &mut self.grid);
        ExtensionResult {
            score: result.score,
            q_ext: result.end_i,
            t_ext: result.end_j,
            pairs: (include_alignment && (result.end_i > 0 || result.end_j > 0))
                .then(|| self.traceback(gotoh, &view, result.end_i, result.end_j)),
        }
    }

    fn traceback(
        &self,
        gotoh: &Gotoh,
        view: &DpView<'_>,
        end_i: usize,
        end_j: usize,
    ) -> SmallVec<[PairClass; 64]> {
        let mut out = SmallVec::new();
        let (mut i, mut j) = (end_i, end_j);
        let mut state = State::Match;

        while i > 0 || j > 0 {
            let next = match state {
                State::Match if i > 0 && j > 0 => {
                    let q_base = view.q_base(i);
                    let t_base = view.t_base(j).complement();
                    out.push(PairClass::from_bases(q_base, t_base));

                    let c = self.grid.get(i, j);
                    let diag = self.grid.get(i - 1, j - 1);

                    let qi_prev = view.q_base(i - 1) as usize;
                    let qi = view.q_base(i) as usize;
                    let tj_prev = view.t_base(j - 1) as usize;
                    let tj = view.t_base(j) as usize;

                    let match_e = gotoh.match_energy(qi_prev, qi, tj_prev, tj);
                    let m_from_bq = gotoh.m_from_bq(qi_prev, qi, tj);
                    let m_from_bt = gotoh.m_from_bt(qi, tj_prev, tj);

                    i -= 1;
                    j -= 1;

                    if is_transition(c.m, diag.m, match_e) {
                        State::Match
                    } else if is_transition(c.m, diag.bq, m_from_bq) {
                        State::GapQ
                    } else if is_transition(c.m, diag.bt, m_from_bt) {
                        State::GapT
                    } else {
                        break;
                    }
                }
                State::GapQ if i > 0 => {
                    out.push(PairClass::QueryBulge);

                    let c = self.grid.get(i, j);
                    let up = self.grid.get(i - 1, j);

                    let qi_prev = view.q_base(i - 1) as usize;
                    let qi = view.q_base(i) as usize;
                    let tj = view.t_base(j) as usize;

                    let bq_open = gotoh.bq_open(qi_prev, qi, tj);
                    let bq_ext = gotoh.bq_extend(qi_prev, qi);

                    i -= 1;

                    if is_transition(c.bq, up.m, bq_open) {
                        State::Match
                    } else if is_transition(c.bq, up.bq, bq_ext) {
                        State::GapQ
                    } else {
                        break;
                    }
                }
                State::GapT if j > 0 => {
                    out.push(PairClass::TargetBulge);

                    let c = self.grid.get(i, j);
                    let left = self.grid.get(i, j - 1);

                    let qi = view.q_base(i) as usize;
                    let tj_prev = view.t_base(j - 1) as usize;
                    let tj = view.t_base(j) as usize;

                    let bt_open = gotoh.bt_open(qi, tj_prev, tj);
                    let bt_ext = gotoh.bt_extend_e(tj_prev, tj);

                    j -= 1;

                    if is_transition(c.bt, left.m, bt_open) {
                        State::Match
                    } else if is_transition(c.bt, left.bt, bt_ext) {
                        State::GapT
                    } else {
                        break;
                    }
                }
                _ => break,
            };
            state = next;
        }

        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Matrix;

    #[test]
    fn extension_result_extents_match_traceback_consumption() {
        let model = ScoringModel::new(Matrix::T04, 0);
        let mut engine = ExtensionEngine::new(8, &model);
        let query = [Base::A, Base::U, Base::G, Base::C];
        let target = [Base::G, Base::C, Base::A, Base::U];

        let result = engine.extend(&query, &target, 0, target.len() - 1, ExtendDir::Right, true);
        let pairs = result.pairs.expect("traceback expected");

        assert!(result.q_ext > 0);
        assert!(result.t_ext > 0);
        assert_eq!(
            pairs.iter().filter(|step| step.consumes_query()).count(),
            result.q_ext
        );
        assert_eq!(
            pairs.iter().filter(|step| step.consumes_target()).count(),
            result.t_ext
        );
    }
}
