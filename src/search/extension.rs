use smallvec::SmallVec;

use crate::adapter::gotoh::GotohModel;
use crate::alignment::PairClass;
use crate::dp::gotoh::Gotoh;
use crate::dp::{DpGrid, DpView, ExtendDir, TraceOp, MAX_EXT};
use crate::dsm::ScoringModel;
use crate::types::Energy;

use super::ExtensionResult;

/// Per-worker extension engine. Each worker owns one; `&mut` is safe because
/// only mutable state (grid, buffers) is mutated.
pub(super) struct ExtensionEngine {
    grid: DpGrid,
    q_buf: [u8; MAX_EXT],
    t_buf: [u8; MAX_EXT],
    gotoh_left: GotohModel,
    gotoh_right: GotohModel,
}

impl ExtensionEngine {
    pub(super) fn new(max_extension: usize, model: &ScoringModel) -> Self {
        let gotoh_right = Gotoh::new(model);
        let gotoh_left = Gotoh::new(&model.transpose());
        Self {
            grid: DpGrid::new(max_extension),
            q_buf: [0u8; MAX_EXT],
            t_buf: [0u8; MAX_EXT],
            gotoh_left,
            gotoh_right,
        }
    }

    pub(super) fn extend(&mut self, view: &DpView<'_>, include_alignment: bool) -> ExtensionResult {
        let gotoh = match view.dir() {
            ExtendDir::Left => &self.gotoh_left,
            ExtendDir::Right => &self.gotoh_right,
        };
        if view.is_empty() {
            return ExtensionResult {
                energy: Energy(
                    gotoh.boundary(view.q_anchor_base().as_u8(), view.t_anchor_base().as_u8()),
                ),
                q_ext: 0,
                t_ext: 0,
                pairs: None,
            };
        }

        let result = gotoh.extend(view, &mut self.grid, &mut self.q_buf, &mut self.t_buf);
        let pairs = if include_alignment && (result.q_idx > 0 || result.t_idx > 0) {
            let ops = gotoh.traceback(
                &self.q_buf,
                &self.t_buf,
                &self.grid,
                result.q_idx,
                result.t_idx,
            );
            Some(map_trace_to_pairs(view, &ops, result.q_idx, result.t_idx))
        } else {
            None
        };
        ExtensionResult {
            energy: Energy(result.energy),
            q_ext: result.q_idx,
            t_ext: result.t_idx,
            pairs,
        }
    }
}

/// Resolve a traceback into pair classes, ordered 5'->3' along the query.
///
/// A traceback runs from `(end_i, end_j)` back to the anchor, which is ascending
/// query coordinate for `Left` but descending for `Right` (`DpView` polarity), so
/// the right-hand walk is reversed to hand `Alignment` one canonical order.
fn map_trace_to_pairs(
    view: &DpView<'_>,
    ops: &[TraceOp],
    end_i: usize,
    end_j: usize,
) -> SmallVec<[PairClass; 64]> {
    let mut out = SmallVec::with_capacity(ops.len());
    let (mut i, mut j) = (end_i, end_j);
    for &op in ops {
        match op {
            TraceOp::Paired => {
                let target = view.t_base(j).complement();
                out.push(PairClass::from_bases(view.q_base(i), target));
                i -= 1;
                j -= 1;
            }
            TraceOp::GapQ => {
                out.push(PairClass::QueryBulge);
                i -= 1;
            }
            TraceOp::GapT => {
                out.push(PairClass::TargetBulge);
                j -= 1;
            }
        }
    }
    if view.dir() == ExtendDir::Right {
        out.reverse();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dsm::DsmRegistry;
    use crate::types::{Base, DsmId, Energy};

    #[test]
    fn extension_result_extents_match_traceback_consumption() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
        let mut engine = ExtensionEngine::new(8, &model);
        let query = [Base::A, Base::U, Base::G, Base::C];
        let target = [Base::G, Base::C, Base::A, Base::U];

        let view = DpView::new(&query, &target, 0, target.len() - 1, ExtendDir::Right, 8);
        let result = engine.extend(&view, true);
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
