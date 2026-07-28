use log::trace;
use smallvec::SmallVec;

use crate::adapter::gotoh::GotohModel;
use crate::alignment::PairClass;
use crate::dp::gotoh::Gotoh;
use crate::dp::{DpGrid, TraceOp, MAX_EXT};
use crate::dsm::ScoringModel;
use crate::types::{Base, Energy};

use super::ExtensionResult;

/// Extension direction — the polarity of the DP window.
///
/// Query and target share one duplex-column frame: query bases run 5'→3' and
/// physical target bases run 3'→5'. Both coordinates therefore decrease during
/// left extension and increase during right extension, so a left flank enters
/// the DP buffers back-to-front.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum ExtendDir {
    Left,
    Right,
}

impl std::fmt::Display for ExtendDir {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Left => write!(f, "[EXT_LEFT]"),
            Self::Right => write!(f, "[EXT_RIGHT]"),
        }
    }
}

/// Per-worker extension engine. Each worker owns one; `&mut` is safe because
/// only mutable state (grid, buffers) is mutated.
///
/// The engine owns the window policy along with the buffers it constrains: the
/// buffers are the reason a window can never exceed [`MAX_EXT`], so nothing
/// upstream has to know that ceiling exists.
pub(super) struct ExtensionEngine {
    grid: DpGrid,
    q_buf: [u8; MAX_EXT],
    t_buf: [u8; MAX_EXT],
    gotoh_left: GotohModel,
    gotoh_right: GotohModel,
    max_window: Option<usize>,
}

impl ExtensionEngine {
    /// `max_window` of `None` is unlimited: each side follows the query.
    pub(super) fn new(max_window: Option<usize>, model: &ScoringModel) -> Self {
        let gotoh_right = Gotoh::new(model);
        let gotoh_left = Gotoh::new(&model.transpose());
        Self {
            // Pre-size to the cap of the longest query the buffers can hold, so
            // no extension pays for a regrow.
            grid: DpGrid::new(Self::window_cap(max_window, MAX_EXT)),
            q_buf: [0u8; MAX_EXT],
            t_buf: [0u8; MAX_EXT],
            gotoh_left,
            gotoh_right,
            max_window,
        }
    }

    /// Per-side window cap for a flank with `query_avail` symbols available.
    ///
    /// An unlimited window follows the query, so the target window never outruns
    /// the query flank it pairs with. [`MAX_EXT`] bounds both, because that is
    /// what the buffers hold — resolving the sentinel here keeps the grid
    /// pre-size in `new` and the per-extension cap on one rule.
    fn window_cap(max_window: Option<usize>, query_avail: usize) -> usize {
        max_window.unwrap_or(query_avail).min(MAX_EXT)
    }

    /// Extend one flank of a seed by DP.
    ///
    /// `query` and `target` are the flanking slices including the anchor column:
    /// they end at the seed boundary for `Left` and start at it for `Right`.
    pub(super) fn extend(
        &mut self,
        query: &[Base],
        target: &[Base],
        dir: ExtendDir,
        include_alignment: bool,
    ) -> ExtensionResult {
        assert!(
            !query.is_empty() && !target.is_empty(),
            "extension flanks must include the anchor column"
        );

        let gotoh = match dir {
            ExtendDir::Left => &self.gotoh_left,
            ExtendDir::Right => &self.gotoh_right,
        };
        let cap = Self::window_cap(self.max_window, query.len());
        let q_len = fill(&mut self.q_buf, query, dir, cap);
        let t_len = fill(&mut self.t_buf, target, dir, cap);
        trace!("{dir} q_len={q_len} t_len={t_len}");

        if q_len == 0 || t_len == 0 {
            return ExtensionResult {
                energy: Energy(gotoh.boundary(anchor(query, dir), anchor(target, dir))),
                q_ext: 0,
                t_ext: 0,
                pairs: None,
            };
        }

        let (q_win, t_win) = (&self.q_buf[..q_len], &self.t_buf[..t_len]);
        let best = gotoh.extend(q_win, t_win, &mut self.grid);
        let pairs = if include_alignment && (best.q_idx > 0 || best.t_idx > 0) {
            let trace = gotoh.traceback(q_win, t_win, &self.grid, best.q_idx, best.t_idx);
            Some(map_trace_to_pairs(
                q_win, t_win, dir, &trace, best.q_idx, best.t_idx,
            ))
        } else {
            None
        };
        ExtensionResult {
            energy: Energy(best.score),
            q_ext: best.q_idx,
            t_ext: best.t_idx,
            pairs,
        }
    }
}

/// Copy a flank into `dst` in DP order — position 0 is the anchor column at the
/// seed boundary — and return the window length.
///
/// `dst.len()` is the hard ceiling, so an unlimited `cap` still cannot overflow
/// the buffer.
fn fill(dst: &mut [u8], src: &[Base], dir: ExtendDir, cap: usize) -> usize {
    let len = src.len().min(cap).min(dst.len());
    match dir {
        ExtendDir::Left => {
            for (dst, base) in dst[..len].iter_mut().zip(src.iter().rev()) {
                *dst = base.as_u8();
            }
        }
        ExtendDir::Right => {
            for (dst, base) in dst[..len].iter_mut().zip(src) {
                *dst = base.as_u8();
            }
        }
    }
    len
}

/// Anchor symbol at the seed boundary, for windows too short to fill.
///
/// Routed through `fill` so the polarity rule has exactly one definition.
fn anchor(src: &[Base], dir: ExtendDir) -> u8 {
    let mut buf = [0u8; 1];
    assert_eq!(
        fill(&mut buf, src, dir, 1),
        1,
        "flank must include the anchor column"
    );
    buf[0]
}

/// Resolve a traceback into pair classes, ordered 5'->3' along the query.
///
/// A traceback runs from `(end_i, end_j)` back to the anchor, which is ascending
/// query coordinate for `Left` but descending for `Right` (window polarity), so
/// the right-hand walk is reversed to hand `Alignment` one canonical order.
fn map_trace_to_pairs(
    q: &[u8],
    t: &[u8],
    dir: ExtendDir,
    ops: &[TraceOp],
    end_i: usize,
    end_j: usize,
) -> SmallVec<[PairClass; 64]> {
    let mut out = SmallVec::with_capacity(ops.len());
    let (mut i, mut j) = (end_i, end_j);
    for &op in ops {
        match op {
            TraceOp::Match => {
                out.push(PairClass::from_bases(
                    Base::from_u8(q[i]),
                    Base::from_u8(t[j]),
                ));
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
    if dir == ExtendDir::Right {
        out.reverse();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dsm::DsmRegistry;
    use crate::types::DsmId;

    #[test]
    fn extension_result_extents_match_traceback_consumption() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
        let mut engine = ExtensionEngine::new(Some(8), &model);
        let query = [Base::A, Base::U, Base::G, Base::C];
        let target = [Base::U, Base::A, Base::C, Base::G];

        let result = engine.extend(&query, &target, ExtendDir::Right, true);
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

    #[test]
    fn fill_walks_query_and_physical_target_together() {
        let query = [Base::Gap, Base::A, Base::C, Base::G, Base::U];
        let target = [Base::N, Base::U, Base::G, Base::C, Base::A, Base::Gap];
        for (q_src, t_src, dir, expected_q, expected_t) in [
            (
                &query[..4],
                &target[..5],
                ExtendDir::Left,
                [Base::G, Base::C, Base::A],
                [Base::A, Base::C, Base::G],
            ),
            (
                &query[1..],
                &target[2..],
                ExtendDir::Right,
                [Base::A, Base::C, Base::G],
                [Base::G, Base::C, Base::A],
            ),
        ] {
            let mut q_buf = [u8::MAX; 3];
            let mut t_buf = [u8::MAX; 3];
            assert_eq!(fill(&mut q_buf, q_src, dir, 3), 3);
            assert_eq!(fill(&mut t_buf, t_src, dir, 3), 3);
            assert_eq!(q_buf, expected_q.map(Base::as_u8));
            assert_eq!(t_buf, expected_t.map(Base::as_u8));
        }
    }

    #[test]
    fn zero_cap_window_scores_the_anchor_column() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
        let mut engine = ExtensionEngine::new(Some(0), &model);
        let query = [Base::A, Base::U, Base::G, Base::C];
        let target = [Base::U, Base::A, Base::C, Base::G];

        // The anchor is the last base for Left and the first for Right, so the
        // two directions score different columns of the same slices.
        let left = engine.extend(&query, &target, ExtendDir::Left, true);
        let right = engine.extend(&query, &target, ExtendDir::Right, true);
        for result in [&left, &right] {
            assert_eq!((result.q_ext, result.t_ext), (0, 0));
            assert!(result.pairs.is_none());
        }
        assert_eq!(
            left.energy,
            Energy(engine.gotoh_left.boundary(Base::C.as_u8(), Base::G.as_u8()))
        );
        assert_eq!(
            right.energy,
            Energy(
                engine
                    .gotoh_right
                    .boundary(Base::A.as_u8(), Base::U.as_u8())
            )
        );
    }

    #[test]
    fn unlimited_window_follows_the_query_not_the_target() {
        let query = [Base::A; 5];
        let long = [Base::A; MAX_EXT + 10];
        let mut buf = [0u8; MAX_EXT];

        // Unlimited: both sides take the query flank, so a target flank orders of
        // magnitude longer than the query cannot outrun the query it pairs with.
        let cap = ExtensionEngine::window_cap(None, query.len());
        assert_eq!(fill(&mut buf, &query, ExtendDir::Right, cap), 5);
        assert_eq!(fill(&mut buf, &long, ExtendDir::Right, cap), 5);

        // Fixed: independent of how much sequence is available on either side.
        let cap = ExtensionEngine::window_cap(Some(20), query.len());
        assert_eq!(fill(&mut buf, &query, ExtendDir::Right, cap), 5);
        assert_eq!(fill(&mut buf, &long, ExtendDir::Right, cap), 20);

        // Buffer-safety invariant: the buffer is the ceiling on every path.
        for max_window in [None, Some(0), Some(20), Some(usize::MAX)] {
            let cap = ExtensionEngine::window_cap(max_window, long.len());
            assert!(fill(&mut buf, &long, ExtendDir::Right, cap) <= MAX_EXT);
        }
    }

    /// Same invariant as above, but through `extend` — the cap has to reach the
    /// DP from the query flank, not merely be computable from it.
    #[test]
    fn unlimited_extension_matches_a_window_fixed_to_the_query_flank() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
        let query = [Base::A, Base::U, Base::G, Base::C];
        // Target flank two orders of magnitude longer than the query, as in a real
        // chromosome search: the window must still follow the query.
        let target = [Base::U; 200];

        let mut unlimited = ExtensionEngine::new(None, &model);
        let mut fixed = ExtensionEngine::new(Some(query.len()), &model);
        for dir in [ExtendDir::Left, ExtendDir::Right] {
            let free = unlimited.extend(&query, &target, dir, true);
            let pinned = fixed.extend(&query, &target, dir, true);
            assert_eq!(
                (free.energy, free.q_ext, free.t_ext),
                (pinned.energy, pinned.q_ext, pinned.t_ext),
                "{dir}: unlimited must equal a window fixed to the query flank"
            );
            assert!(
                free.t_ext < query.len(),
                "{dir}: target extent {} escaped the {}-symbol query window",
                free.t_ext,
                query.len()
            );
        }
    }

    #[test]
    fn traceback_pairs_stay_in_query_order_for_both_directions() {
        let query = [Base::A, Base::G, Base::C];
        let target = [Base::U, Base::U, Base::G];
        let ops = [TraceOp::Match, TraceOp::Match];

        for (dir, expected) in [
            (ExtendDir::Right, [PairClass::Wobble, PairClass::Canonical]),
            (ExtendDir::Left, [PairClass::Canonical, PairClass::Wobble]),
        ] {
            let mut q_buf = [0u8; 3];
            let mut t_buf = [0u8; 3];
            fill(&mut q_buf, &query, dir, 3);
            fill(&mut t_buf, &target, dir, 3);
            assert_eq!(
                map_trace_to_pairs(&q_buf, &t_buf, dir, &ops, 2, 2).as_slice(),
                expected
            );
        }
    }
}
