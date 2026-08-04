use std::ops::Range;

use log::trace;

use crate::adapter::gotoh::GotohModel;
use crate::alignment::AlignColumn;
#[cfg(test)]
use crate::alignment::PairClass;
use crate::dp::gotoh::Gotoh;
use crate::dp::{BestScore, DpGrid, TraceOp, MAX_EXT};
use crate::dsm::ScoringModel;
use crate::seed::SeedHit;
use crate::types::{Base, Energy};

/// A seed grown and scored in both directions: the final half-open ranges and
/// binding energy.
///
/// The target range is in the physical duplex frame the DP works in; converting
/// it to FASTA coordinates is the caller's job.
pub(super) struct SeedExtension {
    pub(super) q_range: Range<usize>,
    pub(super) t_range: Range<usize>,
    pub(super) energy: Energy,
}

/// Extension direction — the polarity of the DP window.
///
/// Query and target share one duplex-column frame: query bases run 5'→3' and
/// physical target bases run 3'→5'. Both coordinates therefore decrease during
/// left extension and increase during right extension, so a left flank enters
/// the DP buffers back-to-front.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ExtendDir {
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

/// Per-worker extension engine.
///
/// The engine owns the window policy along with the buffers it constrains: the
/// buffers are the reason a window can never exceed [`MAX_EXT`]. That ceiling is
/// mirrored upstream — `MAX_EXTENSION`, the clap range on `-l`, and
/// `check_unlimited_fits` — because those reject what the engine would otherwise
/// clamp silently. Raising it means changing all four.
pub(super) struct ExtensionEngine {
    grid: DpGrid,
    q_buf: [u8; MAX_EXT],
    t_buf: [u8; MAX_EXT],
    gotoh_left: GotohModel,
    gotoh_right: GotohModel,
    left_columns: Vec<AlignColumn>,
    right_columns: Vec<AlignColumn>,
    max_window: Option<usize>,
    build_alignment: bool,
}

impl ExtensionEngine {
    /// `max_window` of `None` is unlimited: each side follows the query.
    /// `build_alignment` controls traceback for every seed handled by the worker.
    pub(super) fn new(
        max_window: Option<usize>,
        build_alignment: bool,
        model: &ScoringModel,
    ) -> Self {
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
            // Reused across seeds. Keeping these out of the flank return value
            // leaves the no-alignment hot path with only `BestScore` results.
            left_columns: Vec::new(),
            right_columns: Vec::new(),
            max_window,
            build_alignment,
        }
    }

    /// Canonical-orientation model used for seed and final duplex scoring.
    fn model(&self) -> &ScoringModel {
        &self.gotoh_right.scoring
    }

    /// Per-side window cap for a flank with `query_avail` symbols available.
    ///
    /// An unlimited window follows the query, so the target window never outruns
    /// the query flank it pairs with. [`MAX_EXT`] bounds both, because that is
    /// what the buffers hold — resolving the sentinel here keeps the grid
    /// pre-size in the constructor and the per-extension cap on one rule. Even a
    /// zero extension keeps the anchor column needed for boundary scoring.
    fn window_cap(max_window: Option<usize>, query_avail: usize) -> usize {
        max_window.unwrap_or(query_avail).clamp(1, MAX_EXT)
    }

    /// Extend both flanks of `seed` against the full query and target views and
    /// score the whole duplex.
    ///
    /// Flank slicing, per-flank polarity, and prefix/suffix ordering all stay
    /// inside the engine: the anchor convention has one definition, the
    /// non-empty precondition of [`Self::extend_flank`] holds by construction, and
    /// callers see one extended duplex rather than two halves to reassemble.
    /// When alignment is enabled, the resolved flank columns remain in the
    /// engine until [`Self::materialize_alignment`] or the next extension.
    pub(super) fn extend_seed(
        &mut self,
        query: &[Base],
        target: &[Base],
        seed: &SeedHit,
    ) -> SeedExtension {
        let query_range = seed.query_range();
        let target_range = seed.target_range();
        let seed_query = &query[query_range.clone()];
        let seed_target = &target[target_range.clone()];
        let q_match_end = query_range
            .clone()
            .next_back()
            .expect("SeedHit query range must be non-empty");
        let t_match_end = target_range
            .clone()
            .next_back()
            .expect("SeedHit target range must be non-empty");

        let left = self.extend_flank(
            &query[..=query_range.start],
            &target[..=target_range.start],
            ExtendDir::Left,
        );
        let right = self.extend_flank(
            &query[q_match_end..],
            &target[t_match_end..],
            ExtendDir::Right,
        );

        let q_range = query_range.start - left.q_idx..query_range.end + right.q_idx;
        let t_range = target_range.start - left.t_idx..target_range.end + right.t_idx;
        let model = self.model();
        let score = model.ungapped_duplex_score(seed_query, seed_target)
            + Energy(left.score)
            + Energy(right.score);
        let energy = model.binding_energy(score, q_range.len() + t_range.len());

        SeedExtension {
            q_range,
            t_range,
            energy,
        }
    }

    /// Materialize the most recently extended seed's alignment.
    ///
    /// The flank buffers are overwritten by the next call to [`Self::extend_seed`],
    /// so a retained hit must call this before extension continues.
    pub(super) fn materialize_alignment(
        &self,
        query: &[Base],
        target: &[Base],
        seed: &SeedHit,
    ) -> Option<Box<[AlignColumn]>> {
        self.build_alignment.then(|| {
            assemble_alignment(
                &self.left_columns,
                &query[seed.query_range()],
                &target[seed.target_range()],
                &self.right_columns,
            )
        })
    }

    /// Extend one flank of a seed by DP.
    ///
    /// `query` and `target` are the flanking slices including the anchor column:
    /// they end at the seed boundary for `Left` and start at it for `Right`.
    fn extend_flank(&mut self, query: &[Base], target: &[Base], dir: ExtendDir) -> BestScore {
        assert!(
            !query.is_empty() && !target.is_empty(),
            "extension flanks must include the anchor column"
        );

        let gotoh = match dir {
            ExtendDir::Left => &self.gotoh_left,
            ExtendDir::Right => &self.gotoh_right,
        };
        let cap = Self::window_cap(self.max_window, query.len());
        let q_win = fill(&mut self.q_buf, query, dir, cap);
        let t_win = fill(&mut self.t_buf, target, dir, cap);
        trace!("{dir} q_len={} t_len={}", q_win.len(), t_win.len());

        let best = gotoh.extend(q_win, t_win, &mut self.grid);
        if self.build_alignment {
            let trace = gotoh.traceback(q_win, t_win, &self.grid, best.q_idx, best.t_idx);
            let columns = match dir {
                ExtendDir::Left => &mut self.left_columns,
                ExtendDir::Right => &mut self.right_columns,
            };
            map_trace_to_columns(q_win, t_win, dir, &trace, best.q_idx, best.t_idx, columns);
        }
        best
    }
}

/// Join the two reusable flank buffers around the seed, then compact the
/// finished alignment into its immutable representation.
fn assemble_alignment(
    left: &[AlignColumn],
    seed_query: &[Base],
    seed_target: &[Base],
    right: &[AlignColumn],
) -> Box<[AlignColumn]> {
    debug_assert_eq!(seed_query.len(), seed_target.len());

    let mut columns = Vec::with_capacity(left.len() + seed_query.len() + right.len());
    columns.extend_from_slice(left);
    columns.extend(
        seed_query
            .iter()
            .copied()
            .zip(seed_target.iter().copied())
            .map(|(query, target)| AlignColumn::paired(query, target)),
    );
    columns.extend_from_slice(right);
    columns.into_boxed_slice()
}

/// Copy a flank into `dst` in DP order — position 0 is the anchor column at the
/// seed boundary — and return the filled window.
///
/// `dst.len()` is the hard ceiling, so an unlimited `cap` still cannot overflow
/// the buffer.
///
/// A `Left` flank enters reversed here; `map_trace_to_columns` reverses a `Right`
/// traceback on the way out. Both stem from the same polarity, in opposite
/// directions — change one and check the other.
fn fill<'a>(dst: &'a mut [u8], src: &[Base], dir: ExtendDir, cap: usize) -> &'a [u8] {
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
    &dst[..len]
}

/// Resolve a traceback into columns ordered 5'->3' along the query.
///
/// A traceback runs from `(end_i, end_j)` back to the anchor, which is ascending
/// query coordinate for `Left` but descending for `Right` (window polarity), so
/// the right-hand walk is reversed into the physical duplex order used outside
/// the DP engine.
fn map_trace_to_columns(
    q: &[u8],
    t: &[u8],
    dir: ExtendDir,
    ops: &[TraceOp],
    end_i: usize,
    end_j: usize,
    out: &mut Vec<AlignColumn>,
) {
    out.clear();
    out.reserve(ops.len());
    let (mut i, mut j) = (end_i, end_j);
    for &op in ops {
        // Consumption follows the DP op, not the resolved class: a `Match` over
        // a gap-ranked symbol still consumed both inputs.
        let column = match op {
            TraceOp::Match => {
                let column = AlignColumn::paired(Base::from_u8(q[i]), Base::from_u8(t[j]));
                i -= 1;
                j -= 1;
                column
            }
            TraceOp::GapQ => {
                let column = AlignColumn::query_only(Base::from_u8(q[i]));
                i -= 1;
                column
            }
            TraceOp::GapT => {
                let column = AlignColumn::target_only(Base::from_u8(t[j]));
                j -= 1;
                column
            }
        };
        out.push(column);
    }
    if dir == ExtendDir::Right {
        out.reverse();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::DsmId;

    /// Same constructor the search path uses, so these pin the shipped model.
    fn test_model() -> ScoringModel {
        ScoringModel::load(&DsmId::from("t04"), 37, Energy::from_kcal(0.0)).unwrap()
    }

    #[test]
    fn flank_extents_match_traceback_consumption() {
        let model = test_model();
        let mut engine = ExtensionEngine::new(Some(8), true, &model);
        let query = [Base::A, Base::U, Base::G, Base::C];
        let target = [Base::U, Base::A, Base::C, Base::G];

        let result = engine.extend_flank(&query, &target, ExtendDir::Right);
        let columns = &engine.right_columns;
        assert!(!columns.is_empty(), "traceback expected");

        assert!(result.q_idx > 0);
        assert!(result.t_idx > 0);
        assert_eq!(
            columns
                .iter()
                .filter(|column| column.query() != Base::Gap)
                .count(),
            result.q_idx
        );
        assert_eq!(
            columns
                .iter()
                .filter(|column| column.target() != Base::Gap)
                .count(),
            result.t_idx
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
            assert_eq!(fill(&mut q_buf, q_src, dir, 3), expected_q.map(Base::as_u8));
            assert_eq!(fill(&mut t_buf, t_src, dir, 3), expected_t.map(Base::as_u8));
        }
    }

    #[test]
    fn zero_cap_window_scores_the_anchor_column() {
        let model = test_model();
        let mut engine = ExtensionEngine::new(Some(0), true, &model);
        let query = [Base::A, Base::U, Base::G, Base::C];
        let target = [Base::U, Base::A, Base::C, Base::G];

        // The anchor is the last base for Left and the first for Right, so the
        // two directions score different columns of the same slices.
        let left = engine.extend_flank(&query, &target, ExtendDir::Left);
        let right = engine.extend_flank(&query, &target, ExtendDir::Right);
        for result in [&left, &right] {
            assert_eq!((result.q_idx, result.t_idx), (0, 0));
        }
        assert!(engine.left_columns.is_empty());
        assert!(engine.right_columns.is_empty());
        assert_eq!(
            left.score,
            engine.gotoh_left.boundary(Base::C.as_u8(), Base::G.as_u8())
        );
        assert_eq!(
            right.score,
            engine
                .gotoh_right
                .boundary(Base::A.as_u8(), Base::U.as_u8())
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
        assert_eq!(fill(&mut buf, &query, ExtendDir::Right, cap).len(), 5);
        assert_eq!(fill(&mut buf, &long, ExtendDir::Right, cap).len(), 5);

        // Fixed: independent of how much sequence is available on either side.
        let cap = ExtensionEngine::window_cap(Some(20), query.len());
        assert_eq!(fill(&mut buf, &query, ExtendDir::Right, cap).len(), 5);
        assert_eq!(fill(&mut buf, &long, ExtendDir::Right, cap).len(), 20);

        // Buffer-safety invariant: the buffer is the ceiling on every path.
        for max_window in [None, Some(0), Some(20), Some(usize::MAX)] {
            let cap = ExtensionEngine::window_cap(max_window, long.len());
            assert!(fill(&mut buf, &long, ExtendDir::Right, cap).len() <= MAX_EXT);
        }
    }

    /// Same invariant as above, but through `extend_flank` — the cap has to reach
    /// the DP from the query flank, not merely be computable from it.
    #[test]
    fn unlimited_extension_matches_a_window_fixed_to_the_query_flank() {
        let model = test_model();
        let query = [Base::A, Base::U, Base::G, Base::C];
        // Target flank two orders of magnitude longer than the query, as in a real
        // chromosome search: the window must still follow the query.
        let target = [Base::U; 200];

        let mut unlimited = ExtensionEngine::new(None, true, &model);
        let mut fixed = ExtensionEngine::new(Some(query.len()), true, &model);
        for dir in [ExtendDir::Left, ExtendDir::Right] {
            let free = unlimited.extend_flank(&query, &target, dir);
            let pinned = fixed.extend_flank(&query, &target, dir);
            assert_eq!(
                (free.score, free.q_idx, free.t_idx),
                (pinned.score, pinned.q_idx, pinned.t_idx),
                "{dir}: unlimited must equal a window fixed to the query flank"
            );
            assert!(
                free.t_idx < query.len(),
                "{dir}: target extent {} escaped the {}-symbol query window",
                free.t_idx,
                query.len()
            );
        }
    }

    #[test]
    fn traceback_columns_stay_in_query_order_for_both_directions() {
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
            let mut columns = Vec::new();
            map_trace_to_columns(&q_buf, &t_buf, dir, &ops, 2, 2, &mut columns);
            assert_eq!(
                columns
                    .iter()
                    .map(|column| column.class())
                    .collect::<Vec<_>>(),
                expected
            );
        }
    }
}
