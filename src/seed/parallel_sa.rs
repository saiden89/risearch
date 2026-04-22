//! Parallel suffix-array traversal kernel for seed search.


use std::ops::Range;

use crate::types::Base;

use super::SeedView;

/// A raw seed candidate emitted by the parallel SA traversal.
#[derive(Debug, Clone)]
pub(super) struct SeedMatch {
    pub query_interval: Range<usize>,
    pub target_interval: Range<usize>,
    pub seed_len: usize,
}

pub(super) fn traverse<const WOBBLE: bool, F: FnMut(SeedMatch)>(
    q: SeedView<'_>,
    t: SeedView<'_>,
    min_len: usize,
    max_len: usize,
    max_mm: usize,
    min_prefix: usize,
    min_suffix: usize,
    on_match: &mut F,
) {
    if q.sa_real_len == 0 || t.sa_real_len == 0 || min_len > max_len {
        return;
    }

    let mut traverser = ParallelSaTraverser {
        q,
        t,
        min_len,
        max_len,
        max_mm,
        min_prefix,
        min_suffix,
        on_match,
    };
    let q_range = 0..traverser.q.sa_real_len;
    let t_range = 0..traverser.t.sa_real_len;

    traverser.recurse::<WOBBLE>(q_range, t_range, 0, 0, 0);
}

struct ParallelSaTraverser<'a, F: FnMut(SeedMatch)> {
    q: SeedView<'a>,
    t: SeedView<'a>,
    min_len: usize,
    max_len: usize,
    max_mm: usize,
    min_prefix: usize,
    min_suffix: usize,
    on_match: &'a mut F,
}

impl<F: FnMut(SeedMatch)> ParallelSaTraverser<'_, F> {
    #[inline(always)]
    fn should_emit(&self, depth: usize, match_streak: usize, mm_count: usize) -> bool {
        depth >= self.min_len
            && depth <= self.max_len
            && (mm_count == 0 || (match_streak >= self.min_suffix && match_streak < self.min_len))
    }

    #[inline(always)]
    fn can_reach_suffix(&self, depth: usize, match_streak: usize, mm_count: usize) -> bool {
        mm_count == 0
            || self.min_suffix == 0
            || match_streak + (self.max_len - depth) >= self.min_suffix
    }

    #[inline(always)]
    fn can_mismatch_next(&self, depth: usize, match_streak: usize, mm_count: usize) -> bool {
        let next_depth = depth + 1;
        self.max_mm > 0
            && mm_count < self.max_mm
            && next_depth > self.min_prefix
            && match_streak < self.min_len
            && self.max_len - next_depth >= self.min_suffix
    }
    #[inline(always)]
    fn recurse<const WOBBLE: bool>(
        &mut self,
        query_sa: Range<usize>,
        target_sa: Range<usize>,
        depth: usize,
        match_streak: usize,
        mm_count: usize,
    ) {
        let (query_start, query_end) = (query_sa.start, query_sa.end);
        let (target_start, target_end) = (target_sa.start, target_sa.end);

        if self.should_emit(depth, match_streak, mm_count) {
            (self.on_match)(SeedMatch {
                query_interval: query_start..query_end,
                target_interval: target_start..target_end,
                seed_len: depth,
            });
        }

        if depth >= self.max_len {
            return;
        }

        if !self.can_reach_suffix(depth, match_streak, mm_count) {
            return;
        }

        if query_end - query_start == 1 && target_end - target_start == 1 {
            self.recurse_singleton::<WOBBLE>(
                query_start,
                target_start,
                depth,
                match_streak,
                mm_count,
            );
            return;
        }
        if query_end - query_start == 1 {
            self.recurse_half_singleton::<WOBBLE, true>(
                query_start,
                target_start..target_end,
                depth,
                match_streak,
                mm_count,
            );
            return;
        }
        if target_end - target_start == 1 {
            self.recurse_half_singleton::<WOBBLE, false>(
                target_start,
                query_start..query_end,
                depth,
                match_streak,
                mm_count,
            );
            return;
        }

        let query_partitions = partition(self.q, query_start, query_end, depth);
        let target_partitions = partition(self.t, target_start, target_end, depth);

        if query_partitions[0] == query_end || target_partitions[0] == target_end {
            return;
        }

        let next_depth = depth + 1;
        let ms = match_streak + 1;
        let can_mm = self.can_mismatch_next(depth, match_streak, mm_count);

        for query_slot in [0, 1, 2, 4] {
            if query_partitions[query_slot] >= query_partitions[query_slot + 1] {
                continue;
            }
            let query_base = unsafe { Base::from_idx(query_slot + 1) };

            for target_slot in [0, 1, 2, 4] {
                if target_partitions[target_slot] >= target_partitions[target_slot + 1] {
                    continue;
                }
                let target_base = unsafe { Base::from_idx(target_slot + 1) };

                if query_base.pair_type(target_base).is_match(WOBBLE) {
                    self.recurse::<WOBBLE>(
                        query_partitions[query_slot]..query_partitions[query_slot + 1],
                        target_partitions[target_slot]..target_partitions[target_slot + 1],
                        next_depth,
                        ms,
                        mm_count,
                    );
                } else if can_mm {
                    self.recurse::<WOBBLE>(
                        query_partitions[query_slot]..query_partitions[query_slot + 1],
                        target_partitions[target_slot]..target_partitions[target_slot + 1],
                        next_depth,
                        0,
                        mm_count + 1,
                    );
                }
            }
        }
    }

    /// Half-singleton: one SA interval has a single entry, the other has multiple.
    /// `Q_SINGLETON=true` -> query is the single entry, partition target.
    /// `Q_SINGLETON=false` -> target is the single entry, partition query.
    #[inline(always)]
    fn recurse_half_singleton<const WOBBLE: bool, const Q_SINGLETON: bool>(
        &mut self,
        singleton_sa_idx: usize,
        multi_sa: Range<usize>,
        depth: usize,
        match_streak: usize,
        mm_count: usize,
    ) {
        let (singleton_view, multi_view) = if Q_SINGLETON {
            (self.q, self.t)
        } else {
            (self.t, self.q)
        };

        let singleton_base = unsafe { singleton_view.sa_base_unchecked(singleton_sa_idx, depth) };
        if !singleton_base.is_matchable() {
            return;
        }

        let multi_partitions = partition(multi_view, multi_sa.start, multi_sa.end, depth);

        let d1 = depth + 1;
        let can_mm = self.can_mismatch_next(depth, match_streak, mm_count);
        let singleton_range = singleton_sa_idx..singleton_sa_idx + 1;
        let ms = match_streak + 1;
        let mm1 = mm_count + 1;

        for multi_slot in [0, 1, 2, 4] {
            if multi_partitions[multi_slot] >= multi_partitions[multi_slot + 1] {
                continue;
            }
            let multi_base = unsafe { Base::from_idx(multi_slot + 1) };
            let multi_range = multi_partitions[multi_slot]..multi_partitions[multi_slot + 1];
            let (query_sa, target_sa) = if Q_SINGLETON {
                (singleton_range.clone(), multi_range)
            } else {
                (multi_range, singleton_range.clone())
            };
            if singleton_base.pair_type(multi_base).is_match(WOBBLE) {
                self.recurse::<WOBBLE>(query_sa, target_sa, d1, ms, mm_count);
            } else if can_mm {
                self.recurse::<WOBBLE>(query_sa, target_sa, d1, 0, mm1);
            }
        }
    }

    /// Both-singleton fast path: tight linear scan when both SA intervals have one entry.
    ///
    /// Assumes the caller already handled emission at the current depth, then advances
    /// linearly and emits only for newly reached depths.
    #[inline(always)]
    fn recurse_singleton<const WOBBLE: bool>(
        &mut self,
        query_sa_idx: usize,
        target_sa_idx: usize,
        mut depth: usize,
        mut match_streak: usize,
        mut mm_count: usize,
    ) {
        loop {
            if depth >= self.max_len {
                return;
            }

            if !self.can_reach_suffix(depth, match_streak, mm_count) {
                return;
            }

            let d1 = depth + 1;
            let can_mm = self.can_mismatch_next(depth, match_streak, mm_count);

            let query_base = unsafe { self.q.sa_base_unchecked(query_sa_idx, depth) };
            let target_base = unsafe { self.t.sa_base_unchecked(target_sa_idx, depth) };
            if !query_base.is_matchable() || !target_base.is_matchable() {
                return;
            }

            if query_base.pair_type(target_base).is_match(WOBBLE) {
                depth = d1;
                match_streak += 1;
            } else if can_mm {
                depth = d1;
                mm_count += 1;
                match_streak = 0;
            } else {
                return;
            }

            if self.should_emit(depth, match_streak, mm_count) {
                (self.on_match)(SeedMatch {
                    query_interval: query_sa_idx..query_sa_idx + 1,
                    target_interval: target_sa_idx..target_sa_idx + 1,
                    seed_len: depth,
                });
            }
        }
    }
}

const LINEAR_PARTITION_CUTOFF: usize = 1024;

/// Partition a sorted SA interval by base character at `depth`.
///
/// Returns 6 boundary positions `[A, C, G, N, U, end]` in suffix-array order.
/// Sub-interval for slot `k` is `bounds[k]..bounds[k+1]`.
/// Callers that search only matchable RNA bases should skip the `N` bucket (slot 3).
#[inline(always)]
fn partition(view: SeedView<'_>, start: usize, end: usize, depth: usize) -> [usize; 6] {
    if start >= end {
        return [start; 6];
    }

    let mut out = [0usize; 6];
    if end - start <= LINEAR_PARTITION_CUTOFF {
        let mut i = start;
        for slot in 0..5 {
            let target = (slot + 1) as u8;
            while i < end && unsafe { view.sa_base_unchecked(i, depth) as u8 } < target {
                i += 1;
            }
            out[slot] = i;
        }
    } else {
        for slot in 0..5 {
            let target = (slot + 1) as u8;
            out[slot] = binary_search(view, start, end, depth, target);
        }
    }
    out[5] = end;
    out
}

/// Binary search for leftmost position where character at `depth` >= `target`.
#[inline(always)]
fn binary_search(
    view: SeedView<'_>,
    mut start: usize,
    mut end: usize,
    depth: usize,
    target: u8,
) -> usize {
    let mut half = (end - start) >> 1;
    while start < end {
        let mid = start + half;
        if unsafe { view.sa_base_unchecked(mid, depth) as u8 } >= target {
            end = start + half;
        } else {
            start += if half != 0 { half } else { 1 };
        }
        half >>= 1;
    }
    start
}

#[cfg(test)]
mod tests {
    use crate::config::{MismatchSpec, SeedConfig, SeedSpec};
    use crate::index::sa::SuffixArray;
    use crate::seed::SeedView;
    use crate::seq::Sequence;
    use crate::types::Base;

    use super::*;

    fn build_padded_sa(seq: &Sequence) -> (Vec<u64>, Vec<Base>, usize) {
        let sa = SuffixArray::try_from(&seq[..]).expect("SA construction failed");
        let real_len = sa.len();
        let mut padded: Vec<u64> = sa.into_inner();
        padded.resize(padded.len() + crate::index::store::SA_CHAR_PADDING, 0u64);
        let mut padded_seq: Vec<Base> = seq.iter().copied().collect();
        padded_seq.resize(
            padded_seq.len() + crate::index::store::SA_CHAR_PADDING,
            Base::Gap,
        );
        (padded, padded_seq, real_len)
    }

    #[test]
    fn partition_splits_by_base() {
        let seq_bases = vec![Base::A, Base::G, Base::C, Base::U];
        let seq = Sequence::from(seq_bases);
        let (padded_sa, padded_seq, real_len) = build_padded_sa(&seq);
        let view = SeedView {
            combined_seq: &padded_seq,
            combined_sa: &padded_sa,
            sa_real_len: real_len,
            offsets: &[],
            seq_lens: &[],
        };
        let parts = partition(view, 0, real_len, 0);

        assert_eq!(parts[1] - parts[0], 1); // A
        assert_eq!(parts[2] - parts[1], 1); // C
        assert_eq!(parts[3] - parts[2], 1); // G
        assert_eq!(parts[4] - parts[3], 0); // N
        assert_eq!(parts[5] - parts[4], 1); // U
    }

    #[test]
    fn singleton_handoff_does_not_double_emit() {
        // Distinct first symbols force singleton intervals at depth=1 in both
        // query and target branches, exercising the singleton fast-path handoff.
        let q_seq = Sequence::from(vec![Base::A, Base::C]);
        let t_seq = Sequence::from(vec![Base::U, Base::G]);
        let (q_sa, q_seq_padded, q_sa_len) = build_padded_sa(&q_seq);
        let (t_sa, t_seq_padded, t_sa_len) = build_padded_sa(&t_seq);

        let cfg = SeedConfig::with_wobble(SeedSpec::LengthOnly(1), MismatchSpec::exact(), false);
        let mut seen = std::collections::HashSet::new();
        let mut ctx = ParallelSaTraverser {
            q: SeedView {
                combined_seq: &q_seq_padded,
                combined_sa: &q_sa,
                sa_real_len: q_sa_len,
                offsets: &[],
                seq_lens: &[],
            },
            t: SeedView {
                combined_seq: &t_seq_padded,
                combined_sa: &t_sa,
                sa_real_len: t_sa_len,
                offsets: &[],
                seq_lens: &[],
            },
            min_len: 1,
            max_len: 2,
            max_mm: cfg.mismatch.max_mismatches,
            min_prefix: cfg.mismatch.min_prefix_matches,
            min_suffix: cfg.mismatch.min_suffix_matches,
            on_match: &mut |m| {
                assert!(
                    seen.insert((
                        m.query_interval.start,
                        m.query_interval.end,
                        m.target_interval.start,
                        m.target_interval.end,
                        m.seed_len
                    )),
                    "duplicate seed match emitted: {:?}",
                    m
                );
            },
        };

        let q_range = 0..ctx.q.sa_real_len;
        let t_range = 0..ctx.t.sa_real_len;
        ctx.recurse::<false>(q_range, t_range, 0, 0, 0);
    }
}
