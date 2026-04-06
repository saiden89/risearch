//! Parallel Suffix Array seed search.
//!
//! Recursive traversal of two suffix arrays with base-pairing constraints.
//! Singleton fast-paths avoid partition overhead when one or both SA intervals
//! have a single entry.

use crate::config::SeedConfig;
use crate::types::Base;
use std::ops::Range;

mod partition;
mod singleton;

use partition::partition;
use singleton::{recurse_half_singleton, recurse_singleton};

/// The four matchable RNA bases indexed alongside SLOTS.
const BASES: [Base; 4] = [Base::A, Base::C, Base::G, Base::U];

/// Partition slot for each matchable base in BASES.
/// `partition()` returns `[A=0, C=1, G=2, N=3, U=4, end]`; this mapping
/// intentionally skips the `N` bucket and addresses only searchable A/C/G/U slots.
const SLOTS: [usize; 4] = [0, 1, 2, 4];

pub(crate) trait SeedSaView: Copy {
    fn sa_real_len(&self) -> usize;
    fn sa_suffix_pos(&self, sa_idx: usize) -> usize;
    fn sa_base(&self, sa_idx: usize, offset: usize) -> Base;
}

impl SeedSaView for (&[u64], &[Base], usize) {
    #[inline(always)]
    fn sa_real_len(&self) -> usize {
        self.2
    }

    #[inline(always)]
    fn sa_suffix_pos(&self, sa_idx: usize) -> usize {
        unsafe { *self.0.get_unchecked(sa_idx) as usize }
    }

    #[inline(always)]
    fn sa_base(&self, sa_idx: usize, offset: usize) -> Base {
        let suffix_pos = self.sa_suffix_pos(sa_idx);
        unsafe { *self.1.get_unchecked(suffix_pos + offset) }
    }
}

/// A seed match found by parallel SA search.
#[derive(Debug, Clone)]
pub struct SeedMatch {
    pub(crate) query_interval: Range<usize>,
    pub(crate) target_interval: Range<usize>,
    pub(crate) seed_len: usize,
}

/// Parallel SA searcher.
///
/// Query SA is built on the query as-is.
/// Target SA is built on the COMPLEMENT of the target.
/// Seed matching uses RNA pairing classes against this complemented target alphabet.
pub struct SeedSearcher<'a> {
    q: (&'a [u64], &'a [Base], usize),
    q_sa_start: usize,
    t: (&'a [u64], &'a [Base], usize),
    cfg: &'a SeedConfig,
}

impl<'a> SeedSearcher<'a> {
    pub fn new(
        q: (&'a [u64], &'a [Base], usize),
        q_sa_start: usize,
        t: (&'a [u64], &'a [Base], usize),
        cfg: &'a SeedConfig,
    ) -> Self {
        Self {
            q,
            q_sa_start,
            t,
            cfg,
        }
    }

    pub fn for_each_length_range<F>(&self, min_len: usize, max_len: usize, mut on_match: F)
    where
        F: FnMut(SeedMatch),
    {
        let mut ctx = SeedingContext {
            q: self.q,
            t: self.t,
            min_len,
            max_len,
            max_mm: self.cfg.mismatch.max_mismatches,
            min_prefix: self.cfg.mismatch.min_prefix_matches,
            min_suffix: self.cfg.mismatch.min_suffix_matches,
            on_match: &mut on_match,
        };

        let q = self.q_sa_start..self.q.sa_real_len();
        let s = 0..self.t.sa_real_len();
        if self.cfg.seed_wobble {
            recurse::<_, true>(&mut ctx, q, s, 0, 0, 0);
        } else {
            recurse::<_, false>(&mut ctx, q, s, 0, 0, 0);
        }
    }

    pub fn search_length_range(
        &self,
        min_len: usize,
        max_len: usize,
        results: &mut Vec<SeedMatch>,
    ) {
        self.for_each_length_range(min_len, max_len, |m| results.push(m));
    }
}

pub(super) struct SeedingContext<'a, F: FnMut(SeedMatch)> {
    pub(super) q: (&'a [u64], &'a [Base], usize),
    pub(super) t: (&'a [u64], &'a [Base], usize),
    min_len: usize,
    max_len: usize,
    max_mm: usize,
    min_prefix: usize,
    min_suffix: usize,
    pub(super) on_match: &'a mut F,
}

impl<F: FnMut(SeedMatch)> SeedingContext<'_, F> {
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
}

fn recurse<F: FnMut(SeedMatch), const WOBBLE: bool>(
    ctx: &mut SeedingContext<'_, F>,
    q: Range<usize>,
    s: Range<usize>,
    depth: usize,
    match_streak: usize,
    mm_count: usize,
) {
    let (ql, qr) = (q.start, q.end);
    let (sl, sr) = (s.start, s.end);

    if ctx.should_emit(depth, match_streak, mm_count) {
        (ctx.on_match)(SeedMatch {
            query_interval: ql..qr,
            target_interval: sl..sr,
            seed_len: depth,
        });
    }

    if depth >= ctx.max_len {
        return;
    }

    if !ctx.can_reach_suffix(depth, match_streak, mm_count) {
        return;
    }

    if qr - ql == 1 && sr - sl == 1 {
        recurse_singleton::<F, WOBBLE>(ctx, ql, sl, depth, match_streak, mm_count);
        return;
    }
    if qr - ql == 1 {
        recurse_half_singleton::<F, WOBBLE, true>(ctx, ql, sl..sr, depth, match_streak, mm_count);
        return;
    }
    if sr - sl == 1 {
        recurse_half_singleton::<F, WOBBLE, false>(ctx, sl, ql..qr, depth, match_streak, mm_count);
        return;
    }

    let qi = partition(ctx.q, ql, qr, depth);
    let si = partition(ctx.t, sl, sr, depth);

    if qi[0] == qr || si[0] == sr {
        return;
    }

    let d1 = depth + 1;
    let ms = match_streak + 1;
    let can_mm = ctx.can_mismatch_next(depth, match_streak, mm_count);

    for i in 0..4 {
        let qs = SLOTS[i];
        if qi[qs] >= qi[qs + 1] {
            continue;
        }

        for j in 0..4 {
            let ts = SLOTS[j];
            if si[ts] >= si[ts + 1] {
                continue;
            }

            if BASES[i].pair_type(BASES[j]).is_match(WOBBLE) {
                recurse::<F, WOBBLE>(
                    ctx,
                    qi[qs]..qi[qs + 1],
                    si[ts]..si[ts + 1],
                    d1,
                    ms,
                    mm_count,
                );
            } else if can_mm {
                recurse::<F, WOBBLE>(
                    ctx,
                    qi[qs]..qi[qs + 1],
                    si[ts]..si[ts + 1],
                    d1,
                    0,
                    mm_count + 1,
                );
            }
        }
    }
}

#[cfg(test)]
mod tests;
