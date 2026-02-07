//! Parallel Suffix Array Search for RNA Interaction Discovery
//!
//! This module implements a parallel SA traversal algorithm that searches both
//! query and target suffix arrays simultaneously. Following the C RIsearch2
//! algorithm, the target SA is built on the COMPLEMENT of the target sequence,
//! enabling exact character matching to find complementary RNA base pairs.
//!
//! # Algorithm Overview
//!
//! The C algorithm builds target SA on complement(target). This means:
//! - Query A matching target_comp A → original target has U → A-U pair ✓
//! - Query G matching target_comp A → original target has U → G-U wobble ✓
//!
//! ```text
//! Level 0:  [entire query SA] × [entire target_complement SA]
//!               ↓
//! Level k:  Partition both SA intervals by base (A,C,G,U,N)
//!           Recurse on:
//!             - Same character pairs (canonical base-pairing via complement)
//!             - Wobble pairs: query G with target_comp A, query U with target_comp C
//!               ↓
//! Level s:  When depth == seed_length, collect matches
//! ```

use crate::config::SeedConfig;
use crate::index::sa::SuffixArray;
use crate::seq::Sequence;
use crate::types::{Base, Interval};

const BASES: [Base; 4] = [Base::A, Base::C, Base::G, Base::U];

/// Partitioned intervals for all RNA bases within an SA range.
///
/// After partitioning, the SA interval is divided into contiguous sub-intervals
/// where all suffixes in each sub-interval start with the same base.
///
/// Base discriminant ordering: A(1) < G(2) < C(3) < U(4) < N(5)
/// Boundaries array: [a_start, g_start, c_start, u_start, n_start, end]
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct BaseIntervals {
    /// Boundaries: [a_start, g_start, c_start, u_start, n_start, end]
    /// Reflects discriminant ordering: A < G < C < U < N
    bounds: [usize; 6],
}

impl BaseIntervals {
    /// Create from raw bounds array
    /// Expected order: [a_start, g_start, c_start, u_start, n_start, end]
    pub(crate) fn from_bounds(bounds: [usize; 6]) -> Self {
        Self { bounds }
    }

    /// Get interval for a specific base
    /// Indexing: 0=A, 1=G, 2=C, 3=U, 4=N (matches discriminant order)
    #[inline]
    pub(crate) fn get(&self, base: Base) -> Interval {
        match base {
            Base::Gap => Interval::new(0, 0),
            Base::A => Interval::new(self.bounds[0], self.bounds[1]),
            Base::G => Interval::new(self.bounds[1], self.bounds[2]),
            Base::C => Interval::new(self.bounds[2], self.bounds[3]),
            Base::U => Interval::new(self.bounds[3], self.bounds[4]),
            Base::N => Interval::new(self.bounds[4], self.bounds[5]),
        }
    }


    /// Check if any base has a non-empty interval
    #[inline]
    pub(crate) fn any_non_empty(&self) -> bool {
        self.bounds[5] > self.bounds[0]
    }
}

/// A seed match found by parallel SA search
#[derive(Debug, Clone)]
pub(crate) struct SeedMatch {
    /// Interval in query SA containing matching suffixes
    pub(crate) query_interval: Interval,
    /// Interval in target SA containing matching suffixes
    pub(crate) target_interval: Interval,
    /// Seed length for this match.
    pub(crate) seed_len: usize,
}

/// Search state during parallel SA traversal
#[derive(Debug, Clone, Copy)]
struct SearchState {
    query_interval: Interval,
    target_interval: Interval,
    depth: usize,
    matches_since_mismatch: usize,
    mismatch_count: usize,
}

/// Parallel Suffix Array searcher
/// - Query SA is built on the query sequence as-is
/// - Target SA is built on the COMPLEMENT of the target sequence
/// - Same-character matching finds complementary base pairs
pub(crate) struct SeedSearcher<'a> {
    /// Query suffix array
    query_sa: &'a SuffixArray,
    /// Query sequence (for base lookup)
    query_seq: &'a Sequence,
    /// Target suffix array (built on COMPLEMENT of target)
    target_comp_sa: &'a SuffixArray,
    /// Target complement sequence
    target_comp_seq: &'a Sequence,
    /// Seed configuration (pairing, mismatch spec, etc.)
    seed_config: &'a SeedConfig,
}

impl<'a> SeedSearcher<'a> {
    /// Create a new parallel SA searcher
    ///
    /// IMPORTANT: `target_comp_sa` and `target_comp_seq` should be built on the
    /// COMPLEMENT (not reverse complement) of the target sequence.
    pub(crate) fn new(
        query_sa: &'a SuffixArray,
        query_seq: &'a Sequence,
        target_comp_sa: &'a SuffixArray,
        target_comp_seq: &'a Sequence,
        seed_config: &'a SeedConfig,
    ) -> Self {
        Self {
            query_sa,
            query_seq,
            target_comp_sa,
            target_comp_seq,
            seed_config,
        }
    }

    /// Find seeds for a range of lengths in a single traversal.
    #[inline]
    pub(crate) fn search_length_range(
        &self,
        min_len: usize,
        max_len: usize,
        results: &mut Vec<SeedMatch>,
    ) {
        let initial_state = SearchState {
            query_interval: Interval::new(0, self.query_sa.len()),
            target_interval: Interval::new(0, self.target_comp_sa.len()),
            depth: 0,
            matches_since_mismatch: 0,
            mismatch_count: 0,
        };

        self.recurse_length_range(min_len, max_len, initial_state, results);
    }

    /// Recursive parallel search over a length range (mirrors C's sa_parallel_match_neg)
    #[allow(clippy::too_many_arguments)]
    #[inline(never)] // Keep separate for flamegraph
    fn recurse_length_range(
        &self,
        min_len: usize,
        max_len: usize,
        state: SearchState,
        results: &mut Vec<SeedMatch>,
    ) {
        // Report if depth is within range
        if state.depth >= min_len
            && state.depth <= max_len
            && self.is_valid_match(&state, state.depth)
        {
            results.push(SeedMatch {
                query_interval: state.query_interval,
                target_interval: state.target_interval,
                seed_len: state.depth,
            });
        }

        if state.depth >= max_len {
            return;
        }

        // If a mismatch has occurred, ensure we can still satisfy min_suffix_matches
        if state.mismatch_count > 0 && self.seed_config.mismatch.min_suffix_matches > 0 {
            let max_possible = state.matches_since_mismatch + (max_len - state.depth);
            if max_possible < self.seed_config.mismatch.min_suffix_matches {
                return;
            }
        }

        // If we haven't reached min_len yet, ensure both intervals have at least
        // one suffix long enough to ever reach min_len.
        if state.depth < min_len
            && (!self.has_suffix_len_at_least(
                self.query_sa,
                self.query_seq,
                state.query_interval,
                min_len,
            ) || !self.has_suffix_len_at_least(
                self.target_comp_sa,
                self.target_comp_seq,
                state.target_interval,
                min_len,
            ))
        {
            return;
        }

        // Partition both SA intervals by base at current depth
        let qint = self.partition_interval(
            self.query_sa,
            self.query_seq,
            state.query_interval,
            state.depth,
        );
        let sint = self.partition_interval(
            self.target_comp_sa,
            self.target_comp_seq,
            state.target_interval,
            state.depth,
        );

        // Early exit if either side is empty
        if !qint.any_non_empty() || !sint.any_non_empty() {
            return;
        }

        let next_depth = state.depth + 1;

        // === CANONICAL MATCHES (same character = complementary base pair) ===
        self.explore_canonical_matches(&qint, &sint, next_depth, &state, min_len, max_len, results);

        // === WOBBLE PAIRS ===
        if self.seed_config.allows_wobble() {
            self.explore_wobble_matches(
                &qint, &sint, next_depth, &state, min_len, max_len, results,
            );
        }

        // === MISMATCH EXPLORATION ===
        if self.should_explore_mismatches(&state, max_len) {
            self.explore_mismatches(&qint, &sint, next_depth, &state, min_len, max_len, results);
        }
    }

    /// Explore canonical base pair matches (A-U, C-G, G-C, U-A)
    ///
    /// Always inlined since this is called on every recursion.
    /// For profiling wobble/mismatch overhead, these remain #[inline(never)].
    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    fn explore_canonical_matches(
        &self,
        qint: &BaseIntervals,
        sint: &BaseIntervals,
        next_depth: usize,
        state: &SearchState,
        min_len: usize,
        max_len: usize,
        results: &mut Vec<SeedMatch>,
    ) {
        for &base in BASES.iter() {
            let q = qint.get(base);
            let s = sint.get(base);
            if !q.is_empty() && !s.is_empty() {
                self.recurse_match(q, s, next_depth, state, min_len, max_len, results);
            }
        }
    }

    /// Explore wobble base pair matches (G-U, U-G)
    #[allow(clippy::too_many_arguments)]
    #[inline(never)] // Keep separate for flamegraph
    fn explore_wobble_matches(
        &self,
        qint: &BaseIntervals,
        sint: &BaseIntervals,
        next_depth: usize,
        state: &SearchState,
        min_len: usize,
        max_len: usize,
        results: &mut Vec<SeedMatch>,
    ) {
        // For direct matching (query_RC vs target, both NOT complemented):
        // - G-U wobble: query G (query_RC has C) pairs with target U (T)
        // - U-G wobble: query U (query_RC has A) pairs with target G

        // G-U wobble: query_RC C with target U (stored as T)
        let q_c = qint.get(Base::C);
        let s_u = sint.get(Base::U);
        if !q_c.is_empty() && !s_u.is_empty() {
            self.recurse_match(q_c, s_u, next_depth, state, min_len, max_len, results);
        }

        // U-G wobble: query_RC A with target G
        let q_a = qint.get(Base::A);
        let s_g = sint.get(Base::G);
        if !q_a.is_empty() && !s_g.is_empty() {
            self.recurse_match(q_a, s_g, next_depth, state, min_len, max_len, results);
        }
    }

    /// Recurse with a match (increment match counter)
    #[allow(clippy::too_many_arguments)]
    #[inline]
    fn recurse_match(
        &self,
        q_int: Interval,
        s_int: Interval,
        depth: usize,
        prev_state: &SearchState,
        min_len: usize,
        max_len: usize,
        results: &mut Vec<SeedMatch>,
    ) {
        let new_state = SearchState {
            query_interval: q_int,
            target_interval: s_int,
            depth,
            matches_since_mismatch: prev_state.matches_since_mismatch + 1,
            mismatch_count: prev_state.mismatch_count,
        };
        self.recurse_length_range(min_len, max_len, new_state, results);
    }

    /// Check if we should explore mismatch branches
    #[inline]
    fn should_explore_mismatches(&self, state: &SearchState, seed_len: usize) -> bool {
        let mismatch = &self.seed_config.mismatch;
        mismatch.max_mismatches > 0
            && state.mismatch_count < mismatch.max_mismatches
            && state.depth + 1 > mismatch.min_prefix_matches
            && state.matches_since_mismatch < seed_len
    }

    /// Explore mismatch branches (non-complementary pairs)
    #[allow(clippy::too_many_arguments)]
    #[inline(never)] // Keep separate for flamegraph
    fn explore_mismatches(
        &self,
        qint: &BaseIntervals,
        sint: &BaseIntervals,
        depth: usize,
        state: &SearchState,
        min_len: usize,
        max_len: usize,
        results: &mut Vec<SeedMatch>,
    ) {
        if max_len <= depth || max_len - depth < self.seed_config.mismatch.min_suffix_matches {
            return;
        }

        // For each query base, explore target bases that DON'T form valid pairs
        // This mirrors C's mismatch logic

        for &q_base in BASES.iter() {
            let q_int = qint.get(q_base);
            if q_int.is_empty() {
                continue;
            }

            for &t_base in BASES.iter() {
                // Skip if this is a valid match (canonical or wobble)
                if self.is_valid_pair(q_base, t_base) {
                    continue;
                }

                let t_int = sint.get(t_base);
                if t_int.is_empty() {
                    continue;
                }

                // This is a mismatch
                let new_state = SearchState {
                    query_interval: q_int,
                    target_interval: t_int,
                    depth,
                    matches_since_mismatch: 0, // Reset on mismatch
                    mismatch_count: state.mismatch_count + 1,
                };
                self.recurse_length_range(min_len, max_len, new_state, results);
            }
        }
    }

    /// Check if query base and target_comp base form a valid pair
    #[inline]
    fn is_valid_pair(&self, q_base: Base, t_comp_base: Base) -> bool {
        // Same character = canonical pair (because target is complemented)
        if q_base == t_comp_base {
            return true;
        }

        if self.seed_config.allows_wobble() {
            // G-U wobble: Query_RC C with Target U
            if q_base == Base::C && t_comp_base == Base::U {
                return true;
            }
            // U-G wobble: Query_RC A with Target G
            if q_base == Base::A && t_comp_base == Base::G {
                return true;
            }
        }

        false
    }

    /// Check if current state represents a valid match
    #[inline]
    fn is_valid_match(&self, state: &SearchState, seed_len: usize) -> bool {
        if state.mismatch_count == 0 {
            return true;
        }
        // With mismatches: need sufficient matches after last mismatch
        state.mismatch_count <= self.seed_config.mismatch.max_mismatches
            && state.matches_since_mismatch >= self.seed_config.mismatch.min_suffix_matches
            && state.matches_since_mismatch < seed_len
    }

    /// Partition an SA interval by base at given offset (C's sa_search_interval)
    ///
    /// Uses binary search to find boundaries where bases change.
    /// Returns intervals for [a, g, c, u, n] with end boundary.
    ///
    /// Short suffixes (pos + offset >= seq_len) must be filtered out via linear scan
    /// because they're scattered throughout the SA (sorted by earlier characters).
    fn partition_interval(
        &self,
        sa: &SuffixArray,
        seq: &Sequence,
        interval: Interval,
        offset: usize,
    ) -> BaseIntervals {
        if interval.is_empty() {
            return BaseIntervals::default();
        }

        // Step 1: Find valid suffix range
        let (valid_start, valid_end) = self.find_valid_suffix_range(sa, seq, interval, offset);

        if valid_start >= valid_end {
            return BaseIntervals::default();
        }

        // Step 2: Partition valid range by base (O(log n) binary search)
        self.partition_by_base(sa, seq, valid_start, valid_end, offset)
    }

    /// Find the range of suffixes long enough for this offset
    ///
    /// Optimized: uses two linear scans from both ends with early exit.
    /// For typical workloads, most suffixes are valid, so scans terminate quickly.
    #[inline]
    fn find_valid_suffix_range(
        &self,
        sa: &SuffixArray,
        seq: &Sequence,
        interval: Interval,
        offset: usize,
    ) -> (usize, usize) {
        let start = interval.start;
        let end = interval.end;
        let seq_len = seq.len();

        // For seed lengths (6-22bp) and typical sequences (>1000bp),
        // almost all suffixes are valid, so scan from start until we find first valid
        let mut valid_start = start;
        while valid_start < end {
            let pos = sa[valid_start] as usize;
            if pos + offset < seq_len {
                break;
            }
            valid_start += 1;
        }

        if valid_start >= end {
            return (end, end);
        }

        // Scan from end to find last valid
        let mut valid_end = end;
        while valid_end > valid_start {
            let pos = sa[valid_end - 1] as usize;
            if pos + offset < seq_len {
                break;
            }
            valid_end -= 1;
        }

        (valid_start, valid_end)
    }

    /// Check if any suffix in interval has length >= min_len
    #[inline]
    fn has_suffix_len_at_least(
        &self,
        sa: &SuffixArray,
        seq: &Sequence,
        interval: Interval,
        min_len: usize,
    ) -> bool {
        if min_len == 0 {
            return true;
        }
        let offset = min_len - 1;
        let (valid_start, valid_end) = self.find_valid_suffix_range(sa, seq, interval, offset);
        valid_start < valid_end
    }

    /// Partition a valid SA range by base character (O(log n) binary search)
    ///
    /// Assumes all suffixes in [valid_start..valid_end] have pos + offset < seq.len()
    ///
    /// Base ordering: Gap(0) < A(1) < G(2) < C(3) < U(4) < N(5)
    /// This differs from ASCII ordering: a < c < g < n < t
    #[inline]
    fn partition_by_base(
        &self,
        sa: &SuffixArray,
        seq: &Sequence,
        valid_start: usize,
        valid_end: usize,
        offset: usize,
    ) -> BaseIntervals {
        let sa_slice = &sa[valid_start..valid_end];

        // Find partition points for each base boundary
        // Base discriminant ordering: A(1) < G(2) < C(3) < U(4) < N(5)
        let a_start = valid_start;
        let g_start =
            valid_start + sa_slice.partition_point(|&idx| seq[idx as usize + offset] < Base::G);
        let c_start =
            valid_start + sa_slice.partition_point(|&idx| seq[idx as usize + offset] < Base::C);
        let u_start =
            valid_start + sa_slice.partition_point(|&idx| seq[idx as usize + offset] < Base::U);
        let n_start =
            valid_start + sa_slice.partition_point(|&idx| seq[idx as usize + offset] < Base::N);

        BaseIntervals::from_bounds([a_start, g_start, c_start, u_start, n_start, valid_end])
    }
}

// ============================================================================
// HIGH-LEVEL API
// ============================================================================
#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{MismatchSpec, SeedConfig, SeedSpec};
    use crate::index::sa::SuffixArray;

    /// Create a default SeedConfig for testing with specified wobble policy
    fn test_seed_config(allow_wobble: bool) -> SeedConfig {
        SeedConfig::with_wobble(SeedSpec::LengthOnly(6), MismatchSpec::exact(), allow_wobble)
    }

    #[test]
    fn test_base_complement() {
        use crate::types::Base;
        let seq = [Base::A, Base::C, Base::G, Base::U];
        let comp: Vec<_> = seq.iter().map(|b| b.complement()).collect();
        assert_eq!(comp, vec![Base::U, Base::G, Base::C, Base::A]);

        let seq2 = [Base::A, Base::A, Base::A, Base::A];
        let comp2: Vec<_> = seq2.iter().map(|b| b.complement()).collect();
        assert_eq!(comp2, vec![Base::U, Base::U, Base::U, Base::U]);
    }

    #[test]
    fn test_partition_basic() {
        use crate::types::Base;
        let seq_bases = vec![Base::A, Base::C, Base::G, Base::U];
        let seq = Sequence::from(seq_bases.clone());
        let sa = SuffixArray::try_from(&seq).expect("SA construction failed");
        let seed_args = test_seed_config(true);

        let searcher = SeedSearcher::new(&sa, &seq, &sa, &seq, &seed_args);

        let interval = Interval::new(0, sa.len());
        let parts = searcher.partition_interval(&sa, &seq, interval, 0);

        // Each base should have exactly one entry
        assert_eq!(parts.get(Base::A).len(), 1);
        assert_eq!(parts.get(Base::C).len(), 1);
        assert_eq!(parts.get(Base::G).len(), 1);
        assert_eq!(parts.get(Base::U).len(), 1);
    }

    #[test]
    fn test_homopolymer_debug() {
        use crate::types::Base;
        // Debug test for the aaa/ttt case
        let query_bases = vec![Base::A, Base::A, Base::A];
        let target_bases = vec![Base::U, Base::U, Base::U];

        // Step 1: Check complement
        let target_comp_bases: Vec<_> = target_bases.iter().map(|b| b.complement()).collect();
        eprintln!("Query: {:?}", query_bases);
        eprintln!("Target: {:?}", target_bases);
        eprintln!("Target complement: {:?}", target_comp_bases);
        assert_eq!(
            &target_comp_bases,
            &[Base::A, Base::A, Base::A],
            "complement(UUU) should be AAA"
        );

        let query = Sequence::from(query_bases);
        let _target = Sequence::from(target_bases);
        let target_comp = Sequence::from(target_comp_bases);

        // Step 2: Build SAs
        let q_sa = SuffixArray::try_from(&query).expect("query SA construction failed");
        let t_sa = SuffixArray::try_from(&target_comp).expect("target SA construction failed");
        eprintln!("Query SA: {:?}", q_sa);
        eprintln!("Target comp SA: {:?}", t_sa);

        // Step 3: Create searcher and find seeds
        let seed_args = test_seed_config(false);
        let searcher = SeedSearcher::new(&q_sa, &query, &t_sa, &target_comp, &seed_args);
        let mut matches = Vec::new();
        searcher.search_length_range(3, 3, &mut matches);
        eprintln!("Raw matches: {:?}", matches);

        // The key assertion
        assert!(
            !matches.is_empty(),
            "Should find at least one 3bp match for aaa vs aaa"
        );
    }
}
