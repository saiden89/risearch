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

use crate::args::SeedArgs;
use crate::types::{Base, SeedPairing};
use libsais::SuffixArrayConstruction;

/// Interval in a suffix array [start, end) - half-open range
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SaInterval {
    pub start: usize,
    pub end: usize,
}

impl SaInterval {
    #[inline]
    pub const fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    #[inline]
    pub const fn is_empty(&self) -> bool {
        self.start >= self.end
    }

    #[inline]
    pub const fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }
}

/// Partitioned intervals for all RNA bases within an SA range.
///
/// After partitioning, the SA interval is divided into contiguous sub-intervals
/// where all suffixes in each sub-interval start with the same base.
///
/// The C algorithm uses interval array [6] with boundaries for a,c,g,n,u,end
#[derive(Debug, Clone, Copy, Default)]
pub struct BaseIntervals {
    /// Boundaries: [a_start, c_start, g_start, n_start, u_start, end]
    bounds: [usize; 6],
}

impl BaseIntervals {
    /// Create from raw bounds array (C-style)
    pub fn from_bounds(bounds: [usize; 6]) -> Self {
        Self { bounds }
    }

    /// Get interval for a specific base (using C's indexing: a=0, c=1, g=2, n=3, u=4)
    #[inline]
    pub fn get(&self, base: Base) -> SaInterval {
        match base {
            Base::Gap => SaInterval::new(0, 0),
            Base::A => SaInterval::new(self.bounds[0], self.bounds[1]),
            Base::C => SaInterval::new(self.bounds[1], self.bounds[2]),
            Base::G => SaInterval::new(self.bounds[2], self.bounds[3]),
            Base::N => SaInterval::new(self.bounds[3], self.bounds[4]),
            Base::U => SaInterval::new(self.bounds[4], self.bounds[5]),
        }
    }

    /// Get interval by C-style index (0=a, 1=c, 2=g, 3=n, 4=u)
    #[inline]
    pub fn get_by_idx(&self, idx: usize) -> SaInterval {
        debug_assert!(idx < 5);
        SaInterval::new(self.bounds[idx], self.bounds[idx + 1])
    }

    /// Check if any base has a non-empty interval
    #[inline]
    pub fn any_non_empty(&self) -> bool {
        self.bounds[5] > self.bounds[0]
    }
}

/// A seed match found by parallel SA search
#[derive(Debug, Clone)]
pub struct ParallelSeedMatch {
    /// Interval in query SA containing matching suffixes
    pub query_interval: SaInterval,
    /// Interval in target SA containing matching suffixes
    pub target_interval: SaInterval,
    /// Depth at which match was found (= seed length)
    pub depth: usize,
}

impl ParallelSeedMatch {
    /// Count total number of position pairs in this match
    pub fn pair_count(&self) -> usize {
        self.query_interval.len() * self.target_interval.len()
    }
}

// MismatchSpec deleted - use MismatchSpec from crate::seed instead

/// Search state during parallel SA traversal
#[derive(Debug, Clone, Copy)]
struct SearchState {
    query_interval: SaInterval,
    target_interval: SaInterval,
    depth: usize,
    matches_since_mismatch: usize,
    mismatch_count: usize,
}

/// Parallel Suffix Array searcher (C-style algorithm)
///
/// This searcher uses the C RIsearch2 algorithm where:
/// - Query SA is built on the query sequence as-is
/// - Target SA is built on the COMPLEMENT of the target sequence
/// - Same-character matching finds complementary base pairs
pub struct ParallelSaSearcher<'a> {
    /// Query suffix array
    query_sa: &'a [u32],
    /// Query sequence (for base lookup)
    query_seq: &'a [u8],
    /// Target suffix array (built on COMPLEMENT of target)
    target_comp_sa: &'a [u32],
    /// Target complement sequence (for base lookup)
    target_comp_seq: &'a [u8],
    /// Seed configuration (pairing, mismatch spec, etc.)
    seed_config: &'a SeedArgs,
}

impl<'a> ParallelSaSearcher<'a> {
    /// Create a new parallel SA searcher
    ///
    /// IMPORTANT: `target_comp_sa` and `target_comp_seq` should be built on the
    /// COMPLEMENT (not reverse complement) of the target sequence.
    pub fn new(
        query_sa: &'a [u32],
        query_seq: &'a [u8],
        target_comp_sa: &'a [u32],
        target_comp_seq: &'a [u8],
        seed_config: &'a SeedArgs,
    ) -> Self {
        Self {
            query_sa,
            query_seq,
            target_comp_sa,
            target_comp_seq,
            seed_config,
        }
    }

    /// Find all seed matches of given length
    ///
    /// Returns all (query_interval, target_interval) pairs where suffixes
    /// form valid RNA base-pair complementary matches of length `seed_len`.
    ///
    /// Note: Short suffixes (length < seed_len) are automatically excluded during
    /// partitioning via the sentinel value (255) which sorts after all valid bases.
    pub fn find_seeds(&self, seed_len: usize) -> Vec<ParallelSeedMatch> {
        use log::trace;

        trace!(
            "[PSA] find_seeds seed_len={} q_len={} t_len={} q_seq={:?} t_seq={:?}",
            seed_len,
            self.query_seq.len(),
            self.target_comp_seq.len(),
            String::from_utf8_lossy(self.query_seq),
            String::from_utf8_lossy(self.target_comp_seq)
        );

        let mut results = Vec::new();

        let initial_state = SearchState {
            query_interval: SaInterval::new(0, self.query_sa.len()),
            target_interval: SaInterval::new(0, self.target_comp_sa.len()),
            depth: 0,
            matches_since_mismatch: 0,
            mismatch_count: 0,
        };

        self.search_recursive(seed_len, initial_state, &mut results);

        trace!("[PSA] find_seeds returning {} matches", results.len());
        results
    }

    /// Recursive parallel search (mirrors C's sa_parallel_match_neg)
    fn search_recursive(
        &self,
        seed_len: usize,
        state: SearchState,
        results: &mut Vec<ParallelSeedMatch>,
    ) {
        // Check if we've reached seed length - report match if valid
        if state.depth >= seed_len {
            if self.is_valid_match(&state, seed_len) {
                results.push(ParallelSeedMatch {
                    query_interval: state.query_interval,
                    target_interval: state.target_interval,
                    depth: state.depth,
                });
            }
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
        // C code: query 'a' matches target_comp 'a' → target has 'u' → A-U pair

        // A matches (query A with target_comp A)
        let q_a = qint.get(Base::A);
        let s_a = sint.get(Base::A);
        if !q_a.is_empty() && !s_a.is_empty() {
            self.recurse_match(q_a, s_a, next_depth, &state, seed_len, results);
        }

        // C matches (query C with target_comp C)
        let q_c = qint.get(Base::C);
        let s_c = sint.get(Base::C);
        if !q_c.is_empty() && !s_c.is_empty() {
            self.recurse_match(q_c, s_c, next_depth, &state, seed_len, results);
        }

        // G matches (query G with target_comp G)
        let q_g = qint.get(Base::G);
        let s_g = sint.get(Base::G);
        if !q_g.is_empty() && !s_g.is_empty() {
            self.recurse_match(q_g, s_g, next_depth, &state, seed_len, results);
        }

        // U matches (query U with target_comp U)
        let q_u = qint.get(Base::U);
        let s_u = sint.get(Base::U);
        if !q_u.is_empty() && !s_u.is_empty() {
            self.recurse_match(q_u, s_u, next_depth, &state, seed_len, results);
        }

        // === WOBBLE PAIRS ===
        // For direct matching (query_RC vs target, both NOT complemented):
        // - G-U wobble: query G (query_RC has C) pairs with target U (T)
        //   → match query_RC C with target T (U)
        // - U-G wobble: query U (query_RC has A) pairs with target G
        //   → match query_RC A with target G
        if matches!(self.seed_config.pairing, SeedPairing::AllowWobble) {
            // G-U wobble: query_RC C with target U (stored as T)
            if !q_c.is_empty() && !s_u.is_empty() {
                self.recurse_match(q_c, s_u, next_depth, &state, seed_len, results);
            }

            // U-G wobble: query_RC A with target G
            if !q_a.is_empty() && !s_g.is_empty() {
                self.recurse_match(q_a, s_g, next_depth, &state, seed_len, results);
            }
        }

        // === MISMATCH EXPLORATION ===
        if self.should_explore_mismatches(&state, seed_len) {
            self.explore_mismatches(&qint, &sint, next_depth, &state, seed_len, results);
        }
    }

    /// Recurse with a match (increment match counter)
    #[inline]
    fn recurse_match(
        &self,
        q_int: SaInterval,
        s_int: SaInterval,
        depth: usize,
        prev_state: &SearchState,
        seed_len: usize,
        results: &mut Vec<ParallelSeedMatch>,
    ) {
        let new_state = SearchState {
            query_interval: q_int,
            target_interval: s_int,
            depth,
            matches_since_mismatch: prev_state.matches_since_mismatch + 1,
            mismatch_count: prev_state.mismatch_count,
        };
        self.search_recursive(seed_len, new_state, results);
    }

    /// Check if we should explore mismatch branches
    #[inline]
    fn should_explore_mismatches(&self, state: &SearchState, seed_len: usize) -> bool {
        self.seed_config.mismatch_seed.max_mismatches > 0
            && state.mismatch_count < self.seed_config.mismatch_seed.max_mismatches
            && state.depth + 1 > self.seed_config.mismatch_seed.min_position
            && state.matches_since_mismatch < seed_len
    }

    /// Explore mismatch branches (non-complementary pairs)
    fn explore_mismatches(
        &self,
        qint: &BaseIntervals,
        sint: &BaseIntervals,
        depth: usize,
        state: &SearchState,
        seed_len: usize,
        results: &mut Vec<ParallelSeedMatch>,
    ) {
        // For each query base, explore target bases that DON'T form valid pairs
        // This mirrors C's mismatch logic

        let bases = [Base::A, Base::C, Base::G, Base::U];

        for &q_base in &bases {
            let q_int = qint.get(q_base);
            if q_int.is_empty() {
                continue;
            }

            for &t_base in &bases {
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
                self.search_recursive(seed_len, new_state, results);
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

        // Wobble pairs - must match canonical wobble exploration above!
        // We use query_RC, so:
        // - G-U wobble: query G (query_RC has C) with target U
        // - U-G wobble: query U (query_RC has A) with target G
        if matches!(self.seed_config.pairing, SeedPairing::AllowWobble) {
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
        state.mismatch_count <= self.seed_config.mismatch_seed.max_mismatches
            && state.matches_since_mismatch >= self.seed_config.mismatch_seed.min_matches_after
            && state.matches_since_mismatch < seed_len
    }

    /// Partition an SA interval by base at given offset (C's sa_search_interval)
    ///
    /// Uses binary search to find boundaries where bases change.
    /// Returns intervals for [a, c, g, n, u] with end boundary.
    ///
    /// Suffixes that are too short (OOB at this offset) are assigned sentinel value 255,
    /// which sorts after all valid bases and are excluded from the U interval.
    fn partition_interval(
        &self,
        sa: &[u32],
        seq: &[u8],
        interval: SaInterval,
        offset: usize,
    ) -> BaseIntervals {
        if interval.is_empty() {
            return BaseIntervals::default();
        }

        let start = interval.start;
        let end = interval.end;

        // For suffixes where pos + offset >= seq.len(), the character is undefined.
        // We need to exclude these from valid base intervals.
        //
        // Key insight: SA is sorted by suffix content, not position. Short suffixes
        // can appear BEFORE or AFTER longer ones depending on the sequence content.
        // We use a sentinel approach: OOB positions map to byte 0 (sorts before 'a'),
        // so they get excluded from the A interval start.
        //
        // Find where valid entries start (first suffix long enough for this offset)
        let valid_start = start
            + (start..end)
                .position(|i| {
                    let pos = sa[i] as usize;
                    pos + offset < seq.len()
                })
                .unwrap_or(end - start);

        // Find where valid entries end (last valid + 1)
        let valid_end = start
            + (start..end)
                .rposition(|i| {
                    let pos = sa[i] as usize;
                    pos + offset < seq.len()
                })
                .map(|p| p + 1)
                .unwrap_or(0);

        if valid_start >= valid_end {
            return BaseIntervals::default();
        }

        let sa_slice = &sa[valid_start..valid_end];

        // Now partition only the valid range
        let a_start = valid_start;
        let c_start = valid_start
            + sa_slice.partition_point(|&idx| {
                let c = seq[idx as usize + offset];
                c < b'c'
            });
        let g_start = valid_start
            + sa_slice.partition_point(|&idx| {
                let c = seq[idx as usize + offset];
                c < b'g'
            });
        let n_start = valid_start
            + sa_slice.partition_point(|&idx| {
                let c = seq[idx as usize + offset];
                c < b'n'
            });
        let u_start = valid_start
            + sa_slice.partition_point(|&idx| {
                let c = seq[idx as usize + offset];
                c < b't'
            });

        BaseIntervals::from_bounds([a_start, c_start, g_start, n_start, u_start, valid_end])
    }
}

// ============================================================================
// HIGH-LEVEL API
// ============================================================================
/// Compute the complement of a sequence (A<->U, C<->G)
pub fn complement_sequence(seq: &[u8]) -> Vec<u8> {
    use crate::types::COMPLEMENT;
    seq.iter().map(|&b| COMPLEMENT[b as usize]).collect()
}

/// Build suffix array for a sequence
pub fn build_suffix_array(seq: &[u8]) -> Vec<u32> {
    SuffixArrayConstruction::for_text(seq)
        .in_owned_buffer()
        .single_threaded()
        .run()
        .expect("SA construction should not fail for valid sequences")
        .into_vec()
        .into_iter()
        .map(|x: i64| x as u32)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::SeedArgs;
    use crate::seed::{MismatchSpec, SeedSpec};

    /// Create a default SeedArgs for testing with specified pairing mode
    fn test_seed_args(pairing: SeedPairing) -> SeedArgs {
        SeedArgs {
            seed: SeedSpec::Length(6),
            pairing,
            mismatch_seed: MismatchSpec::exact(),
        }
    }

    #[test]
    fn test_complement_sequence() {
        assert_eq!(complement_sequence(b"acgt"), b"tgca");
        assert_eq!(complement_sequence(b"ACGU"), b"tgca");
        assert_eq!(complement_sequence(b"aaaa"), b"tttt");
    }

    #[test]
    fn test_partition_basic() {
        let seq = b"acgt";
        let sa = build_suffix_array(seq);
        let seed_args = test_seed_args(SeedPairing::AllowWobble);

        let searcher = ParallelSaSearcher::new(&sa, seq, &sa, seq, &seed_args);

        let interval = SaInterval::new(0, sa.len());
        let parts = searcher.partition_interval(&sa, seq, interval, 0);

        // Each base should have exactly one entry
        assert_eq!(parts.get(Base::A).len(), 1);
        assert_eq!(parts.get(Base::C).len(), 1);
        assert_eq!(parts.get(Base::G).len(), 1);
        assert_eq!(parts.get(Base::U).len(), 1);
    }

    #[test]
    fn test_homopolymer_debug() {
        // Debug test for the aaa/ttt case
        let query = b"aaa";
        let target = b"ttt";

        // Step 1: Check complement
        let target_comp = complement_sequence(target);
        eprintln!("Query: {:?}", String::from_utf8_lossy(query));
        eprintln!("Target: {:?}", String::from_utf8_lossy(target));
        eprintln!(
            "Target complement: {:?}",
            String::from_utf8_lossy(&target_comp)
        );
        assert_eq!(&target_comp, b"aaa", "complement(ttt) should be aaa");

        // Step 2: Build SAs
        let q_sa = build_suffix_array(query);
        let t_sa = build_suffix_array(&target_comp);
        eprintln!("Query SA: {:?}", q_sa);
        eprintln!("Target comp SA: {:?}", t_sa);

        // Step 3: Create searcher and find seeds
        let seed_args = test_seed_args(SeedPairing::Strict);
        let searcher = ParallelSaSearcher::new(&q_sa, query, &t_sa, &target_comp, &seed_args);
        let matches = searcher.find_seeds(3);
        eprintln!("Raw matches: {:?}", matches);

        // The key assertion
        assert!(
            !matches.is_empty(),
            "Should find at least one 3bp match for aaa vs aaa"
        );
    }
}
