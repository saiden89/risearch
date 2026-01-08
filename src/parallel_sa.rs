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

use crate::seed::MismatchSpec;
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
    /// Whether to allow G-U wobble pairs in seeds
    allow_wobble: bool,
    /// Mismatch configuration
    mismatch_config: MismatchSpec,
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
        pairing: SeedPairing,
    ) -> Self {
        Self {
            query_sa,
            query_seq,
            target_comp_sa,
            target_comp_seq,
            allow_wobble: matches!(pairing, SeedPairing::AllowWobble),
            mismatch_config: MismatchSpec::exact(),
        }
    }

    /// Set mismatch configuration
    pub fn with_mismatches(mut self, config: MismatchSpec) -> Self {
        self.mismatch_config = config;
        self
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
        if self.allow_wobble {
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
        self.mismatch_config.max_mismatches > 0
            && state.mismatch_count < self.mismatch_config.max_mismatches
            && state.depth + 1 > self.mismatch_config.min_position
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

        // Wobble pairs
        if self.allow_wobble {
            // Query G with target_comp A (target has U) → G-U wobble
            if q_base == Base::G && t_comp_base == Base::A {
                return true;
            }
            // Query U with target_comp C (target has G) → U-G wobble
            if q_base == Base::U && t_comp_base == Base::C {
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
        state.mismatch_count <= self.mismatch_config.max_mismatches
            && state.matches_since_mismatch >= self.mismatch_config.min_matches_after
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
///
/// This produces the sequence needed for the target SA in C-style parallel search.
pub fn complement_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'a' | b'A' => b't',               // A -> U (stored as t)
            b't' | b'T' | b'u' | b'U' => b'a', // U/T -> A
            b'c' | b'C' => b'g',               // C -> G
            b'g' | b'G' => b'c',               // G -> C
            b'n' | b'N' => b'n',               // N stays N
            _ => b'n',                         // Unknown -> N
        })
        .collect()
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

/// High-level seed finder that handles SA construction
///
/// This is the main entry point for finding all complementary seeds between
/// a query and target sequence.
pub struct SeedFinder {
    /// Whether to allow G-U wobble pairs
    pub allow_wobble: bool,
    /// Mismatch configuration
    pub mismatch_config: MismatchSpec,
}

impl Default for SeedFinder {
    fn default() -> Self {
        Self {
            allow_wobble: true,
            mismatch_config: MismatchSpec::exact(),
        }
    }
}

impl SeedFinder {
    pub fn new(pairing: SeedPairing) -> Self {
        Self {
            allow_wobble: matches!(pairing, SeedPairing::AllowWobble),
            mismatch_config: MismatchSpec::exact(),
        }
    }

    pub fn with_mismatches(mut self, config: MismatchSpec) -> Self {
        self.mismatch_config = config;
        self
    }

    /// Find all seeds between query and target
    ///
    /// Returns (query_pos, target_pos) tuples for all matches.
    pub fn find_all_seeds(
        &self,
        query_seq: &[u8],
        target_seq: &[u8],
        seed_len: usize,
    ) -> Vec<(usize, usize)> {
        // Build query SA
        let query_sa = build_suffix_array(query_seq);

        // Build target complement and its SA
        let target_comp = complement_sequence(target_seq);
        let target_comp_sa = build_suffix_array(&target_comp);

        // Create searcher and find matches
        let pairing = if self.allow_wobble {
            SeedPairing::AllowWobble
        } else {
            SeedPairing::Strict
        };

        let searcher =
            ParallelSaSearcher::new(&query_sa, query_seq, &target_comp_sa, &target_comp, pairing)
                .with_mismatches(self.mismatch_config);

        let matches = searcher.find_seeds(seed_len);

        // Expand matches to concrete positions
        // Filter out positions where suffix is shorter than seed_len
        let mut positions = Vec::new();
        for m in &matches {
            for qi in m.query_interval.start..m.query_interval.end {
                let q_pos = query_sa[qi] as usize;
                // Skip if query suffix too short
                if q_pos + seed_len > query_seq.len() {
                    continue;
                }
                for ti in m.target_interval.start..m.target_interval.end {
                    let t_pos = target_comp_sa[ti] as usize;
                    // Skip if target suffix too short
                    if t_pos + seed_len > target_comp.len() {
                        continue;
                    }
                    positions.push((q_pos, t_pos));
                }
            }
        }

        positions
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

        let searcher = ParallelSaSearcher::new(&sa, seq, &sa, seq, SeedPairing::AllowWobble);

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
        let searcher =
            ParallelSaSearcher::new(&q_sa, query, &t_sa, &target_comp, SeedPairing::Strict);
        let matches = searcher.find_seeds(3);
        eprintln!("Raw matches: {:?}", matches);

        // The key assertion
        assert!(
            !matches.is_empty(),
            "Should find at least one 3bp match for aaa vs aaa"
        );
    }

    #[test]
    fn test_find_au_pair() {
        // Query has A at position 0
        // Target has U at position 0
        // With complement: target_comp has A at position 0
        // Query A should match target_comp A
        let query = b"a";
        let target = b"t"; // U stored as t

        let finder = SeedFinder::new(SeedPairing::Strict);
        let matches = finder.find_all_seeds(query, target, 1);

        assert_eq!(matches.len(), 1, "Should find A-U pair");
        assert_eq!(matches[0], (0, 0));
    }

    #[test]
    fn test_find_gc_pair() {
        let query = b"g";
        let target = b"c";

        let finder = SeedFinder::new(SeedPairing::Strict);
        let matches = finder.find_all_seeds(query, target, 1);

        assert_eq!(matches.len(), 1, "Should find G-C pair");
    }

    #[test]
    fn test_wobble_gu() {
        let query = b"g";
        let target = b"t"; // U

        // With wobble
        let finder_wobble = SeedFinder::new(SeedPairing::AllowWobble);
        let matches_wobble = finder_wobble.find_all_seeds(query, target, 1);

        // Without wobble
        let finder_strict = SeedFinder::new(SeedPairing::Strict);
        let matches_strict = finder_strict.find_all_seeds(query, target, 1);

        assert_eq!(matches_wobble.len(), 1, "Wobble should find G-U");
        assert_eq!(matches_strict.len(), 0, "Strict should NOT find G-U");
    }

    #[test]
    fn test_longer_seed() {
        let query = b"acgt"; // ACGU
        let target = b"tgca"; // This will be complemented to ACGT

        let finder = SeedFinder::new(SeedPairing::Strict);
        let matches = finder.find_all_seeds(query, target, 4);

        // Query ACGU at pos 0 should match target TGCA (which pairs as ACGU)
        assert!(!matches.is_empty(), "Should find 4bp complementary seed");
    }

    #[test]
    fn test_multiple_matches() {
        // Query has AA at positions 0-1 and 2-3
        let query = b"aaaa";
        // Target has UU (TT) at multiple positions
        let target = b"tttt";

        let finder = SeedFinder::new(SeedPairing::Strict);
        let matches = finder.find_all_seeds(query, target, 2);

        // Should find multiple A-U seed matches
        assert!(matches.len() > 1, "Should find multiple 2bp seeds");
    }

    #[test]
    fn test_sa_interval() {
        let interval = SaInterval::new(5, 10);
        assert_eq!(interval.len(), 5);
        assert!(!interval.is_empty());

        let empty = SaInterval::new(5, 5);
        assert!(empty.is_empty());
        assert_eq!(empty.len(), 0);
    }

    #[test]
    fn test_no_match() {
        // Query A, Target A -> complement T
        // Query A doesn't pair with Target A
        let query = b"a";
        let target = b"a"; // complement is T, query A != target_comp T

        let finder = SeedFinder::new(SeedPairing::Strict);
        let matches = finder.find_all_seeds(query, target, 1);

        assert!(matches.is_empty(), "A should not pair with A");
    }

    // =========================================================================
    // GROUND TRUTH PARAMETERIZED TESTS (rstest)
    // =========================================================================
    //
    // These tests verify seed finding correctness using hardcoded expected results.
    // No C comparison needed - the expected positions are the ground truth.

    use rstest::rstest;

    /// Watson-Crick complement (no reverse) for parallel SA matching.
    /// The algorithm uses same-character matching on complement:
    /// Query A matches target_comp A → target has U → A-U pair
    fn wc_complement(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .map(|&c| match c {
                b'a' | b'A' => b't',
                b't' | b'T' | b'u' | b'U' => b'a',
                b'c' | b'C' => b'g',
                b'g' | b'G' => b'c',
                _ => c,
            })
            .collect()
    }

    /// Test all 16 canonical dinucleotide seed combinations.
    /// Query XY should find seed at (0,0) when target is complement of XY.
    #[rstest]
    fn test_dinuc_canonical(
        #[values(b'a', b'c', b'g', b't')] b1: u8,
        #[values(b'a', b'c', b'g', b't')] b2: u8,
    ) {
        let query = vec![b1, b2];
        let target = wc_complement(&query);

        let finder = SeedFinder::new(SeedPairing::Strict);
        let seeds = finder.find_all_seeds(&query, &target, 2);

        // Ground truth: complementary 2bp query/target should find exactly one seed at (0,0)
        assert_eq!(
            seeds.len(),
            1,
            "Dinuc {:?}{:?} should find exactly 1 seed, got {:?}",
            b1 as char,
            b2 as char,
            seeds
        );
        assert_eq!(
            seeds[0],
            (0, 0),
            "Dinuc {:?}{:?} seed should be at (0,0)",
            b1 as char,
            b2 as char
        );
    }

    /// Test all 64 trinucleotide combinations (3bp seeds).
    #[rstest]
    fn test_trinuc_canonical(
        #[values(b'a', b'c', b'g', b't')] b1: u8,
        #[values(b'a', b'c', b'g', b't')] b2: u8,
        #[values(b'a', b'c', b'g', b't')] b3: u8,
    ) {
        let query = vec![b1, b2, b3];
        let target = wc_complement(&query);

        let finder = SeedFinder::new(SeedPairing::Strict);
        let seeds = finder.find_all_seeds(&query, &target, 3);

        assert_eq!(
            seeds.len(),
            1,
            "Trinuc {:?}{:?}{:?} should find exactly 1 seed",
            b1 as char,
            b2 as char,
            b3 as char
        );
        assert_eq!(seeds[0], (0, 0));
    }

    /// Test non-complementary pairs should find NO seeds.
    #[rstest]
    fn test_non_complementary(#[values(b'a', b'c', b'g', b't')] base: u8) {
        // Same base on both sides = not complementary
        let query = vec![base, base];
        let target = vec![base, base]; // NOT the complement

        let finder = SeedFinder::new(SeedPairing::Strict);
        let seeds = finder.find_all_seeds(&query, &target, 2);

        assert!(
            seeds.is_empty(),
            "Non-complementary {:?}{:?} should find 0 seeds, got {:?}",
            base as char,
            base as char,
            seeds
        );
    }

    /// Test G-U wobble pairs are found with AllowWobble mode.
    #[rstest]
    fn test_wobble_gu_pairs(#[values((b'g', b't'), (b't', b'g'))] pair: (u8, u8)) {
        let (q, t) = pair;
        let query = vec![q];
        let target = vec![t];

        // With wobble enabled
        let finder_wobble = SeedFinder::new(SeedPairing::AllowWobble);
        let seeds_wobble = finder_wobble.find_all_seeds(&query, &target, 1);
        assert_eq!(seeds_wobble.len(), 1, "Wobble mode should find G-U pair");

        // Without wobble - should NOT find
        let finder_strict = SeedFinder::new(SeedPairing::Strict);
        let seeds_strict = finder_strict.find_all_seeds(&query, &target, 1);
        assert!(
            seeds_strict.is_empty(),
            "Strict mode should NOT find G-U pair"
        );
    }

    /// Test multiple seeds in longer sequences.
    /// Query "aa" in target "tttt" should find seeds at multiple positions.
    #[test]
    fn test_multiple_seed_positions() {
        let query = b"aa";
        let target = b"tttt"; // Complement is "aaaa", has 3 positions for 2bp seed

        let finder = SeedFinder::new(SeedPairing::Strict);
        let mut seeds = finder.find_all_seeds(query, target, 2);
        seeds.sort();

        // Ground truth: query "aa" at pos 0, target positions 0, 1, 2
        assert_eq!(seeds.len(), 3, "Should find 3 seed positions");
        assert_eq!(seeds, vec![(0, 0), (0, 1), (0, 2)]);
    }
}
