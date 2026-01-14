use std::str::FromStr;

use crate::types::Strand;

/// A candidate seed match found during suffix array search.
#[derive(Debug, Clone)]
pub struct SeedCandidate {
    /// Position in query sequence (0-based)
    pub query_pos: usize,
    /// Index of target sequence in the index
    pub target_idx: usize,
    /// Start position in target sequence (0-based)
    pub target_start: usize,
    /// Length of the seed match
    pub len: usize,
    /// Strand of the match (forward or reverse)
    pub strand: Strand,
}

/// Representation of the `-m` (mismatch) flag.
///
/// Format: `c:p` where:
/// - `c` = max number of mismatches allowed in seed
/// - `p` = min number of consecutive matches required at seed start/end
///
/// Example: `-m 1:3` allows 1 mismatch with 3 consecutive matches at ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MismatchSpec {
    /// Maximum number of mismatches allowed in seed
    pub max_mismatches: usize,
    /// Minimum position (1-indexed) where mismatch can occur
    pub min_position: usize,
    /// Minimum consecutive matches required after last mismatch
    pub min_matches_after: usize,
}

impl MismatchSpec {
    /// No mismatches allowed (exact matching)
    pub const fn exact() -> Self {
        Self {
            max_mismatches: 0,
            min_position: 1,
            min_matches_after: 0,
        }
    }

    /// Create with specific parameters
    pub const fn new(max: usize, min_pos: usize, min_after: usize) -> Self {
        Self {
            max_mismatches: max,
            min_position: min_pos,
            min_matches_after: min_after,
        }
    }
}

impl FromStr for MismatchSpec {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        if s.is_empty() {
            return Err("empty mismatch spec".into());
        }

        let parts: Vec<&str> = s.split(':').collect();

        let (max_mismatches, min_consecutive) = match parts.len() {
            1 => {
                // Allow "c" as shorthand for "c:0" (matches C behavior)
                let max = parts[0]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid max mismatches: {}", e))?;
                (max, 0)
            }
            2 => {
                let max = parts[0]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid max mismatches: {}", e))?;
                let min = parts[1]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid min consecutive: {}", e))?;
                (max, min)
            }
            _ => {
                return Err(format!(
                    "invalid mismatch spec '{}': expected 'c' or 'c:p' format",
                    s
                ));
            }
        };

        // CLI format c:p maps to: max=c, min_position=p, min_after=p
        Ok(MismatchSpec {
            max_mismatches,
            min_position: min_consecutive,
            min_matches_after: min_consecutive,
        })
    }
}

/// Representation of the `-s` flag:
/// - `-s l`              => SeedSpec::Length(l)
/// - `-s m:n`            => SeedSpec::Interval { start: m, end: n, length: None }
/// - `-s m:n/l`          => SeedSpec::Interval { start: m, end: n, length: Some(l) }
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeedSpec {
    Length(i64),
    Interval {
        start: i64, // TODO: consider using RangeInclusive<i64>
        end: i64,
        length: Option<i64>, // TODO: enforce strictly positive via newtype
    },
}

impl FromStr for SeedSpec {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        if s.is_empty() {
            return Err("empty seed spec".into());
        }

        // Parse helper for i64 with context
        let parse_i64 = |val: &str, field: &str| -> Result<i64, String> {
            val.parse::<i64>()
                .map_err(|e| format!("bad {}: {}", field, e))
        };

        match s.find(':') {
            None => {
                // Length-only format: "l"
                let len = parse_i64(s, "length")?;
                Ok(SeedSpec::Length(len))
            }
            Some(colon_idx) => {
                let (start_str, rest) = s.split_at(colon_idx);
                let rest = &rest[1..]; // Skip the colon

                if rest.is_empty() {
                    return Err("missing end in interval".into());
                }

                let start = parse_i64(start_str, "start")?;

                // Check for optional length after '/'
                match rest.find('/') {
                    None => {
                        // Format: "start:end"
                        let end = parse_i64(rest, "end")?;
                        Ok(SeedSpec::Interval {
                            start,
                            end,
                            length: None,
                        })
                    }
                    Some(slash_idx) => {
                        // Format: "start:end/length"
                        let (end_str, len_part) = rest.split_at(slash_idx);
                        let len_str = &len_part[1..]; // Skip the slash

                        if end_str.is_empty() {
                            return Err("missing end in interval".into());
                        }
                        if len_str.is_empty() {
                            return Err("missing length after '/'".into());
                        }

                        let end = parse_i64(end_str, "end")?;
                        let length = parse_i64(len_str, "length")?;

                        Ok(SeedSpec::Interval {
                            start,
                            end,
                            length: Some(length),
                        })
                    }
                }
            }
        }
    }
}

impl SeedSpec {
    /// Normalize to canonical form: (start, end, length)
    ///
    /// - Uses 1-based indexing for start/end (to match the original C implementation).
    /// - query_len is the total length (n). Returned indices are in 1..=n.
    /// - Returns Err(String) for invalid specs.
    pub fn normalize(&self, query_len: usize) -> Result<(usize, usize, usize), String> {
        let n = query_len as i64;
        if n <= 0 {
            return Err("query length must be positive".into());
        }

        match *self {
            SeedSpec::Length(l) if l <= 0 => Err("Invalid seed length".into()),
            SeedSpec::Length(l) => {
                let length = l.min(n) as usize;
                Ok((1, query_len, length))
            }
            SeedSpec::Interval { start, end, length } => {
                // Validate sign consistency
                let (s_pos, e_pos) = match (start.signum(), end.signum()) {
                    // Both positive
                    (1, 1) if end >= start && end <= n => (start as usize, end as usize),
                    // Both negative
                    (-1, -1) if end >= start && start >= -n => {
                        let to_pos = |v: i64| -> usize {
                            let x = v + 1; // Adjust negative index
                            let r = (x + n) % n;
                            if r == 0 { n as usize } else { r as usize }
                        };
                        (to_pos(start), to_pos(end))
                    }
                    // Mixed signs or invalid
                    _ => return Err("Invalid seed interval: mixed sign".into()),
                };

                // Validate positions
                if s_pos == 0 || e_pos == 0 {
                    return Err("Invalid seed interval".into());
                }

                let interval_len = e_pos.saturating_sub(s_pos) + 1;
                if interval_len == 0 {
                    return Err("Invalid seed interval: empty".into());
                }

                // Determine final length
                let final_len = match length {
                    Some(l) if l < 0 => return Err("Invalid seed length".into()),
                    Some(l) if (l as usize) > interval_len => {
                        return Err("Invalid seed length (exceeds interval)".into());
                    }
                    Some(l) => l as usize,
                    None => interval_len,
                };

                Ok((s_pos, e_pos, final_len))
            }
        }
    }
}

// =============================================================================
// SEED ALIGNMENT BUILDING
// =============================================================================

use crate::seq::Seq;
use crate::types::Pairing;

/// Build alignment for the seed region.
///
/// Creates a vector of Pairing entries representing the base-pair interactions
/// in the seed region. Used by both extension and seed-only paths.
///
/// # Arguments
/// * `query` - Query sequence wrapper
/// * `target` - Target sequence wrapper
/// * `q_pos` - Starting position in query (0-based)
/// * `t_match_end` - Ending position in target (0-based, antiparallel)
/// * `len` - Length of seed
pub fn build_seed_alignment(
    query: &Seq<'_>,
    target: &Seq<'_>,
    q_pos: usize,
    t_match_end: usize,
    len: usize,
) -> smallvec::SmallVec<[Pairing; 64]> {
    let mut seed_alignment = smallvec::SmallVec::new();
    for n in 0..len {
        let q_idx = q_pos + n;
        let t_idx = t_match_end.saturating_sub(n);
        let q_b = query.base(q_idx);
        let t_b = target.base(t_idx);
        seed_alignment.push(Pairing::from_bases(q_b, t_b));
    }
    seed_alignment
}

// =============================================================================
// SEED SEARCH API - Two-level design for parallelization
// =============================================================================

use crate::args::SeedArgs;
use crate::parallel_sa::{ParallelSaSearcher, ParallelSeedMatch, build_suffix_array};
use crate::sa::{SaIndexFile, SequenceIndex};
use crate::seq::reverse_complement_dna;

/// Pre-computed query data to avoid rebuilding per-target.
pub struct QueryPrep {
    /// Normalized query (lowercase DNA, U->T)
    pub q_norm: Vec<u8>,
    /// Reverse complement of normalized query
    pub q_rc: Vec<u8>,
    /// Suffix array of reverse complement
    pub q_rc_sa: Vec<u32>,
    /// Seed interval start (0-based)
    pub start0: usize,
    /// Seed interval end (1-based, exclusive)
    pub end1: usize,
    /// Minimum seed length
    pub mi_len: usize,
    /// Bitmap: has_n[i] = true if position i contains 'n'
    /// Used for O(1) N-checking instead of O(seed_len) scan
    has_n: Vec<bool>,
}

impl QueryPrep {
    /// Build query preprocessing data (called once per query)
    pub fn new(query: &[u8], config: &SeedArgs) -> Option<Self> {
        let q_len = query.len();

        // Normalize query to lowercase DNA (t not u)
        let q_norm: Vec<u8> = query
            .iter()
            .map(|&b| {
                let lower = b.to_ascii_lowercase();
                if lower == b'u' { b't' } else { lower }
            })
            .collect();

        // Precompute N positions for O(1) checking
        let has_n: Vec<bool> = q_norm.iter().map(|&b| b == b'n').collect();

        // Get seed interval bounds
        let (start1, end1, mi_len) = config.seed.normalize(q_len).ok()?;

        // Build query RC and its SA once
        let q_rc = reverse_complement_dna(&q_norm);
        let q_rc_sa = build_suffix_array(&q_rc);

        Some(Self {
            q_norm,
            q_rc,
            q_rc_sa,
            start0: start1 - 1,
            end1,
            mi_len,
            has_n,
        })
    }

    /// Check if any position in range [start, start+len) contains 'n'
    #[inline]
    pub fn contains_n(&self, start: usize, len: usize) -> bool {
        self.has_n[start..start + len].iter().any(|&x| x)
    }
}

/// Find all seed matches between a query and a single target sequence.
///
/// This is the low-level, parallelizable unit of work. Callers can use Rayon
/// to parallelize over queries or targets as needed.
///
/// Appends candidates to the provided Vec (avoids allocation per target).
pub fn find_seeds_in_target_into(
    prep: &QueryPrep,
    target: &SequenceIndex,
    target_idx: usize,
    config: &SeedArgs,
    candidates: &mut Vec<SeedCandidate>,
    matches: &mut Vec<ParallelSeedMatch>,
) {
    let q_len = prep.q_norm.len();
    let start0 = prep.start0;
    let end1 = prep.end1;
    let mi_len = prep.mi_len;

    // Forward strand: use pre-built forward_sa
    let t_sa = &target.forward_sa;
    // Iterate all seed lengths from mi_len to q_len; position filter handles interval
    for seed_len in mi_len..=q_len {
        matches.clear(); // Reuse allocation
        let searcher = ParallelSaSearcher::new(&prep.q_rc_sa, &prep.q_rc, t_sa, &target.sequence, config);
        searcher.find_seeds_into(seed_len, matches);

        for m in matches.iter() {
            for &q_rc_pos_i32 in &prep.q_rc_sa[m.query_interval.start..m.query_interval.end] {
                let q_rc_pos = q_rc_pos_i32 as usize;
                if q_rc_pos + seed_len > q_len {
                    continue;
                }
                let q_pos = q_len - q_rc_pos - seed_len;
                // Check: seed must start >= start0 AND end <= end0
                // seed_end = q_pos + seed_len - 1, so check q_pos + seed_len <= end0 + 1 = end1
                if q_pos < start0 || q_pos + seed_len > end1 {
                    continue;
                }
                if prep.contains_n(q_pos, seed_len) {
                    continue;
                }

                for &t_pos_i32 in &t_sa[m.target_interval.start..m.target_interval.end] {
                    let t_pos = t_pos_i32 as usize;
                    if t_pos + seed_len > target.sequence.len() {
                        continue;
                    }
                    candidates.push(SeedCandidate {
                        query_pos: q_pos,
                        target_idx,
                        target_start: t_pos,
                        len: seed_len,
                        strand: Strand::Forward,
                    });
                }
            }
        }
    }

    // Reverse strand: use pre-built reverse_sa and sequence_rc
    let t_rc_sa = &target.reverse_sa;
    let t_rc = &target.sequence_rc;
    for seed_len in mi_len..=q_len {
        matches.clear(); // Reuse allocation
        let searcher = ParallelSaSearcher::new(&prep.q_rc_sa, &prep.q_rc, t_rc_sa, t_rc, config);
        searcher.find_seeds_into(seed_len, matches);

        for m in matches.iter() {
            for &q_rc_pos_i32 in &prep.q_rc_sa[m.query_interval.start..m.query_interval.end] {
                let q_rc_pos = q_rc_pos_i32 as usize;
                if q_rc_pos + seed_len > q_len {
                    continue;
                }
                let q_pos = q_len - q_rc_pos - seed_len;
                // Check: seed must start >= start0 AND end <= end0
                if q_pos < start0 || q_pos + seed_len > end1 {
                    continue;
                }
                if prep.contains_n(q_pos, seed_len) {
                    continue;
                }

                for &t_pos_i32 in &t_rc_sa[m.target_interval.start..m.target_interval.end] {
                    let t_pos = t_pos_i32 as usize;
                    if t_pos + seed_len > t_rc.len() {
                        continue;
                    }
                    candidates.push(SeedCandidate {
                        query_pos: q_pos,
                        target_idx,
                        target_start: t_pos,
                        len: seed_len,
                        strand: Strand::Reverse,
                    });
                }
            }
        }
    }
}

/// Find all seed matches between a query and a single target sequence.
/// Convenience wrapper that returns a new Vec.
pub fn find_seeds_in_target(
    prep: &QueryPrep,
    target: &SequenceIndex,
    target_idx: usize,
    config: &SeedArgs,
) -> Vec<SeedCandidate> {
    let mut candidates = Vec::new();
    let mut matches = Vec::new();
    find_seeds_in_target_into(prep, target, target_idx, config, &mut candidates, &mut matches);
    candidates
}

/// Find all seed matches between a query and all targets in an index.
///
/// Convenience wrapper that iterates over all targets. For parallelization,
/// use `find_seeds_in_target` directly with Rayon.
pub fn find_seeds(query: &[u8], index: &SaIndexFile, config: &SeedArgs) -> Vec<SeedCandidate> {
    let mut candidates = Vec::new();
    let mut matches = Vec::new();
    find_seeds_into(query, index, config, &mut candidates, &mut matches);
    candidates
}

/// Find all seed matches, reusing provided Vecs to avoid allocation.
///
/// Clears `candidates` and `matches` before filling.
pub fn find_seeds_into(
    query: &[u8],
    index: &SaIndexFile,
    config: &SeedArgs,
    candidates: &mut Vec<SeedCandidate>,
    matches: &mut Vec<ParallelSeedMatch>,
) {
    candidates.clear();
    matches.clear();

    // Pre-compute query data once (SA, RC, etc.)
    let Some(prep) = QueryPrep::new(query, config) else {
        return;
    };

    for (idx, target) in index.sequences.iter().enumerate() {
        find_seeds_in_target_into(&prep, target, idx, config, candidates, matches);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // helper: parse and normalize, returning the tuple or panic with message
    fn parse_and_norm(s: &str, qlen: usize) -> Result<(usize, usize, usize), String> {
        let spec = SeedSpec::from_str(s).map_err(|e| format!("parse err: {}", e))?;
        spec.normalize(qlen)
    }

    #[test]
    fn test_length_only() {
        let res = parse_and_norm("10", 100).expect("should parse");
        assert_eq!(res, (1, 100, 10));
    }

    #[test]
    fn test_interval_only() {
        let res = parse_and_norm("10:20", 100).expect("should parse");
        // interval 10..20 inclusive -> length 11
        assert_eq!(res, (10, 20, 11));
    }

    #[test]
    fn test_interval_with_length() {
        let res = parse_and_norm("10:20/5", 100).expect("should parse");
        assert_eq!(res, (10, 20, 5));
    }

    #[test]
    fn test_negative_interval() {
        let res = parse_and_norm("-5:-1", 100).expect("should parse");
        // -5 -> 96, -1 -> 100  (1-based), length = 5
        assert_eq!(res, (96, 100, 5));
    }

    #[test]
    fn test_invalid_length_zero() {
        let err = parse_and_norm("0", 100).unwrap_err();
        assert!(err.contains("Invalid seed length"));
    }

    #[test]
    fn test_mixed_sign_interval() {
        let err = parse_and_norm("5:-1", 100).unwrap_err();
        assert!(err.contains("Invalid seed interval"));
    }

    #[test]
    fn test_length_exceeds_interval() {
        let err = parse_and_norm("10:12/5", 100).unwrap_err();
        assert!(err.contains("exceeds interval"));
    }

    #[test]
    fn test_bad_parse() {
        assert!(SeedSpec::from_str("abc").is_err());
    }
}
