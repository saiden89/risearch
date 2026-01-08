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
    pub max_mismatches: u32,
    /// Minimum consecutive matches required at seed start/end
    pub min_consecutive: u32,
}

impl FromStr for MismatchSpec {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        if s.is_empty() {
            return Err("empty mismatch spec".into());
        }

        let parts: Vec<&str> = s.split(':').collect();
        if parts.len() != 2 {
            return Err(format!(
                "invalid mismatch spec '{}': expected 'c:p' format",
                s
            ));
        }

        let max_mismatches = parts[0]
            .parse::<u32>()
            .map_err(|e| format!("invalid max mismatches: {}", e))?;
        let min_consecutive = parts[1]
            .parse::<u32>()
            .map_err(|e| format!("invalid min consecutive: {}", e))?;

        Ok(MismatchSpec {
            max_mismatches,
            min_consecutive,
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
) -> Vec<Pairing> {
    let mut seed_alignment = Vec::with_capacity(len);
    for n in 0..len {
        let q_idx = q_pos + n;
        let t_idx = t_match_end.saturating_sub(n);
        let q_b = query.base(q_idx);
        let t_b = target.base(t_idx);
        seed_alignment.push(Pairing::from_bases(q_b, t_b));
    }
    seed_alignment
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
