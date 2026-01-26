use std::str::FromStr;

/// Representation of the `-m` (mismatch) flag.
///
/// Format: `c[:ps[:pe]]` where:
/// - `c` = max number of mismatches allowed in seed
/// - `ps` = min number of consecutive matches required at seed start
/// - `pe` = min number of consecutive matches required at seed end
///
/// Examples:
/// - `-m 1`     allows 1 mismatch with 1-match protected ends (C-compatible)
/// - `-m 1:3`   allows 1 mismatch with 3 matches at both ends
/// - `-m 1:3:5` allows 1 mismatch with 3 matches at start and 5 at end
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MismatchSpec {
    /// Maximum number of mismatches allowed in seed
    pub max_mismatches: usize,
    /// Minimum consecutive matches at seed start (prefix, 5')
    pub min_prefix_matches: usize,
    /// Minimum consecutive matches at seed end (suffix, 3')
    pub min_suffix_matches: usize,
}

impl MismatchSpec {
    // TODO: Consider validating against seed length (e.g., prefix/suffix > seed_len)
    // to avoid configurations that yield zero hits, while preserving C compatibility.
    /// No mismatches allowed (exact matching)
    pub const fn exact() -> Self {
        Self {
            max_mismatches: 0,
            min_prefix_matches: 1,
            min_suffix_matches: 0,
        }
    }

    /// Create with specific parameters
    pub const fn new(max: usize, min_prefix: usize, min_suffix: usize) -> Self {
        Self {
            max_mismatches: max,
            min_prefix_matches: min_prefix,
            min_suffix_matches: min_suffix,
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

        let (max_mismatches, min_start, min_end) = match parts.len() {
            1 => {
                // Allow "c" as shorthand for "c:c:c" (matches C behavior)
                let max = parts[0]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid max mismatches: {}", e))?;
                (max, max, max)
            }
            2 => {
                let max = parts[0]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid max mismatches: {}", e))?;
                let min = parts[1]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid min consecutive: {}", e))?;
                (max, min, min)
            }
            3 => {
                let max = parts[0]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid max mismatches: {}", e))?;
                let min_start = parts[1]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid min start matches: {}", e))?;
                let min_end = parts[2]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid min end matches: {}", e))?;
                (max, min_start, min_end)
            }
            _ => {
                return Err(format!(
                    "invalid mismatch spec '{}': expected 'c', 'c:p', or 'c:ps:pe' format",
                    s
                ));
            }
        };

        // CLI format c:p / c:ps:pe maps to: max=c, min_prefix_matches=ps, min_suffix_matches=pe
        Ok(MismatchSpec {
            max_mismatches,
            min_prefix_matches: min_start,
            min_suffix_matches: min_end,
        })
    }
}

/// Representation of the `-s` flag:
/// - `-s l`              => SeedSpec::LengthOnly(l)
/// - `-s m:n`            => SeedSpec::Interval { start: m, end: n }
/// - `-s m:n/l`          => SeedSpec::IntervalWithLength { start: m, end: n, length: l }
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeedSpec {
    LengthOnly(i64),
    Interval {
        start: i64, // TODO: consider using RangeInclusive<i64>
        end: i64,
    },
    IntervalWithLength {
        start: i64, // TODO: consider using RangeInclusive<i64>
        end: i64,
        length: i64, // TODO: enforce strictly positive via newtype
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
                Ok(SeedSpec::LengthOnly(len))
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
                        Ok(SeedSpec::Interval { start, end })
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

                        Ok(SeedSpec::IntervalWithLength {
                            start,
                            end,
                            length,
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
            SeedSpec::LengthOnly(l) if l <= 0 => Err("Invalid seed length".into()),
            SeedSpec::LengthOnly(l) => {
                let length = l.min(n) as usize;
                Ok((1, query_len, length))
            }
            SeedSpec::Interval { start, end } | SeedSpec::IntervalWithLength { start, end, .. } => {
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
                let length_opt = match *self {
                    SeedSpec::IntervalWithLength { length, .. } => Some(length),
                    SeedSpec::Interval { .. } => None,
                    SeedSpec::LengthOnly(_) => None,
                };

                let final_len = match length_opt {
                    Some(length) if length < 0 => {
                        return Err("Invalid seed length".into());
                    }
                    Some(length) if (length as usize) > interval_len => {
                        return Err("Invalid seed length (exceeds interval)".into());
                    }
                    Some(length) => length as usize,
                    None => interval_len,
                };

                Ok((s_pos, e_pos, final_len))
            }
        }
    }
}
