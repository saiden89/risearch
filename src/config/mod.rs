use clap::ValueEnum;

use crate::types::{DsmId, Energy};

// =============================================================================
// SEED SPECIFICATION TYPES (moved from seed::spec)
// =============================================================================

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

/// Representation of the seed specification:
/// - length only         => SeedSpec::LengthOnly(l)
/// - interval            => SeedSpec::Interval { start, end, length: None }
/// - interval + length   => SeedSpec::Interval { start, end, length: Some(l) }
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeedSpec {
    LengthOnly(i64),
    Interval {
        start: i64, // TODO: consider using RangeInclusive<i64>
        end: i64,
        length: Option<i64>, // None = full interval; Some = constrained seed length
    },
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
                            if r == 0 {
                                n as usize
                            } else {
                                r as usize
                            }
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

                let final_len = match length {
                    Some(length) if length <= 0 => {
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

// =============================================================================
// ENUMS (shared by config and CLI via clap derives)
// =============================================================================


#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq, Default)]
#[clap(rename_all = "lowercase")]
pub enum OutputFormat {
    /// Detailed format with alignment (C: -p1)
    Detailed,
    /// CIGAR-like interaction structure (C: -p2)
    Cigar,
    /// Binding site with flanks for CRISPR (C: -p3)
    BindingSite,
    /// Minimal: target, position, strand, energy (C: -p4)
    #[default]
    Minimal,
}

/// CLI-facing codec selector — what the `--compress` flag parses into.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
#[clap(rename_all = "lowercase")]
pub enum OutputCodec {
    None,
    #[value(alias = "gz")]
    Gzip,
    #[value(alias = "zst")]
    Zstd,
}

impl From<&std::path::Path> for OutputCodec {
    /// Infer codec from file extension. Unrecognised or absent → `None`.
    fn from(path: &std::path::Path) -> Self {
        match path.extension().and_then(|e| e.to_str()) {
            Some(ext) => match ext.to_ascii_lowercase().as_str() {
                "gz" | "gzip" => Self::Gzip,
                "zst" | "zstd" => Self::Zstd,
                _ => Self::None,
            },
            None => Self::None,
        }
    }
}

/// Config-facing compression value — codec and level bound together.
/// Illegal states (level without codec) are unrepresentable.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum OutputCompression {
    #[default]
    None,
    Gzip(u8),  // level 0–9
    Zstd(i32), // level -7..22
}

// =============================================================================
// CONFIG TYPES
// =============================================================================

/// Arguments for seed generation
#[derive(Debug, Clone)]
pub struct SeedConfig {
    /// Seed specification (length or interval)
    pub seed: SeedSpec,

    /// Allow G-U wobble pairs when locating and maximizing seeds.
    pub seed_wobble: bool,

    /// Mismatch specification used during seed search.
    pub mismatch: MismatchSpec,
}

impl SeedConfig {
    pub fn with_wobble(seed: SeedSpec, mismatch: MismatchSpec, seed_wobble: bool) -> Self {
        Self {
            seed,
            seed_wobble,
            mismatch,
        }
    }
}

/// Global scoring model shared by seed scoring and DP extension.
#[derive(Debug, Clone)]
pub struct ScoreConfig {
    pub dsm_id: DsmId,
    pub penalty: Energy,
    pub temperature: i32,
}

/// Arguments for extension strategy.
#[derive(Debug, Clone)]
pub struct ExtendConfig {
    /// Max extension length on the seed (do DP for max this length up- and downstream of seed)
    pub max_extension: u8,

    // Placeholder flags from C implementation - not yet implemented
    // (see cli/args for detailed documentation)
    pub band: Option<u32>,
}

/// Hit acceptance and pruning policies.
#[derive(Debug, Clone)]
pub struct FilterConfig {
    /// Set deltaG energy threshold (in kcal/mol) to filter predictions
    pub delta_g: Energy,

    /// Energy per length threshold that filters seeds
    pub seed_energy: Energy,

    /// Disable maximality check (allows redundant seeds)
    pub no_max_prune: bool,
}

/// Options that apply to the `search` subcommand
#[derive(Debug, Clone)]
pub struct SearchConfig {
    pub seed: SeedConfig,
    pub score: ScoreConfig,
    pub extend: ExtendConfig,
    pub filter: FilterConfig,
    pub output: OutputConfig,

    // Placeholder flags from C implementation - not yet implemented
    // (see cli/args for detailed documentation)
    pub one_vs_one: bool,
    pub three_prime_match: Option<String>,
    pub five_prime_match: Option<String>,
}

#[derive(Debug, Clone)]
pub struct OutputConfig {
    /// Output format
    pub format: OutputFormat,

    /// Compression codec and level, bound together (resolved at arg boundary)
    pub compress: OutputCompression,

    /// Write one output file per query (directory mode)
    pub multifile: bool,
}

#[cfg(test)]
mod spec_tests;
