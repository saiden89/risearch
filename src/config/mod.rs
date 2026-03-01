use clap::ValueEnum;

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

                // Determine final length
                let length_opt = match *self {
                    SeedSpec::IntervalWithLength { length, .. } => Some(length),
                    SeedSpec::Interval { .. } => None,
                    SeedSpec::LengthOnly(_) => None,
                };

                let final_len = match length_opt {
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
pub enum Matrix {
    /// Turner 1999 RNA-RNA parameters
    T99,
    /// Turner 2004 RNA-RNA parameters (default)
    #[default]
    T04,
}

impl std::fmt::Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Matrix::T99 => write!(f, "t99"),
            Matrix::T04 => write!(f, "t04"),
        }
    }
}

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

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq, Default)]
#[clap(rename_all = "lowercase")]
pub enum OutputCompression {
    #[default]
    None,
    #[value(alias = "gz")]
    Gzip,
    #[value(alias = "zst")]
    Zstd,
}

// =============================================================================
// CONFIG TYPES
// =============================================================================

use crate::types::SeedPairingMode;

/// Arguments for seed generation
#[derive(Debug, Clone)]
pub struct SeedConfig {
    /// Seed specification (length or interval)
    pub seed: SeedSpec,

    // TODO: remove after legacy -U flag is dropped — redundant with `pairing`
    /// DEPRECATED (will be removed in a future release): disable G-U wobble pairs within the seed
    pub no_guseed: bool,

    /// Seed pairing mode (allow_wobble or strict)
    pub pairing: SeedPairingMode,

    // TODO: collapse mismatch + mismatch_{max,prefix,suffix} into a single
    // MismatchSpec after the legacy -m flag is removed
    /// Mismatch specification (legacy -m or named overrides)
    pub mismatch: MismatchSpec,

    /// Max number of mismatches allowed in the seed (preferred)
    pub mismatch_max: Option<usize>,

    /// Min consecutive matches at seed start (prefix / 5')
    pub mismatch_prefix: Option<usize>,

    /// Min consecutive matches at seed end (suffix / 3')
    pub mismatch_suffix: Option<usize>,
}

impl SeedConfig {
    pub fn with_wobble(seed: SeedSpec, mismatch: MismatchSpec, allow_wobble: bool) -> Self {
        let pairing = if allow_wobble {
            SeedPairingMode::AllowWobble
        } else {
            SeedPairingMode::Strict
        };
        Self {
            seed,
            no_guseed: !allow_wobble,
            pairing,
            mismatch,
            mismatch_max: None,
            mismatch_prefix: None,
            mismatch_suffix: None,
        }
    }

    #[inline]
    pub fn allows_wobble(&self) -> bool {
        matches!(self.pairing, SeedPairingMode::AllowWobble)
    }

    pub fn apply_pairing_overrides(&mut self, explicit_pairing: bool) {
        if self.no_guseed && !explicit_pairing {
            self.pairing = SeedPairingMode::Strict;
        }
    }

    pub fn apply_mismatch_overrides(&mut self) {
        if let Some(max) = self.mismatch_max {
            self.mismatch.max_mismatches = max;
        }
        if let Some(start) = self.mismatch_prefix {
            self.mismatch.min_prefix_matches = start;
        }
        if let Some(end) = self.mismatch_suffix {
            self.mismatch.min_suffix_matches = end;
        }
    }

    pub fn has_named_mismatch(&self) -> bool {
        self.mismatch_max.is_some()
            || self.mismatch_prefix.is_some()
            || self.mismatch_suffix.is_some()
    }
}

/// Global scoring model shared by seed scoring and DP extension.
#[derive(Debug, Clone)]
pub struct ScoreConfig {
    /// Energy matrix for RNA-RNA duplexes
    pub matrix: Matrix,

    /// Per-nucleotide penalty (in kcal/mol), applied wherever score is computed
    pub penalty: f64,

    // Placeholder flags from C implementation - not yet implemented
    // (see cli/args for detailed documentation)
    pub matrix2: Option<String>,
    pub matpath: Option<String>,
    pub temperature: Option<String>,
    pub weights: Option<String>,
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
    pub delta_g: f64,

    /// Energy per length threshold that filters seeds
    pub seed_energy: f64,

    /// Disable maximality check (allows redundant seeds)
    pub no_max_prune: bool,

    /// Disable shadow dedup filtering (keep hits contained by better hits)
    /// Default: true (filters contained hits). Set --no-dedup-shadow for C-compatible behavior.
    pub dedup_shadow: bool,
}

/// Options that apply to the `search` subcommand
#[derive(Debug, Clone)]
pub struct SearchArgs {
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

    /// Output compression codec (overrides file extension inference; gzip/gz, zstd/zst accepted)
    pub compress: Option<OutputCompression>,

    /// Output compression level (codec-specific)
    pub level: Option<i32>,

    /// Write one output file per query (directory mode)
    pub multifile: bool,
}

#[cfg(test)]
mod spec_tests;
