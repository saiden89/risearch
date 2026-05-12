use clap::ValueEnum;

use crate::types::{DsmId, Energy};

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

impl OutputCompression {
    /// Return a suitable file extension (including leading dot) for this compression
    /// configuration, assuming TSV format, e.g. `".tsv"`, `".tsv.gz"`, `".tsv.zst"`.
    pub fn extension(&self) -> &'static str {
        match self {
            Self::None => ".tsv",
            Self::Gzip(_) => ".tsv.gz",
            Self::Zstd(_) => ".tsv.zst",
        }
    }
}


/// Arguments for seed generation.
#[derive(Debug, Clone, Default)]
pub struct SeedConfig {
    /// Seed interval start (1-based, can be negative). `None` = full query.
    pub seed_start: Option<i64>,
    /// Seed interval end (1-based, can be negative). `None` = full query.
    pub seed_end: Option<i64>,
    /// Seed length. In length-only mode this is the desired seed length.
    /// In interval mode, `None` means use the full interval width.
    pub seed_length: Option<i64>,

    /// Allow G-U wobble pairs when locating and maximizing seeds.
    pub seed_wobble: bool,
    
    /// Maximum number of mismatches allowed in seed.
    pub max_mismatches: usize,
    /// Minimum consecutive matches at seed start (prefix, 5').
    pub min_prefix_matches: usize,
    /// Minimum consecutive matches at seed end (suffix, 3').
    pub min_suffix_matches: usize,
}

/// Default seed length when none is specified.
pub const DEFAULT_SEED_LEN: i64 = 6;

impl SeedConfig {
    /// Resolve this config's seed bounds against a specific query length.
    pub fn resolve(&self, query_len: usize) -> Result<(usize, usize, usize), String> {
        let n = query_len as i64;
        if n <= 0 {
            return Err("query length must be positive".into());
        }

        match (self.seed_start, self.seed_end) {
            (None, None) => {
                let effective_length = self.seed_length.unwrap_or(DEFAULT_SEED_LEN);
                if effective_length <= 0 {
                    return Err("Invalid seed length".into());
                }
                let length = effective_length.min(n) as usize;
                Ok((1, query_len, length))
            }
            (Some(start), Some(end)) => {
                let (s_pos, e_pos) = match (start.signum(), end.signum()) {
                    (1, 1) if end >= start && end <= n => (start as usize, end as usize),
                    (-1, -1) if end >= start && start >= -n => {
                        let to_pos = |v: i64| -> usize {
                            let x = v + 1;
                            let r = (x + n) % n;
                            if r == 0 {
                                n as usize
                            } else {
                                r as usize
                            }
                        };
                        (to_pos(start), to_pos(end))
                    }
                    _ => return Err("Invalid seed interval: mixed sign".into()),
                };

                if s_pos == 0 || e_pos == 0 {
                    return Err("Invalid seed interval".into());
                }

                let interval_len = e_pos.saturating_sub(s_pos) + 1;
                if interval_len == 0 {
                    return Err("Invalid seed interval: empty".into());
                }

                let final_len = match self.seed_length {
                    None => interval_len,
                    // 0 or negative length explicitly requested -> use full interval
                    Some(l) if l <= 0 => interval_len,
                    Some(l) if (l as usize) > interval_len => {
                        return Err("Invalid seed length (exceeds interval)".into())
                    }
                    Some(l) => l as usize,
                };

                Ok((s_pos, e_pos, final_len))
            }
            _ => Err("use seed_start + seed_end together, or neither".into()),
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
