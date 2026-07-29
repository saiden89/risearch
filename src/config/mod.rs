use clap::ValueEnum;

use crate::dp::MAX_EXT;
use crate::dsm::DsmRegistry;
use crate::types::{DsmId, Energy};

// `Default` implementations below own frontend default policy. Keep constants
// here only for effective fallbacks and shared validity bounds or sentinels.
pub const DEFAULT_SEED_LEN: i64 = 6;

pub const MIN_PENALTY_KCAL: f64 = 0.0;
pub const MAX_PENALTY_KCAL: f64 = 50.0;
pub const MIN_TEMPERATURE_C: i32 = 0;
pub const MAX_TEMPERATURE_C: i32 = 100;

pub const UNLIMITED_EXTENSION: i32 = -1;
pub const MAX_EXTENSION: i32 = MAX_EXT as i32;

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

impl OutputFormat {
    /// Whether this format prints pairing data, and so needs the extension to
    /// run traceback. The single source of truth for
    /// [`ExtendConfig::build_alignment`].
    pub const fn needs_alignment(self) -> bool {
        !matches!(self, Self::Minimal)
    }
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
#[derive(Debug, Clone)]
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

    /// Emit non-maximal seeds too — those that could be grown by one more
    /// pairing column, and so are shorter copies of a longer match.
    pub no_max_prune: bool,

    /// Maximum number of mismatches allowed in seed.
    pub max_mismatches: usize,
    /// Minimum consecutive matches at seed start (prefix, 5').
    pub min_prefix_matches: usize,
    /// Minimum consecutive matches at seed end (suffix, 3').
    pub min_suffix_matches: usize,
}

impl Default for SeedConfig {
    fn default() -> Self {
        Self {
            seed_start: None,
            seed_end: None,
            seed_length: None,
            seed_wobble: true,
            no_max_prune: false,
            max_mismatches: 0,
            min_prefix_matches: 1,
            min_suffix_matches: 0,
        }
    }
}

impl SeedConfig {
    /// Validate query-independent seed invariants.
    ///
    /// Bounds against a particular sequence length, and a requested length
    /// against that normalized interval, are checked by [`Self::resolve`].
    pub fn validate(&self) -> Result<(), String> {
        match (self.seed_start, self.seed_end) {
            (None, None) => {
                if self.seed_length.is_some_and(|length| length <= 0) {
                    return Err("Invalid seed length: must be positive without an interval".into());
                }
            }
            (Some(start), Some(end)) => {
                if start == 0 || end == 0 {
                    return Err(
                        "Invalid seed interval: coordinates are 1-based and cannot be zero".into(),
                    );
                }
                if start.signum() != end.signum() {
                    return Err("Invalid seed interval: mixed sign".into());
                }
                if end < start {
                    return Err("Invalid seed interval: end precedes start".into());
                }

                // In interval mode, zero or a negative length means the full
                // interval, preserving the legacy CLI's established behavior.
                if let Some(length) = self.seed_length.filter(|length| *length > 0) {
                    let interval_len = i128::from(end) - i128::from(start) + 1;
                    if i128::from(length) > interval_len {
                        return Err("Invalid seed length (exceeds interval)".into());
                    }
                }
            }
            _ => return Err("use seed_start + seed_end together, or neither".into()),
        }

        Ok(())
    }

    /// Resolve this config's seed bounds against a specific query length.
    pub fn resolve(&self, query_len: usize) -> Result<(usize, usize, usize), String> {
        self.validate()?;

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

impl Default for ScoreConfig {
    fn default() -> Self {
        Self {
            dsm_id: DsmId::from("t04"),
            penalty: Energy::default(),
            temperature: 37,
        }
    }
}

impl ScoreConfig {
    pub fn validate(&self) -> Result<(), String> {
        DsmRegistry::parse_id(self.dsm_id.0.as_str()).map_err(|err| err.to_string())?;

        let penalty = self.penalty.to_kcal();
        if !(MIN_PENALTY_KCAL..=MAX_PENALTY_KCAL).contains(&penalty) {
            return Err(format!(
                "penalty must be between {MIN_PENALTY_KCAL} and {MAX_PENALTY_KCAL}, got {penalty}"
            ));
        }
        if !(MIN_TEMPERATURE_C..=MAX_TEMPERATURE_C).contains(&self.temperature) {
            return Err(format!(
                "temperature must be between {MIN_TEMPERATURE_C} and {MAX_TEMPERATURE_C}, got {}",
                self.temperature
            ));
        }

        Ok(())
    }
}

/// Arguments for extension strategy.
#[derive(Debug, Clone)]
pub struct ExtendConfig {
    /// Max extension length on the seed (do DP for max this length up- and downstream of seed).
    /// Negative = unlimited: extend across the whole query. Queries longer than
    /// MAX_EXT per side are rejected rather than silently clamped.
    pub max_extension: i32,

    /// Record the pairing alignment during DP traceback, populating
    /// [`SearchHit::alignment`](crate::SearchHit). Off skips that work entirely.
    ///
    /// Dedup's tie-break compares fingerprints, so with this off a different
    /// member of an energy-tied bounding box may survive. Tied hits share a box
    /// and an energy, so every other field is identical either way.
    pub build_alignment: bool,
}

impl Default for ExtendConfig {
    fn default() -> Self {
        Self {
            max_extension: 20,
            build_alignment: true,
        }
    }
}

impl ExtendConfig {
    pub fn validate(&self) -> Result<(), String> {
        if !(UNLIMITED_EXTENSION..=MAX_EXTENSION).contains(&self.max_extension) {
            return Err(format!(
                "max extension must be {UNLIMITED_EXTENSION} (unlimited) or between 0 and {MAX_EXTENSION}, got {}",
                self.max_extension
            ));
        }
        Ok(())
    }
}

/// Hit acceptance and pruning policies.
#[derive(Debug, Clone)]
pub struct FilterConfig {
    /// Set deltaG energy threshold (in kcal/mol) to filter predictions
    pub delta_g: Energy,

    /// Energy per length threshold that filters seeds
    pub seed_energy: Energy,

    /// Report every maximal seed as its own hit. When unset (the default),
    /// hits whose extension resolves to the same final bounding box are
    /// collapsed to the single lowest-energy alignment.
    pub no_dedup: bool,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            delta_g: Energy::from_kcal(-20.0),
            seed_energy: Energy::default(),
            no_dedup: false,
        }
    }
}

/// Options that apply to the `search` subcommand
#[derive(Debug, Clone, Default)]
pub struct SearchConfig {
    pub seed: SeedConfig,
    pub score: ScoreConfig,
    pub extend: ExtendConfig,
    pub filter: FilterConfig,
}

impl SearchConfig {
    pub fn validate(&self) -> Result<(), String> {
        self.seed.validate()?;
        self.score.validate()?;
        self.extend.validate()?;
        Ok(())
    }
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
