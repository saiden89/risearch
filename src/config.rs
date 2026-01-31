use clap::ValueEnum;

use crate::seed::{MismatchSpec, SeedSpec};
use crate::types::SeedPairingMode;

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

/// Arguments for seed generation
#[derive(Debug, Clone)]
pub struct SeedConfig {
    /// Seed specification (length or interval)
    pub seed: SeedSpec,

    /// DEPRECATED (will be removed in a future release): disable G-U wobble pairs within the seed
    pub no_guseed: bool,

    /// Seed pairing mode (allow_wobble or strict)
    pub pairing: SeedPairingMode,

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

/// Arguments for seed extension and scoring
#[derive(Debug, Clone)]
pub struct ExtendConfig {
    /// Max extension length on the seed (do DP for max this length up- and downstream of seed)
    pub max_extension: u8,

    /// Set deltaG energy threshold (in kcal/mol) to filter predictions
    pub delta_g: f64,

    /// Energy matrix for RNA-RNA duplexes
    pub matrix: Matrix,

    /// Per-nucleotide extension penalty (in kcal/mol)
    pub penalty: f64,

    /// Energy per length threshold that filters seeds
    pub seed_energy: f64,

    /// Disable maximality check (allows redundant seeds)
    pub no_max_prune: bool,

    /// Disable shadow dedup filtering (keep hits contained by better hits)
    /// Default: true (filters contained hits). Set --no-dedup-shadow for C-compatible behavior.
    pub dedup_shadow: bool,

    // Placeholder flags from C implementation - not yet implemented
    // (see cli_args.rs for detailed documentation)
    pub band: Option<u32>,
    pub matrix2: Option<String>,
    pub matpath: Option<String>,
    pub temperature: Option<String>,
    pub weights: Option<String>,
}

/// Options that apply to the `search` subcommand
#[derive(Debug, Clone)]
pub struct SearchArgs {
    pub seed: SeedConfig,
    pub extend: ExtendConfig,
    pub output: OutputConfig,

    // Placeholder flags from C implementation - not yet implemented
    // (see cli_args.rs for detailed documentation)
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
}
