use clap::ValueEnum;

pub use crate::types::{SeedPairing, Strand};

#[derive(ValueEnum, Clone, Debug)]
#[clap(rename_all = "snake_case")]
pub enum Matrix {
    /// Turner 1999 RNA-RNA parameters
    T99,
    /// Turner 2004 RNA-RNA parameters (default)
    T04,

    // ========================================================================
    // TODO: Placeholder matrix types from C implementation - not yet implemented
    // ========================================================================
    /// TODO: SantaLucia 1995 RNA-DNA duplex parameters
    #[value(name = "su95")]
    Su95,

    /// TODO: SantaLucia 1995 RNA-DNA modified for CRISPRoff2
    #[value(name = "su95c2")]
    Su95c2,

    /// TODO: SantaLucia 1995 RNA-DNA with mismatches as loop size 2
    #[value(name = "su95wk11")]
    Su95wk11,

    /// TODO: SantaLucia 1995 RNA-DNA without G-U wobble pairs
    #[value(name = "su95_nogu")]
    Su95NoGU,

    /// TODO: SantaLucia 2004 DNA-DNA without G-T wobble pairs
    #[value(name = "sl04_nogu")]
    Sl04NoGU,
}

#[derive(clap::ValueEnum, Clone, Copy, Debug)]
#[value(rename_all = "lowercase")]
pub enum OutputFormat {
    /// Report predictions in detailed format (C: -p or -p1)
    Detailed,
    /// Report predictions in a simple format together with CIGAR-like string for interaction structure (C: -p2)
    Cigar,
    /// Report predictions in a simple format together with binding site (3'->5'), flanking 5'end (3'->5') and flanking 3'end (5'->3') sequences of the target (required for post-processing of CRISPR off-target predictions) (C: -p3)
    BindingSite,

    // ========================================================================
    // TODO: Placeholder output format from C implementation - not yet implemented
    // ========================================================================
    /// TODO: Minimal format - target, start, strand, and energy only (C: -p4)
    Minimal,
}

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
#[clap(rename_all = "lowercase")]
pub enum OutputCompression {
    /// No compression (default)
    None,
    /// Gzip compression
    Gzip,
    /// Zstandard compression
    Zstd,
}

#[derive(ValueEnum, Clone, Debug)]
#[clap(rename_all = "lowercase")]
pub enum Backend {
    /// FM-Index backend (requires 'fm-index' feature)
    Fm,
    /// Suffix Array backend
    Sa,
}

use crate::seed::{MismatchSpec, SeedSpec};

/// Arguments for seed generation
#[derive(clap::Args, Debug, Clone)]
pub struct SeedArgs {
    /// Set seed length (-s l = length only; -s n:m = full interval)
    #[arg(
        short = 's',
        long = "seed",
        value_name = "start:end/length",
        default_value = "6"
    )]
    pub seed: SeedSpec,

    /// DEPRECATED (will be removed in a future release): disable G-U wobble pairs within the seed
    #[arg(
        short = 'U',
        long = "no-guseed",
        alias = "noGUseed",
        action = clap::ArgAction::SetTrue,
        help_heading = "Deprecated"
    )]
    pub no_guseed: bool,

    /// Allow G-U wobble pairs within the seed (default behavior)
    #[arg(short = 'w', long = "wobble", action = clap::ArgAction::SetTrue)]
    pub wobble_legacy: bool,

    /// Seed pairing mode (allow_wobble or strict)
    #[arg(
        long = "seed-pairing",
        value_enum,
        default_value_t = SeedPairing::AllowWobble
    )]
    pub pairing: SeedPairing,

    /// DEPRECATED (will be removed in a future release): legacy mismatch shorthand
    /// Set max mismatches (c) and min consecutive matches at seed start/end (p)
    /// These seeds will not overlap with perfect complementary seeds.
    #[arg(
        short = 'm',
        long = "mismatch",
        value_name = "c[:ps[:pe]]",
        default_value = "0:0",
        help_heading = "Deprecated"
    )]
    pub mismatch_seed: MismatchSpec,

    /// Max number of mismatches allowed in the seed (preferred)
    #[arg(long = "mismatch-max", value_name = "C")]
    pub mismatch_max: Option<usize>,

    /// Min consecutive matches at seed start (prefix / 5')
    #[arg(long = "mismatch-prefix", value_name = "PS")]
    pub mismatch_prefix: Option<usize>,

    /// Min consecutive matches at seed end (suffix / 3')
    #[arg(long = "mismatch-suffix", value_name = "PE")]
    pub mismatch_suffix: Option<usize>,
}

impl SeedArgs {
    pub fn apply_pairing_overrides(&mut self, explicit_pairing: bool) {
        if self.no_guseed {
            self.pairing = SeedPairing::Strict;
            return;
        }
        if self.wobble_legacy && !explicit_pairing {
            self.pairing = SeedPairing::AllowWobble;
        }
    }

    pub fn apply_mismatch_overrides(&mut self) {
        if let Some(max) = self.mismatch_max {
            self.mismatch_seed.max_mismatches = max;
        }
        if let Some(start) = self.mismatch_prefix {
            self.mismatch_seed.min_position = start;
        }
        if let Some(end) = self.mismatch_suffix {
            self.mismatch_seed.min_matches_after = end;
        }
    }

    pub fn has_named_mismatch(&self) -> bool {
        self.mismatch_max.is_some() || self.mismatch_prefix.is_some() || self.mismatch_suffix.is_some()
    }
}

/// Arguments for seed extension and scoring
#[derive(clap::Args, Debug, Clone)]
pub struct ExtendArgs {
    /// Max extension length on the seed (do DP for max this length up- and downstream of seed)
    #[arg(
        short = 'l',
        long = "extension",
        value_name = "LENGTH",
        default_value_t = 20
    )]
    pub max_extension: u8, // TODO: should be strictly positive

    /// Set deltaG energy threshold (in kcal/mol) to filter predictions
    #[arg(
            short = 'e',
            long = "energy",
            value_name = "dG",
            default_value_t = -20.0,
            allow_hyphen_values = true
        )]
    pub delta_g: f64,

    /// Energy matrix for RNA-RNA duplexes
    #[arg(short = 'z', long = "matrix", value_name = "MATRIX", default_value_t = Matrix::T04, value_enum)]
    pub matrix: Matrix,

    /// Per-nucleotide extension penalty (in kcal/mol)
    #[arg(
        short = 'p',
        long = "penalty",
        value_name = "PENALTY",
        default_value_t = 0.0
    )]
    pub penalty: f64,

    /// Energy per length threshold that filters seeds
    #[arg(long = "seed-energy", value_name = "THRESHOLD", default_value_t = 0.0)]
    pub seed_energy: f64,

    /// Disable maximality check (allows redundant seeds)
    #[arg(long = "no-max-prune", action = clap::ArgAction::SetTrue)]
    pub no_max_prune: bool,

    /// Disable shadow dedup filtering (keep hits contained by better hits)
    /// Default: true (filters contained hits). Set --no-dedup-shadow for C-compatible behavior.
    #[arg(long = "no-dedup-shadow", action = clap::ArgAction::SetFalse, default_value_t = true)]
    pub dedup_shadow: bool,

    // ========================================================================
    // TODO: Placeholder flags from C implementation - not yet implemented
    // ========================================================================
    /// TODO: Banded search - limits the search for bulged matches.
    /// In C: `-b band, --band=band` - Integer size of bands limiting bulge search.
    /// The minimum size is 1; use seed option to avoid any bulge.
    #[arg(long = "band", value_name = "BAND", hide = true)]
    pub band: Option<u32>,

    /// TODO: Secondary energy matrix for custom energy parameters.
    /// In C: `-y mat2, --matrix2=mat2` - Only needed for custom energy matrices.
    #[arg(long = "matrix2", value_name = "MATRIX2", hide = true)]
    pub matrix2: Option<String>,

    /// TODO: Path to directory holding custom energy matrices.
    /// In C: `-M PATH, --matpath=PATH` - Directory with energy matrix files.
    #[arg(long = "matpath", value_name = "PATH", hide = true)]
    pub matpath: Option<String>,

    /// TODO: Temperature scaling for energy calculations.
    /// In C: `-K T1[,T2,T3], --temperature=T0[,T1,T2]` - Temperatures in Kelvin.
    /// T0 is the target temperature; T1/T2 only needed for custom energy parameters.
    #[arg(long = "temperature", value_name = "T1[,T2,T3]", hide = true)]
    pub temperature: Option<String>,

    /// TODO: CRISPR weighting for gRNA-target interactions.
    /// In C: `-w arr, --weights=arr` - Use "CRISPR_gRNApPAM" to weight by CRISPR/Cas9 impact.
    #[arg(long = "weights", value_name = "WEIGHTS", hide = true)]
    pub weights: Option<String>,
}

/// Options that apply to the `search` subcommand
#[derive(clap::Args, Debug, Clone)]
pub struct SearchArgs {
    #[command(flatten)]
    pub seed: SeedArgs,

    #[command(flatten)]
    pub extend: ExtendArgs,

    /// Output format
    #[arg(
        short = 'f',
        long = "format",
        value_name = "LEVEL",
        default_missing_value = "detailed",
        default_value = "detailed",
        value_enum
    )]
    pub report_format: Option<OutputFormat>,

    /// Output compression codec (overrides file extension inference)
    #[arg(long = "output-compress", value_enum)]
    pub output_compress: Option<OutputCompression>,

    /// Output compression level (codec-specific)
    #[arg(long = "output-level", value_name = "LEVEL")]
    pub output_level: Option<i32>,

    // ========================================================================
    // TODO: Placeholder flags from C implementation - not yet implemented
    // ========================================================================
    /// TODO: One-vs-one mode - only print results where query name matches target name.
    /// In C: `-1, --one_vs_one` - Filters results to matching query/target names.
    #[arg(short = '1', long = "one-vs-one", action = clap::ArgAction::SetTrue, hide = true)]
    pub one_vs_one: bool,

    /// TODO: 3' PAM filter - report only predictions matching a 3' PAM pattern.
    /// In C: `-3 <reg>, --three_prime_match=PC` - Regex for 3' PAM forward complement.
    /// Example for cas9: `^(.cc|.uc|.cu)` matching NGG/NAG/NGA 3' PAMs.
    /// Requires output format 3 or 4 (binding_site).
    #[arg(
        short = '3',
        long = "three-prime-match",
        value_name = "REGEX",
        hide = true
    )]
    pub three_prime_match: Option<String>,

    /// TODO: 5' PAM filter - report only predictions matching a 5' PAM pattern.
    /// In C: `-5 <reg>, --five_prime_match=PC` - Regex for 5' PAM reverse complement.
    /// Example for cas12a: `^([^a]aaa)` matching TTTV 5' PAMs.
    /// Requires output format 3 or 4 (binding_site).
    #[arg(
        short = '5',
        long = "five-prime-match",
        value_name = "REGEX",
        hide = true
    )]
    pub five_prime_match: Option<String>,
}
