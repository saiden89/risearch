use clap::ValueEnum;

pub use crate::types::{SeedPairing, Strand};

#[derive(ValueEnum, Clone, Debug)]
#[clap(rename_all = "snake_case")]
pub enum Matrix {
    T99,
    T04,
}

#[derive(clap::ValueEnum, Clone, Debug)]
#[value(rename_all = "lowercase")]
pub enum OutputFormat {
    /// Report predictions in detailed format
    Detailed,
    /// Report predictions in a simple format together with CIGAR-like string for interaction structure
    Cigar,
    /// Report predictions in a simple format together with binding site (3'->5'), flanking 5'end (3'->5') and flanking 3'end (5'->3') sequences of the target x(required for post-processing of CRISPR off-target predictions)
    BindingSite,
}

#[derive(ValueEnum, Clone, Debug)]
#[clap(rename_all = "lowercase")]
pub enum Backend {
    /// FM-Index backend (requires 'fm-index' feature)
    Fm,
    /// Suffix Array backend
    Sa,
}

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
    pub seed: Option<String>, // TODO: this should be a SeedSpec

    /// Consider G-U wobble pairs as mismatch within the seed (only for locating seeds, energy model is not affected)
    #[arg(
        short = 'w',
        long = "wobble",
        value_enum,
        default_value_t = SeedPairing::AllowWobble,
        default_missing_value = "strict",
        num_args = 0
    )]
    pub pairing: SeedPairing,

    /// Introduce mismatched seeds
    /// Set the max num of mismatches (c) allowed in the seed and min num of consecutive matches required at seed start/end (p)
    ///These seeds will not overlap with perfect complementary seeds.
    #[arg(
        short = 'm',
        long = "mismatch",
        value_name = "c:p",
        default_value = "0:0"
    )]
    pub mismatch_seed: Option<String>, // TODO: this should be a MismatchSpec
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
        require_equals = true,
        value_enum
    )]
    pub report_format: Option<OutputFormat>,
}
