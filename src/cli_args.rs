use log::warn;

use crate::config::{self, Matrix, OutputCompression, OutputFormat};
use crate::seed::{MismatchSpec, SeedSpec};
use crate::types::SeedPairingMode;

/// Arguments for seed generation
#[derive(clap::Args, Debug, Clone)]
pub struct SeedConfig {
    /// DEPRECATED (will be removed in a future release): legacy seed spec
    /// Formats: "l", "m:n", "m:n/l"
    #[arg(
        short = 's',
        long = "seed",
        value_name = "start:end/length",
        default_value = "6",
        help_heading = "Deprecated"
    )]
    pub seed_legacy: SeedSpec,

    /// Seed interval start (1-based, can be negative)
    /// TODO: Consider explicit one-sided bounds (e.g. --seed-to-end/--seed-from-start)
    /// instead of inferring missing start/end. Keep strict parsing for now.
    #[arg(
        long = "seed-start",
        value_name = "START",
        allow_hyphen_values = true,
        requires = "seed_end"
    )]
    pub seed_start: Option<i64>,

    /// Seed interval end (1-based, can be negative)
    #[arg(
        long = "seed-end",
        value_name = "END",
        allow_hyphen_values = true,
        requires = "seed_start"
    )]
    pub seed_end: Option<i64>,

    /// Seed length (use alone, or with seed-start/seed-end to constrain interval)
    #[arg(long = "seed-length", value_name = "LENGTH")]
    pub seed_length: Option<i64>,

    /// DEPRECATED (will be removed in a future release): disable G-U wobble pairs within the seed
    #[arg(
        short = 'U',
        long = "no-guseed",
        alias = "noGUseed",
        action = clap::ArgAction::SetTrue,
        help_heading = "Deprecated"
    )]
    pub no_guseed: bool,

    /// Seed pairing mode (allow_wobble or strict)
    #[arg(
        long = "seed-pairing",
        value_enum,
        default_value_t = SeedPairingMode::Strict
    )]
    pub pairing: SeedPairingMode,

    /// DEPRECATED (will be removed in a future release): legacy mismatch shorthand
    /// Set max mismatches (c) and min consecutive matches at seed start/end (p)
    /// These seeds will not overlap with perfect complementary seeds.
    /// Prefer --mismatch-max/--mismatch-prefix/--mismatch-suffix.
    #[arg(
        short = 'm',
        long = "mismatch",
        value_name = "c[:ps[:pe]]",
        default_value = "0:0",
        help_heading = "Deprecated"
    )]
    pub mismatch_legacy: MismatchSpec,

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

impl From<SeedConfig> for config::SeedConfig {
    fn from(value: SeedConfig) -> Self {
        let seed = if value.seed_start.is_some()
            || value.seed_end.is_some()
            || value.seed_length.is_some()
        {
            match (value.seed_start, value.seed_end) {
                (Some(start), Some(end)) => {
                    if let Some(length) = value.seed_length {
                        SeedSpec::IntervalWithLength { start, end, length }
                    } else {
                        SeedSpec::Interval { start, end }
                    }
                }
                (None, None) => SeedSpec::LengthOnly(value.seed_length.unwrap_or(6)),
                _ => value.seed_legacy,
            }
        } else {
            value.seed_legacy
        };
        config::SeedConfig {
            seed,
            no_guseed: value.no_guseed,
            pairing: value.pairing,
            mismatch: value.mismatch_legacy,
            mismatch_max: value.mismatch_max,
            mismatch_prefix: value.mismatch_prefix,
            mismatch_suffix: value.mismatch_suffix,
        }
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
    pub max_extension: u8,

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
        short = 'd',
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

impl From<ExtendArgs> for config::ExtendConfig {
    fn from(value: ExtendArgs) -> Self {
        config::ExtendConfig {
            max_extension: value.max_extension,
            delta_g: value.delta_g,
            matrix: value.matrix,
            penalty: value.penalty,
            seed_energy: value.seed_energy,
            no_max_prune: value.no_max_prune,
            dedup_shadow: value.dedup_shadow,
            band: value.band,
            matrix2: value.matrix2,
            matpath: value.matpath,
            temperature: value.temperature,
            weights: value.weights,
        }
    }
}

/// Options that apply to the `search` subcommand
#[derive(clap::Args, Debug, Clone)]
pub struct SearchArgs {
    #[command(flatten)]
    pub seed: SeedConfig,

    #[command(flatten)]
    pub extend: ExtendArgs,

    /// Output format
    #[arg(
        short = 'f',
        long = "format",
        value_name = "LEVEL",
        default_missing_value = "detailed",
        value_enum
    )]
    pub report_format: Option<OutputFormat>,

    /// DEPRECATED: Legacy argument for output format (1=detailed, 2=cigar, 3=binding_site, 4=minimal)
    #[arg(
        short = 'p',
        long = "report-alignment",
        value_name = "MODE",
        num_args = 0..=1,
        default_missing_value = "1",
        help_heading = "Deprecated"
    )]
    pub report_legacy: Option<u8>,

    /// Output compression codec (overrides file extension inference; gzip/gz, zstd/zst accepted)
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

impl From<SearchArgs> for config::SearchArgs {
    fn from(value: SearchArgs) -> Self {
        let format = if let Some(f) = value.report_format {
            f
        } else if let Some(legacy_mode) = value.report_legacy {
            warn!("-p/--report-alignment is deprecated. Use -f/--format instead.");
            match legacy_mode {
                1 => config::OutputFormat::Detailed,
                2 => config::OutputFormat::Cigar,
                3 => config::OutputFormat::BindingSite,
                4 => config::OutputFormat::Minimal,
                _ => config::OutputFormat::Detailed,
            }
        } else {
            config::OutputFormat::Detailed
        };

        config::SearchArgs {
            seed: value.seed.into(),
            extend: value.extend.into(),
            output: config::OutputConfig {
                format,
                compress: value.output_compress,
                level: value.output_level,
            },
            one_vs_one: value.one_vs_one,
            three_prime_match: value.three_prime_match,
            five_prime_match: value.five_prime_match,
        }
    }
}
