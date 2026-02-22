use crate::config::{self, OutputCompression, OutputFormat, SearchAxis};

use super::{ExtendArgs, FilterArgs, ScoreArgs, SeedConfig};

/// Options that apply to the `search` subcommand
#[derive(clap::Args, Debug, Clone)]
pub struct SearchArgs {
    #[command(flatten)]
    pub seed: SeedConfig,

    #[command(flatten)]
    pub score: ScoreArgs,

    #[command(flatten)]
    pub extend: ExtendArgs,

    #[command(flatten)]
    pub filter: FilterArgs,

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

    /// Search parallelization axis override.
    #[arg(long = "search-axis", value_enum, default_value_t = SearchAxis::Auto)]
    pub search_axis: SearchAxis,

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
            score: value.score.into(),
            extend: value.extend.into(),
            filter: value.filter.into(),
            output: config::OutputConfig {
                format,
                compress: value.output_compress,
                level: value.output_level,
            },
            axis: value.search_axis,
            one_vs_one: value.one_vs_one,
            three_prime_match: value.three_prime_match,
            five_prime_match: value.five_prime_match,
        }
    }
}
