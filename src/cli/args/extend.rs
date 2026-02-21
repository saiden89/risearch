use crate::config;

/// Arguments for seed extension strategy
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

    // ========================================================================
    // TODO: Placeholder flags from C implementation - not yet implemented
    // ========================================================================
    /// TODO: Banded search - limits the search for bulged matches.
    /// In C: `-b band, --band=band` - Integer size of bands limiting bulge search.
    /// The minimum size is 1; use seed option to avoid any bulge.
    #[arg(long = "band", value_name = "BAND", hide = true)]
    pub band: Option<u32>,
}

impl From<ExtendArgs> for config::ExtendConfig {
    fn from(value: ExtendArgs) -> Self {
        config::ExtendConfig {
            max_extension: value.max_extension,
            band: value.band,
        }
    }
}
