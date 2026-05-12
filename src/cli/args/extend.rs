use crate::config;
use config::ExtendConfig;

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
}

impl From<ExtendArgs> for ExtendConfig {
    fn from(value: ExtendArgs) -> Self {
        ExtendConfig {
            max_extension: value.max_extension,
        }
    }
}
