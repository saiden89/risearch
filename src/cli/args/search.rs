use crate::config;
use anyhow::Error;

use super::{ExtendArgs, FilterArgs, InputArgs, OutputArgs, ScoreArgs, SeedArgs};

/// Arguments that apply to the `search` subcommand
#[derive(clap::Args, Debug, Clone)]
pub struct SearchArgs {
    #[command(flatten)]
    pub input: InputArgs,

    #[command(flatten)]
    pub seed: SeedArgs,

    #[command(flatten)]
    pub score: ScoreArgs,

    #[command(flatten)]
    pub extend: ExtendArgs,

    #[command(flatten)]
    pub filter: FilterArgs,

    #[command(flatten)]
    pub output: OutputArgs,
}

impl TryFrom<SearchArgs> for config::SearchConfig {
    type Error = Error;

    fn try_from(value: SearchArgs) -> Result<Self, Self::Error> {
        Ok(config::SearchConfig {
            seed: value.seed.try_into()?,
            score: value.score.into(),
            extend: value.extend.into(),
            filter: value.filter.into(),
            output: value.output.try_into()?,
        })
    }
}
