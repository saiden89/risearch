use crate::config::{ExtendConfig, OutputConfig, OutputFormat, SearchConfig};
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

impl SearchArgs {
    /// Resolve CLI-only output policy alongside the frontend-independent search
    /// configuration. Output format controls whether traceback is worth doing,
    /// but the output settings themselves do not belong in [`SearchConfig`].
    pub fn try_into_configs(self) -> Result<(SearchConfig, OutputConfig), Error> {
        let output: OutputConfig = self.output.try_into()?;
        let mut extend: ExtendConfig = self.extend.into();
        // Only the alignment-printing formats read the alignment back.
        extend.build_alignment = output.format != OutputFormat::Minimal;

        let search = SearchConfig {
            seed: self.seed.try_into()?,
            score: self.score.into(),
            extend,
            filter: self.filter.into(),
        };
        search.validate().map_err(Error::msg)?;

        Ok((search, output))
    }
}
