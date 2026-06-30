use crate::config::FilterConfig;
use crate::types::Energy;

/// Arguments for filtering and pruning policies
#[derive(clap::Args, Debug, Clone)]
pub struct FilterArgs {
    /// Set deltaG energy threshold (in kcal/mol) to filter predictions
    #[arg(
        short = 'e',
        long = "energy",
        value_name = "dG",
        default_value = "-20.0",
        allow_hyphen_values = true
    )]
    pub total_energy: Energy,

    /// Energy per length threshold that filters seeds (in kcal/mol)
    #[arg(long = "seed-energy", value_name = "THRESHOLD", default_value = "0.0")]
    pub seed_energy: Energy,

    /// Disable maximality check (allows redundant seeds)
    #[arg(long = "no-max-prune", action = clap::ArgAction::SetTrue)]
    pub no_max_prune: bool,

    /// Report every maximal seed as its own hit instead of collapsing hits
    /// that share a final bounding box to the lowest-energy alignment
    #[arg(long = "no-dedup", action = clap::ArgAction::SetTrue)]
    pub no_dedup: bool,
}

impl From<FilterArgs> for FilterConfig {
    fn from(value: FilterArgs) -> Self {
        FilterConfig {
            delta_g: value.total_energy,
            seed_energy: value.seed_energy,
            no_max_prune: value.no_max_prune,
            no_dedup: value.no_dedup,
        }
    }
}
