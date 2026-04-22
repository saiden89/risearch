use crate::config;

use config::FilterConfig;
/// Arguments for filtering and pruning policies
#[derive(clap::Args, Debug, Clone)]
pub struct FilterArgs {
    /// Set deltaG energy threshold (in kcal/mol) to filter predictions
    #[arg(
        short = 'e',
        long = "energy",
        value_name = "dG",
        default_value_t = -20.0,
        allow_hyphen_values = true
    )]
    pub total_energy: f64,

    /// Energy per length threshold that filters seeds (in kcal/mol)
    #[arg(long = "seed-energy", value_name = "THRESHOLD", default_value_t = 0.0)]
    pub seed_energy: f64,

    /// Disable maximality check (allows redundant seeds)
    #[arg(long = "no-max-prune", action = clap::ArgAction::SetTrue)]
    pub no_max_prune: bool,
}

impl From<FilterArgs> for FilterConfig {
    fn from(value: FilterArgs) -> Self {
        FilterConfig {
            delta_g: value.total_energy,
            seed_energy: value.seed_energy,
            no_max_prune: value.no_max_prune,
        }
    }
}
