use crate::config;

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
    pub delta_g: f64,

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

impl From<FilterArgs> for config::FilterConfig {
    fn from(value: FilterArgs) -> Self {
        config::FilterConfig {
            delta_g: value.delta_g,
            seed_energy: value.seed_energy,
            no_max_prune: value.no_max_prune,
            dedup_shadow: value.dedup_shadow,
        }
    }
}
