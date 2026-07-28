use clap::builder::{PossibleValuesParser, TypedValueParser};

use crate::config::{
    ScoreConfig, MAX_PENALTY_KCAL, MAX_TEMPERATURE_C, MIN_PENALTY_KCAL, MIN_TEMPERATURE_C,
};
use crate::dsm::DsmRegistry;
use crate::types::{DsmId, Energy};

fn parse_penalty(s: &str) -> Result<Energy, String> {
    let v: f64 = s.parse().map_err(|e| format!("{e}"))?;
    if !(MIN_PENALTY_KCAL..=MAX_PENALTY_KCAL).contains(&v) {
        return Err(format!(
            "penalty must be between {} and {}, got {v}",
            MIN_PENALTY_KCAL, MAX_PENALTY_KCAL
        ));
    }
    Energy::try_from(v)
}

fn parse_temperature(s: &str) -> Result<i32, String> {
    let v: i32 = s.parse().map_err(|e| format!("{e}"))?;
    if !(MIN_TEMPERATURE_C..=MAX_TEMPERATURE_C).contains(&v) {
        return Err(format!(
            "temperature must be between {} and {}, got {v}",
            MIN_TEMPERATURE_C, MAX_TEMPERATURE_C
        ));
    }
    Ok(v)
}

/// Arguments for global scoring model
#[derive(clap::Args, Debug, Clone)]
pub struct ScoreArgs {
    /// Dinucleotide stacking model for energy calculations
    #[arg(
        short = 'z',
        long = "matrix",
        value_name = "MATRIX",
        default_value_t = ScoreConfig::default().dsm_id,
        value_parser = PossibleValuesParser::new(DsmRegistry::all_names())
            .map(|name| DsmId(name))
    )]
    pub dsm_id: DsmId,

    /// Per-nucleotide penalty used by the scoring model (in kcal/mol, 0–50)
    #[arg(
        short = 'd',
        long = "penalty",
        value_name = "PENALTY",
        default_value_t = ScoreConfig::default().penalty,
        value_parser = parse_penalty
    )]
    pub penalty: Energy,

    /// Temperature for energy calculations (degrees Celsius, 0–100)
    #[arg(
        short = 'T',
        long = "temperature",
        value_name = "TEMP",
        default_value_t = ScoreConfig::default().temperature,
        value_parser = parse_temperature
    )]
    pub temperature: i32,
}

impl From<ScoreArgs> for ScoreConfig {
    fn from(value: ScoreArgs) -> Self {
        ScoreConfig {
            dsm_id: value.dsm_id,
            penalty: value.penalty,
            temperature: value.temperature,
        }
    }
}
