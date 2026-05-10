use clap::builder::PossibleValuesParser;

use crate::config;
use crate::dsm::DSM_IDS;
use crate::types::{DsmId, Energy};

fn parse_penalty(s: &str) -> Result<Energy, String> {
    let v: f64 = s.parse().map_err(|e| format!("{e}"))?;
    if !(0.0..=50.0).contains(&v) {
        return Err(format!("penalty must be between 0 and 50, got {v}"));
    }
    Ok(Energy::from(v))
}

fn parse_temperature(s: &str) -> Result<i32, String> {
    let v: i32 = s.parse().map_err(|e| format!("{e}"))?;
    if !(0..=100).contains(&v) {
        return Err(format!("temperature must be between 0 and 100, got {v}"));
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
        default_value = "t04",
        value_parser = PossibleValuesParser::new(DSM_IDS)
    )]
    pub dsm_id: String,

    /// Per-nucleotide penalty used by the scoring model (in kcal/mol, 0–50)
    #[arg(
        short = 'd',
        long = "penalty",
        value_name = "PENALTY",
        default_value = "0.0",
        value_parser = parse_penalty
    )]
    pub penalty: Energy,

    /// Temperature for energy calculations (degrees Celsius, 0–100)
    #[arg(
        short = 'T',
        long = "temperature",
        value_name = "TEMP",
        default_value = "37",
        value_parser = parse_temperature
    )]
    pub temperature: i32,
}

impl From<ScoreArgs> for config::ScoreConfig {
    fn from(value: ScoreArgs) -> Self {
        config::ScoreConfig {
            // PossibleValuesParser guarantees dsm_id is a valid DSM_IDS member
            dsm_id: DsmId::try_from(value.dsm_id.as_str()).expect("clap validated"),
            penalty: value.penalty,
            temperature: value.temperature,
        }
    }
}
