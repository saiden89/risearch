use crate::config::{self, MismatchSpec, SeedSpec};
use crate::types::SeedPairingMode;
use std::str::FromStr;

const DEFAULT_SEED_LEN: i64 = 6;

/// Legacy CLI parser boundary for `-m/--mismatch`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CliMismatchSpec(pub MismatchSpec);

impl FromStr for CliMismatchSpec {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        if s.is_empty() {
            return Err("empty mismatch spec".into());
        }

        let parts: Vec<&str> = s.split(':').collect();
        let (max_mismatches, min_prefix, min_suffix) = match parts.len() {
            1 => {
                let max = parts[0]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid max mismatches: {}", e))?;
                (max, max, max)
            }
            2 => {
                let max = parts[0]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid max mismatches: {}", e))?;
                let min = parts[1]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid min consecutive: {}", e))?;
                (max, min, min)
            }
            3 => {
                let max = parts[0]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid max mismatches: {}", e))?;
                let min_prefix = parts[1]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid min prefix matches: {}", e))?;
                let min_suffix = parts[2]
                    .parse::<usize>()
                    .map_err(|e| format!("invalid min suffix matches: {}", e))?;
                (max, min_prefix, min_suffix)
            }
            _ => {
                return Err(format!(
                    "invalid mismatch spec '{}': expected 'c', 'c:p', or 'c:ps:pe' format",
                    s
                ));
            }
        };

        Ok(Self(MismatchSpec {
            max_mismatches,
            min_prefix_matches: min_prefix,
            min_suffix_matches: min_suffix,
        }))
    }
}

impl From<CliMismatchSpec> for MismatchSpec {
    fn from(value: CliMismatchSpec) -> Self {
        value.0
    }
}

/// Legacy CLI parser boundary for `-s/--seed`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CliSeedSpec(pub SeedSpec);

impl FromStr for CliSeedSpec {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        if s.is_empty() {
            return Err("empty seed spec".into());
        }

        let parse_i64 = |val: &str, field: &str| -> Result<i64, String> {
            val.parse::<i64>()
                .map_err(|e| format!("bad {}: {}", field, e))
        };

        match s.find(':') {
            None => {
                let len = parse_i64(s, "length")?;
                Ok(Self(SeedSpec::LengthOnly(len)))
            }
            Some(colon_idx) => {
                let (start_str, rest) = s.split_at(colon_idx);
                let rest = &rest[1..];

                if rest.is_empty() {
                    return Err("missing end in interval".into());
                }

                let start = parse_i64(start_str, "start")?;

                match rest.find('/') {
                    None => {
                        let end = parse_i64(rest, "end")?;
                        Ok(Self(SeedSpec::Interval { start, end }))
                    }
                    Some(slash_idx) => {
                        let (end_str, len_part) = rest.split_at(slash_idx);
                        let len_str = &len_part[1..];

                        if end_str.is_empty() {
                            return Err("missing end in interval".into());
                        }
                        if len_str.is_empty() {
                            return Err("missing length after '/'".into());
                        }

                        let end = parse_i64(end_str, "end")?;
                        let length = parse_i64(len_str, "length")?;
                        Ok(Self(SeedSpec::IntervalWithLength { start, end, length }))
                    }
                }
            }
        }
    }
}

impl From<CliSeedSpec> for SeedSpec {
    fn from(value: CliSeedSpec) -> Self {
        value.0
    }
}

/// Arguments for seed generation
#[derive(clap::Args, Debug, Clone)]
pub struct SeedConfig {
    /// DEPRECATED (will be removed in a future release): legacy seed spec
    /// Formats: "l", "m:n", "m:n/l"
    #[arg(
        short = 's',
        long = "seed",
        value_name = "start:end/length",
        default_value = "6",
        help_heading = "Deprecated"
    )]
    pub seed_legacy: CliSeedSpec,

    /// Seed interval start (1-based, can be negative)
    /// TODO: Consider explicit one-sided bounds (e.g. --seed-to-end/--seed-from-start)
    /// instead of inferring missing start/end. Keep strict parsing for now.
    #[arg(
        long = "seed-start",
        value_name = "START",
        allow_hyphen_values = true,
        requires = "seed_end"
    )]
    pub seed_start: Option<i64>,

    /// Seed interval end (1-based, can be negative)
    #[arg(
        long = "seed-end",
        value_name = "END",
        allow_hyphen_values = true,
        requires = "seed_start"
    )]
    pub seed_end: Option<i64>,

    /// Seed length (use alone, or with seed-start/seed-end to constrain interval)
    #[arg(long = "seed-length", value_name = "LENGTH")]
    pub seed_length: Option<i64>,

    /// DEPRECATED (will be removed in a future release): disable G-U wobble pairs within the seed
    #[arg(
        short = 'U',
        long = "no-guseed",
        alias = "noGUseed",
        action = clap::ArgAction::SetTrue,
        help_heading = "Deprecated"
    )]
    pub no_guseed: bool,

    /// Seed pairing mode (allow_wobble or strict)
    #[arg(
        long = "seed-pairing",
        value_enum,
        default_value_t = SeedPairingMode::Strict
    )]
    pub pairing: SeedPairingMode,

    /// DEPRECATED (will be removed in a future release): legacy mismatch shorthand
    /// Set max mismatches (c) and min consecutive matches at seed start/end (p)
    /// These seeds will not overlap with perfect complementary seeds.
    /// Prefer --mismatch-max/--mismatch-prefix/--mismatch-suffix.
    #[arg(
        short = 'm',
        long = "mismatch",
        value_name = "c[:ps[:pe]]",
        default_value = "0:0",
        help_heading = "Deprecated"
    )]
    pub mismatch_legacy: CliMismatchSpec,

    /// Max number of mismatches allowed in the seed (preferred)
    #[arg(long = "mismatch-max", value_name = "C")]
    pub mismatch_max: Option<usize>,

    /// Min consecutive matches at seed start (prefix / 5')
    #[arg(long = "mismatch-prefix", value_name = "PS")]
    pub mismatch_prefix: Option<usize>,

    /// Min consecutive matches at seed end (suffix / 3')
    #[arg(long = "mismatch-suffix", value_name = "PE")]
    pub mismatch_suffix: Option<usize>,

}

impl From<SeedConfig> for config::SeedConfig {
    fn from(value: SeedConfig) -> Self {
        let legacy_seed: SeedSpec = value.seed_legacy.into();
        let legacy_mismatch: MismatchSpec = value.mismatch_legacy.into();
        let seed = if value.seed_start.is_some()
            || value.seed_end.is_some()
            || value.seed_length.is_some()
        {
            match (value.seed_start, value.seed_end) {
                (Some(start), Some(end)) => {
                    if let Some(length) = value.seed_length {
                        SeedSpec::IntervalWithLength { start, end, length }
                    } else {
                        SeedSpec::Interval { start, end }
                    }
                }
                (None, None) => SeedSpec::LengthOnly(value.seed_length.unwrap_or(DEFAULT_SEED_LEN)),
                _ => legacy_seed,
            }
        } else {
            legacy_seed
        };

        config::SeedConfig {
            seed,
            no_guseed: value.no_guseed,
            pairing: value.pairing,
            mismatch: legacy_mismatch,
            mismatch_max: value.mismatch_max,
            mismatch_prefix: value.mismatch_prefix,
            mismatch_suffix: value.mismatch_suffix,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{CliMismatchSpec, CliSeedSpec};
    use crate::config::{MismatchSpec, SeedSpec};
    use std::str::FromStr;

    #[test]
    fn parses_seed_cli_specs() {
        assert_eq!(
            CliSeedSpec::from_str("10").expect("parse").0,
            SeedSpec::LengthOnly(10)
        );
        assert_eq!(
            CliSeedSpec::from_str("10:20").expect("parse").0,
            SeedSpec::Interval { start: 10, end: 20 }
        );
        assert_eq!(
            CliSeedSpec::from_str("10:20/5").expect("parse").0,
            SeedSpec::IntervalWithLength {
                start: 10,
                end: 20,
                length: 5
            }
        );
        assert!(CliSeedSpec::from_str("abc").is_err());
    }

    #[test]
    fn parses_mismatch_cli_specs() {
        assert_eq!(
            CliMismatchSpec::from_str("1").expect("parse").0,
            MismatchSpec::new(1, 1, 1)
        );
        assert_eq!(
            CliMismatchSpec::from_str("1:3").expect("parse").0,
            MismatchSpec::new(1, 3, 3)
        );
        assert_eq!(
            CliMismatchSpec::from_str("1:3:5").expect("parse").0,
            MismatchSpec::new(1, 3, 5)
        );
        assert!(CliMismatchSpec::from_str("1:2:3:4").is_err());
    }
}
