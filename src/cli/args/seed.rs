use crate::config::{self, MismatchSpec, SeedSpec};
use std::str::FromStr;

const DEFAULT_SEED_LEN: i64 = 6;

pub fn seed_spec_from_args(
    seed_start: Option<i64>,
    seed_end: Option<i64>,
    seed_length: Option<i64>,
) -> Result<SeedSpec, String> {
    match (seed_start, seed_end, seed_length) {
        (Some(start), Some(end), length) => Ok(SeedSpec::Interval { start, end, length }),
        (None, None, Some(length)) => Ok(SeedSpec::LengthOnly(length)),
        (None, None, None) => Ok(SeedSpec::LengthOnly(DEFAULT_SEED_LEN)),
        _ => Err(
            "use seed_length alone, or seed_start + seed_end (optionally with seed_length)".into(),
        ),
    }
}

/// Legacy CLI parser boundary for `-m/--mismatch`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LegacyMismatchSpec(pub MismatchSpec);

impl FromStr for LegacyMismatchSpec {
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

/// Legacy CLI parser boundary for `-s/--seed`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LegacySeedSpec(pub SeedSpec);

impl FromStr for LegacySeedSpec {
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
                        Ok(Self(SeedSpec::Interval { start, end, length: None }))
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
                        Ok(Self(SeedSpec::Interval { start, end, length: Some(length) }))
                    }
                }
            }
        }
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
        conflicts_with_all = ["seed_start", "seed_end", "seed_length"],
        help_heading = "Deprecated"
    )]
    pub seed_legacy: Option<LegacySeedSpec>,

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

    /// Disable G-U wobble pairs when locating and maximizing seeds
    #[arg(long = "no-seed-wobble", action = clap::ArgAction::SetTrue)]
    pub no_seed_wobble: bool,

    /// DEPRECATED: legacy C flag for --no-seed-wobble
    #[arg(long = "noGUseed", hide = true, action = clap::ArgAction::SetTrue)]
    pub no_guseed_legacy: bool,

    /// DEPRECATED (will be removed in a future release): legacy mismatch shorthand
    /// Set max mismatches (c) and min consecutive matches at seed start/end (p)
    /// These seeds will not overlap with perfect complementary seeds.
    /// Prefer --mismatch-max/--mismatch-prefix/--mismatch-suffix.
    #[arg(
        short = 'm',
        long = "mismatch",
        value_name = "c[:ps[:pe]]",
        conflicts_with_all = ["mismatch_max", "mismatch_prefix", "mismatch_suffix"],
        help_heading = "Deprecated"
    )]
    pub mismatch_legacy: Option<LegacyMismatchSpec>,

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
        let seed = if value.seed_start.is_some()
            || value.seed_end.is_some()
            || value.seed_length.is_some()
        {
            seed_spec_from_args(value.seed_start, value.seed_end, value.seed_length)
                .expect("seed arguments should be validated by clap")
        } else {
            value
                .seed_legacy
                .map(|s| s.0)
                .unwrap_or(SeedSpec::LengthOnly(DEFAULT_SEED_LEN))
        };
        let mismatch = value.mismatch_legacy.map_or_else(
            || {
                MismatchSpec::new(
                    value.mismatch_max.unwrap_or(0),
                    value.mismatch_prefix.unwrap_or(0),
                    value.mismatch_suffix.unwrap_or(0),
                )
            },
            |m| m.0,
        );

        config::SeedConfig {
            seed,
            seed_wobble: !(value.no_seed_wobble || value.no_guseed_legacy),
            mismatch,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{seed_spec_from_args, LegacyMismatchSpec, LegacySeedSpec};
    use crate::config::{MismatchSpec, SeedSpec};
    use std::str::FromStr;

    #[test]
    fn parses_seed_cli_specs() {
        assert_eq!(
            LegacySeedSpec::from_str("10").expect("parse").0,
            SeedSpec::LengthOnly(10)
        );
        assert_eq!(
            LegacySeedSpec::from_str("10:20").expect("parse").0,
            SeedSpec::Interval { start: 10, end: 20, length: None }
        );
        assert_eq!(
            LegacySeedSpec::from_str("10:20/5").expect("parse").0,
            SeedSpec::Interval { start: 10, end: 20, length: Some(5) }
        );
        assert!(LegacySeedSpec::from_str("abc").is_err());
    }

    #[test]
    fn parses_mismatch_cli_specs() {
        assert_eq!(
            LegacyMismatchSpec::from_str("1").expect("parse").0,
            MismatchSpec::new(1, 1, 1)
        );
        assert_eq!(
            LegacyMismatchSpec::from_str("1:3").expect("parse").0,
            MismatchSpec::new(1, 3, 3)
        );
        assert_eq!(
            LegacyMismatchSpec::from_str("1:3:5").expect("parse").0,
            MismatchSpec::new(1, 3, 5)
        );
        assert!(LegacyMismatchSpec::from_str("1:2:3:4").is_err());
    }

    #[test]
    fn builds_seed_specs_from_named_args() {
        assert_eq!(
            seed_spec_from_args(None, None, None).expect("parse"),
            SeedSpec::LengthOnly(6)
        );
        assert_eq!(
            seed_spec_from_args(Some(3), Some(12), None).expect("parse"),
            SeedSpec::Interval { start: 3, end: 12, length: None }
        );
        assert_eq!(
            seed_spec_from_args(Some(3), Some(12), Some(7)).expect("parse"),
            SeedSpec::Interval { start: 3, end: 12, length: Some(7) }
        );
    }

    #[test]
    fn rejects_partial_named_seed_args() {
        assert!(seed_spec_from_args(Some(3), None, None).is_err());
        assert!(seed_spec_from_args(None, Some(12), None).is_err());
    }
}
