#[cfg(test)]
mod tests {
    use crate::config::{
        ExtendConfig, ScoreConfig, SearchConfig, SeedConfig, DEFAULT_SEED_LEN, MAX_EXTENSION,
    };
    use crate::types::{DsmId, Energy};

    #[test]
    fn defaults_are_the_canonical_frontend_defaults() {
        let cfg = SearchConfig::default();

        assert_eq!(cfg.seed.seed_length, None);
        assert_eq!(DEFAULT_SEED_LEN, 6);
        assert!(!cfg.seed.seed_wobble);
        assert_eq!(cfg.seed.max_mismatches, 0);
        assert_eq!(cfg.seed.min_prefix_matches, 1);
        assert_eq!(cfg.seed.min_suffix_matches, 0);
        assert_eq!(cfg.score.dsm_id, DsmId::from("t04"));
        assert_eq!(cfg.score.penalty, Energy::default());
        assert_eq!(cfg.score.temperature, 37);
        assert_eq!(cfg.extend.max_extension, 20);
        assert!(cfg.extend.build_alignment);
        cfg.validate().expect("default config");
    }

    #[test]
    fn length_only() {
        let cfg = SeedConfig {
            seed_length: Some(10),
            ..Default::default()
        };
        let res = cfg.resolve(100).expect("should normalize");
        assert_eq!(res, (1, 100, 10));
    }

    #[test]
    fn interval_only() {
        // seed_length = None in interval mode → use full interval width
        let cfg = SeedConfig {
            seed_start: Some(10),
            seed_end: Some(20),
            ..Default::default()
        };
        let res = cfg.resolve(100).expect("should normalize");
        // interval 10..20 inclusive -> length 11
        assert_eq!(res, (10, 20, 11));
    }

    #[test]
    fn interval_with_length() {
        let cfg = SeedConfig {
            seed_start: Some(10),
            seed_end: Some(20),
            seed_length: Some(5),
            ..Default::default()
        };
        let res = cfg.resolve(100).expect("should normalize");
        assert_eq!(res, (10, 20, 5));
    }

    #[test]
    fn interval_with_length_zero_uses_full_interval() {
        let cfg = SeedConfig {
            seed_start: Some(10),
            seed_end: Some(20),
            seed_length: Some(0),
            ..Default::default()
        };
        let res = cfg.resolve(100).expect("should normalize");
        assert_eq!(res, (10, 20, 11));
    }

    #[test]
    fn negative_interval() {
        let cfg = SeedConfig {
            seed_start: Some(-5),
            seed_end: Some(-1),
            ..Default::default()
        };
        let res = cfg.resolve(100).expect("should normalize");
        // -5 -> 96, -1 -> 100  (1-based), length = 5
        assert_eq!(res, (96, 100, 5));
    }

    #[test]
    fn length_zero_invalid() {
        let cfg = SeedConfig {
            seed_length: Some(0),
            ..Default::default()
        };
        let err = cfg.resolve(100).unwrap_err();
        assert!(err.contains("Invalid seed length"));
    }

    #[test]
    fn mixed_sign_interval_rejected() {
        let cfg = SeedConfig {
            seed_start: Some(5),
            seed_end: Some(-1),
            ..Default::default()
        };
        let err = cfg.resolve(100).unwrap_err();
        assert!(err.contains("Invalid seed interval"));
    }

    #[test]
    fn length_exceeds_interval_rejected() {
        let cfg = SeedConfig {
            seed_start: Some(10),
            seed_end: Some(12),
            seed_length: Some(5),
            ..Default::default()
        };
        let err = cfg.resolve(100).unwrap_err();
        assert!(err.contains("exceeds interval"));
    }

    #[test]
    fn query_independent_seed_invariants_are_validated() {
        let partial = SeedConfig {
            seed_start: Some(1),
            ..Default::default()
        };
        assert!(partial.validate().unwrap_err().contains("together"));

        let zero_coordinate = SeedConfig {
            seed_start: Some(0),
            seed_end: Some(5),
            ..Default::default()
        };
        assert!(zero_coordinate
            .validate()
            .unwrap_err()
            .contains("cannot be zero"));

        let reversed = SeedConfig {
            seed_start: Some(5),
            seed_end: Some(1),
            ..Default::default()
        };
        assert!(reversed.validate().unwrap_err().contains("precedes start"));
    }

    #[test]
    fn score_validation_does_not_depend_on_clap() {
        let cfg = ScoreConfig {
            penalty: Energy::from_kcal(-0.1),
            ..Default::default()
        };
        assert!(cfg.validate().unwrap_err().contains("penalty"));

        let cfg = ScoreConfig {
            penalty: Energy::from_kcal(50.1),
            ..Default::default()
        };
        assert!(cfg.validate().unwrap_err().contains("penalty"));

        let cfg = ScoreConfig {
            temperature: -1,
            ..Default::default()
        };
        assert!(cfg.validate().unwrap_err().contains("temperature"));

        let cfg = ScoreConfig {
            temperature: 101,
            ..Default::default()
        };
        assert!(cfg.validate().unwrap_err().contains("temperature"));

        let cfg = ScoreConfig {
            dsm_id: DsmId::from("not-a-dsm"),
            ..Default::default()
        };
        assert!(cfg.validate().unwrap_err().contains("unknown DSM"));
    }

    #[test]
    fn extension_validation_accepts_only_the_documented_sentinel_and_range() {
        for value in [-1, 0, MAX_EXTENSION] {
            ExtendConfig {
                max_extension: value,
                ..Default::default()
            }
            .validate()
            .expect("valid extension");
        }

        for value in [-2, MAX_EXTENSION + 1] {
            assert!(ExtendConfig {
                max_extension: value,
                ..Default::default()
            }
            .validate()
            .is_err());
        }
    }
}
