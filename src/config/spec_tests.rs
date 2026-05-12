#[cfg(test)]
mod tests {
    use crate::config::SeedConfig;


    #[test]
    fn length_only() {
        let cfg = SeedConfig { seed_length: Some(10), ..Default::default() };
        let res = cfg.resolve(100).expect("should normalize");
        assert_eq!(res, (1, 100, 10));
    }

    #[test]
    fn interval_only() {
        // seed_length = None in interval mode → use full interval width
        let cfg = SeedConfig { seed_start: Some(10), seed_end: Some(20), ..Default::default() };
        let res = cfg.resolve(100).expect("should normalize");
        // interval 10..20 inclusive -> length 11
        assert_eq!(res, (10, 20, 11));
    }

    #[test]
    fn interval_with_length() {
        let cfg = SeedConfig { seed_start: Some(10), seed_end: Some(20), seed_length: Some(5), ..Default::default() };
        let res = cfg.resolve(100).expect("should normalize");
        assert_eq!(res, (10, 20, 5));
    }

    #[test]
    fn interval_with_length_zero_uses_full_interval() {
        let cfg = SeedConfig { seed_start: Some(10), seed_end: Some(20), seed_length: Some(0), ..Default::default() };
        let res = cfg.resolve(100).expect("should normalize");
        assert_eq!(res, (10, 20, 11));
    }

    #[test]
    fn negative_interval() {
        let cfg = SeedConfig { seed_start: Some(-5), seed_end: Some(-1), ..Default::default() };
        let res = cfg.resolve(100).expect("should normalize");
        // -5 -> 96, -1 -> 100  (1-based), length = 5
        assert_eq!(res, (96, 100, 5));
    }

    #[test]
    fn length_zero_invalid() {
        let cfg = SeedConfig { seed_length: Some(0), ..Default::default() };
        let err = cfg.resolve(100).unwrap_err();
        assert!(err.contains("Invalid seed length"));
    }

    #[test]
    fn mixed_sign_interval_rejected() {
        let cfg = SeedConfig { seed_start: Some(5), seed_end: Some(-1), ..Default::default() };
        let err = cfg.resolve(100).unwrap_err();
        assert!(err.contains("Invalid seed interval"));
    }

    #[test]
    fn length_exceeds_interval_rejected() {
        let cfg = SeedConfig { seed_start: Some(10), seed_end: Some(12), seed_length: Some(5), ..Default::default() };
        let err = cfg.resolve(100).unwrap_err();
        assert!(err.contains("exceeds interval"));
    }
}
