#[cfg(test)]
mod tests {
    use super::super::SeedSpec;

    #[test]
    fn test_length_only() {
        let res = SeedSpec::LengthOnly(10).normalize(100).expect("should normalize");
        assert_eq!(res, (1, 100, 10));
    }

    #[test]
    fn test_interval_only() {
        let res = SeedSpec::Interval { start: 10, end: 20 }
            .normalize(100)
            .expect("should normalize");
        // interval 10..20 inclusive -> length 11
        assert_eq!(res, (10, 20, 11));
    }

    #[test]
    fn test_interval_with_length() {
        let res = SeedSpec::IntervalWithLength {
            start: 10,
            end: 20,
            length: 5,
        }
        .normalize(100)
        .expect("should normalize");
        assert_eq!(res, (10, 20, 5));
    }

    #[test]
    fn test_negative_interval() {
        let res = SeedSpec::Interval { start: -5, end: -1 }
            .normalize(100)
            .expect("should normalize");
        // -5 -> 96, -1 -> 100  (1-based), length = 5
        assert_eq!(res, (96, 100, 5));
    }

    #[test]
    fn test_invalid_length_zero() {
        let err = SeedSpec::LengthOnly(0).normalize(100).unwrap_err();
        assert!(err.contains("Invalid seed length"));
    }

    #[test]
    fn test_mixed_sign_interval() {
        let err = SeedSpec::Interval { start: 5, end: -1 }
            .normalize(100)
            .unwrap_err();
        assert!(err.contains("Invalid seed interval"));
    }

    #[test]
    fn test_length_exceeds_interval() {
        let err = SeedSpec::IntervalWithLength {
            start: 10,
            end: 12,
            length: 5,
        }
        .normalize(100)
        .unwrap_err();
        assert!(err.contains("exceeds interval"));
    }
}
