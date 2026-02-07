#[cfg(test)]
mod tests {
    use super::super::SeedSpec;
    use std::str::FromStr;

    // helper: parse and normalize, returning the tuple or panic with message
    fn parse_and_norm(s: &str, qlen: usize) -> Result<(usize, usize, usize), String> {
        let spec = SeedSpec::from_str(s).map_err(|e| format!("parse err: {}", e))?;
        spec.normalize(qlen)
    }

    #[test]
    fn test_length_only() {
        let res = parse_and_norm("10", 100).expect("should parse");
        assert_eq!(res, (1, 100, 10));
    }

    #[test]
    fn test_interval_only() {
        let res = parse_and_norm("10:20", 100).expect("should parse");
        // interval 10..20 inclusive -> length 11
        assert_eq!(res, (10, 20, 11));
    }

    #[test]
    fn test_interval_with_length() {
        let res = parse_and_norm("10:20/5", 100).expect("should parse");
        assert_eq!(res, (10, 20, 5));
    }

    #[test]
    fn test_negative_interval() {
        let res = parse_and_norm("-5:-1", 100).expect("should parse");
        // -5 -> 96, -1 -> 100  (1-based), length = 5
        assert_eq!(res, (96, 100, 5));
    }

    #[test]
    fn test_invalid_length_zero() {
        let err = parse_and_norm("0", 100).unwrap_err();
        assert!(err.contains("Invalid seed length"));
    }

    #[test]
    fn test_mixed_sign_interval() {
        let err = parse_and_norm("5:-1", 100).unwrap_err();
        assert!(err.contains("Invalid seed interval"));
    }

    #[test]
    fn test_length_exceeds_interval() {
        let err = parse_and_norm("10:12/5", 100).unwrap_err();
        assert!(err.contains("exceeds interval"));
    }

    #[test]
    fn test_bad_parse() {
        assert!(SeedSpec::from_str("abc").is_err());
    }
}
