//! Status types for parity comparison.
//!
//! Defines the outcome states when comparing Rust vs C hits.

// =============================================================================
// HIT STATUS
// =============================================================================

/// Status of a hit in parity comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HitStatus {
    /// Exact match between Rust and C (markers already stripped)
    Identical,
    /// Equal energy, different trace (co-optimal alignment)
    CoOptimal,
    /// Rust has better (lower) energy
    RustBetter,
    /// Rust has worse (higher) energy
    RustWorse,
    /// Only in Rust output
    #[allow(dead_code)] // Part of complete status enum
    Extra, // Only in RustC output
    /// Only in C output
    Missing,
}

impl std::fmt::Display for HitStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Identical => write!(f, "IDENTICAL"),
            Self::CoOptimal => write!(f, "CO-OPTIMAL"),
            Self::RustBetter => write!(f, "RUST BETTER"),
            Self::RustWorse => write!(f, "RUST WORSE"),
            Self::Extra => write!(f, "EXTRA"),
            Self::Missing => write!(f, "MISSING"),
        }
    }
}

// =============================================================================
// MISSING REASON
// =============================================================================

/// Why a C hit was not found in Rust output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MissingReason {
    /// Rust found overlapping hit with better (lower) energy
    BetterEnergy,
    /// Rust found overlapping hit with equal energy (co-optimal)
    EqualEnergy,
    /// Rust found overlapping hit with worse (higher) energy
    WorseEnergy,
    /// No overlapping Rust hit - completely missed
    NoOverlap,
}

impl std::fmt::Display for MissingReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::BetterEnergy => write!(f, "COVERED_BETTER"),
            Self::EqualEnergy => write!(f, "COVERED_EQUAL"),
            Self::WorseEnergy => write!(f, "COVERED_WORSE"),
            Self::NoOverlap => write!(f, "NO_OVERLAP"),
        }
    }
}

impl MissingReason {
    /// Check if this reason is acceptable given the parity mode.
    pub fn is_acceptable(&self, mode: ParityMode) -> bool {
        match mode {
            ParityMode::Absolute => false,
            ParityMode::Strict => false,
            ParityMode::Relaxed => matches!(self, Self::BetterEnergy | Self::EqualEnergy),
        }
    }
}

// =============================================================================
// PARITY MODE
// =============================================================================

/// Comparison mode for parity tests.
///
/// | Mode     | Identical | CoOptimal | RustBetter | RustWorse | Extra | Missing |
/// |----------|-----------|-----------|------------|-----------|-------|---------|
/// | Absolute | ✓         | ✗         | ✗          | ✗         | ✗     | ✗       |
/// | Strict   | ✓         | ✓         | ✗          | ✗         | ✗     | ✗       |
/// | Relaxed  | ✓         | ✓         | ✓          | ✗         | ✓     | ✗*      |
///
/// *Missing is acceptable in Relaxed only if covered by better/equal energy hit
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ParityMode {
    /// 100% identical to C - no differences at all
    Absolute,
    /// Allow co-optimal alignments (same energy, different trace)
    #[default]
    Strict,
    /// Accept improvements: co-optimal, better energy, extras
    Relaxed,
}

/// Single source of truth: the parity mode used for all tests.
/// Change this to adjust what level of parity is required across the test suite.
pub const TEST_PARITY_MODE: ParityMode = ParityMode::Strict;

// =============================================================================
// UNIT TESTS
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_missing_reason_acceptable_absolute() {
        // Absolute mode: nothing is acceptable
        assert!(!MissingReason::BetterEnergy.is_acceptable(ParityMode::Absolute));
        assert!(!MissingReason::EqualEnergy.is_acceptable(ParityMode::Absolute));
        assert!(!MissingReason::WorseEnergy.is_acceptable(ParityMode::Absolute));
        assert!(!MissingReason::NoOverlap.is_acceptable(ParityMode::Absolute));
    }

    #[test]
    fn test_missing_reason_acceptable_strict() {
        // Strict mode: nothing is acceptable (co-optimal is about hits, not missings)
        assert!(!MissingReason::BetterEnergy.is_acceptable(ParityMode::Strict));
        assert!(!MissingReason::EqualEnergy.is_acceptable(ParityMode::Strict));
        assert!(!MissingReason::WorseEnergy.is_acceptable(ParityMode::Strict));
        assert!(!MissingReason::NoOverlap.is_acceptable(ParityMode::Strict));
    }

    #[test]
    fn test_missing_reason_acceptable_relaxed() {
        // Relaxed mode: better/equal energy covered is acceptable
        assert!(MissingReason::BetterEnergy.is_acceptable(ParityMode::Relaxed));
        assert!(MissingReason::EqualEnergy.is_acceptable(ParityMode::Relaxed));
        assert!(!MissingReason::WorseEnergy.is_acceptable(ParityMode::Relaxed));
        assert!(!MissingReason::NoOverlap.is_acceptable(ParityMode::Relaxed));
    }

    #[test]
    fn test_hit_status_display() {
        assert_eq!(format!("{}", HitStatus::RustBetter), "RUST BETTER");
        assert_eq!(format!("{}", HitStatus::Missing), "MISSING");
        assert_eq!(format!("{}", HitStatus::CoOptimal), "CO-OPTIMAL");
    }

    #[test]
    fn test_parity_mode_default() {
        assert_eq!(ParityMode::default(), ParityMode::Strict);
    }
}
