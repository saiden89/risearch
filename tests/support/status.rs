//! Status types for parity comparison.
//!
//! Defines the outcome states when comparing Rust vs C hits.

use std::sync::LazyLock;

// =============================================================================
// HIT STATUS
// =============================================================================

/// Status of a hit in parity comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum HitStatus {
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
pub(crate) enum MissingReason {
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
    pub(crate) fn is_acceptable(&self, mode: ParityMode) -> bool {
        match mode {
            ParityMode::Absolute => false,
            ParityMode::Strict => false,
            ParityMode::Balanced => false, // No missing allowed in balanced mode
            ParityMode::Relaxed => matches!(self, Self::BetterEnergy | Self::EqualEnergy),
        }
    }
}

// =============================================================================
// PARITY MODE
// =============================================================================

/// Comparison mode for parity tests.
///
/// | Mode     | Identical | CoOptimal | RustBetter | Extras | Missing |
/// |----------|-----------|-----------|------------|--------|---------|
/// | Absolute | ✓         | ✗         | ✗          | ✗      | ✗       |
/// | Strict   | ✓         | ✓         | ✗          | ✗      | ✗       |
/// | Balanced | ✓         | ✓         | ✗          | ✓      | ✗       |
/// | Relaxed  | ✓         | ✓         | ✓          | ✓      | ✗*      |
///
/// *Missing is acceptable in Relaxed only if covered by better/equal energy hit
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) enum ParityMode {
    /// 100% identical to C - no differences at all
    Absolute,
    /// Allow co-optimal alignments (same energy, different trace)
    Strict,
    /// Allow co-optimal and extras, but no missing or rust-better
    #[default]
    Balanced,
    /// Accept improvements: co-optimal, better energy, extras, covered missings
    Relaxed,
}

/// Single source of truth: the parity mode used for all tests.
///
/// Set via environment variable `PARITY_MODE`:
/// - `absolute` - 100% identical to C
/// - `strict` (default) - allow co-optimal alignments
/// - `balanced` - allow co-optimal and extras, but no missing
/// - `relaxed` - accept all improvements including covered missings
///
/// Example: `PARITY_MODE=balanced cargo test --test c_parity`
pub(crate) static TEST_PARITY_MODE: LazyLock<ParityMode> = LazyLock::new(|| {
    match std::env::var("PARITY_MODE").as_deref() {
        Ok("absolute") => ParityMode::Absolute,
        Ok("balanced") => ParityMode::Balanced,
        Ok("relaxed") => ParityMode::Relaxed,
        _ => ParityMode::Balanced, // Default: allow extras, reject missing
    }
});
