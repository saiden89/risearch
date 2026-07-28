//! Core domain types used across the codebase

use serde::{Deserialize, Serialize};

/// Nucleotide/gap representation for DSM indexing and sequence operations.
///
/// The discriminants are the canonical internal rank. DSM tables, DP scoring,
/// and suffix-array construction/partitioning all depend on this exact order.
#[repr(u8)]
#[derive(
    Clone,
    Copy,
    Debug,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Default,
    Serialize,
    Deserialize,
    rkyv::Archive,
    rkyv::Serialize,
    rkyv::Deserialize,
)]
pub enum Base {
    #[default]
    Gap = 0,
    A = 1,
    C = 2,
    G = 3,
    N = 4,
    U = 5,
}

/// Base → uppercase ASCII byte
static BASE_TO_UPPER: [u8; 6] = *b"-ACGNU";

/// Base → lowercase ASCII byte (for to_byte())
static BASE_TO_BYTE: [u8; 6] = *b"-acgnu";

impl TryFrom<char> for Base {
    type Error = String;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c.to_ascii_uppercase() {
            'A' => Ok(Base::A),
            'C' => Ok(Base::C),
            'G' => Ok(Base::G),
            'U' | 'T' => Ok(Base::U),
            'N' => Ok(Base::N),
            '-' => Ok(Base::Gap),
            _ => Err(format!("invalid base '{c}'")),
        }
    }
}

impl Base {
    /// Convert to array index
    #[inline]
    pub const fn as_usize(self) -> usize {
        self as usize
    }

    /// Convert raw u8 byte representation back to Base.
    #[inline]
    pub fn from_u8(v: u8) -> Self {
        match v {
            0 => Base::Gap,
            1 => Base::A,
            2 => Base::C,
            3 => Base::G,
            4 => Base::N,
            5 => Base::U,
            _ => panic!("invalid base rank {v}"),
        }
    }

    /// Convert raw u8 byte representation.
    #[inline]
    pub const fn as_u8(self) -> u8 {
        self as u8
    }

    /// Get standard uppercase ASCII byte (A, G, C, U, N)
    #[inline]
    pub fn to_u8_upper(self) -> u8 {
        BASE_TO_UPPER[self.as_usize()]
    }

    /// Convert Base to lowercase ASCII byte (a, g, c, t, n, -)
    /// Used for display and output formatting.
    #[inline]
    pub fn to_byte(self) -> u8 {
        BASE_TO_BYTE[self.as_usize()]
    }

    /// Watson-Crick complement (A <-> U/T, G <-> C)
    #[inline(always)]
    pub const fn complement(self) -> Self {
        match self {
            Base::Gap => Base::Gap,
            Base::A => Base::U,
            Base::G => Base::C,
            Base::C => Base::G,
            Base::U => Base::A,
            Base::N => Base::N,
        }
    }

    /// True for the four matchable bases (A, C, G, U); false for N and Gap.
    #[inline(always)]
    pub const fn is_matchable(self) -> bool {
        matches!(self, Base::A | Base::C | Base::G | Base::U)
    }

    /// Classify pairing between two bases.
    ///
    /// Watson-Crick: A↔U, C↔G.  Wobble: G↔U.  Everything else: mismatch.
    #[inline(always)]
    pub const fn pair_type(self, other: Base) -> PairType {
        match (self, other) {
            (Base::A, Base::U) | (Base::U, Base::A) | (Base::C, Base::G) | (Base::G, Base::C) => {
                PairType::Canonical
            }
            (Base::G, Base::U) | (Base::U, Base::G) => PairType::Wobble,
            _ => PairType::Mismatch,
        }
    }
}

/// Chemistry of a base pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairType {
    /// Watson-Crick: A↔U, C↔G
    Canonical,
    /// G↔U wobble
    Wobble,
    /// Non-pairing
    Mismatch,
}

impl PairType {
    /// Whether this pair type counts as a match under the given wobble mode.
    #[inline(always)]
    pub const fn is_match(self, wobble: bool) -> bool {
        match self {
            PairType::Canonical => true,
            PairType::Wobble => wobble,
            PairType::Mismatch => false,
        }
    }
}

/// Number of nucleotide types (Gap, A, C, G, N, U)
pub const BASE_COUNT: usize = 6;

/// Constant for Gap index used in array indexing and DSM lookups.
pub const GAP: u8 = Base::Gap.as_u8();

/// Strand direction for search

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

impl From<char> for Strand {
    fn from(c: char) -> Self {
        match c {
            '+' => Strand::Forward,
            _ => Strand::Reverse,
        }
    }
}

impl From<Strand> for char {
    fn from(s: Strand) -> char {
        match s {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        }
    }
}

/// Represents thermodynamic energy, internally stored as an integer (scaling kcal/mol by 10,000).
/// All conversions from floating-point values expect inputs in kcal/mol.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct Energy(pub(crate) i32);
impl Energy {
    pub const SCALE: f64 = 10000.0;

    /// Most-negative threshold. `energy <= Energy::MIN` is ~never true, so used
    /// as a filter bound it records nothing — the integer-domain stand-in for
    /// a `-inf` threshold (which `from_kcal` rightly rejects as non-finite).
    pub const MIN: Energy = Energy(i32::MIN);

    /// Create Energy from a value in kcal/mol. Panics if the value is non-finite or overflows.
    /// Use this for trusted internal constants and tests.
    pub fn from_kcal(kcal: f64) -> Self {
        Self::try_from(kcal).expect("invalid energy value")
    }

    /// Convert the internal integer representation back to kcal/mol.
    pub fn to_kcal(self) -> f64 {
        f64::from(self.0) / Self::SCALE
    }
}

impl TryFrom<f64> for Energy {
    type Error = String;

    /// Create Energy from a value in kcal/mol.
    /// Returns an error if the value is non-finite or would overflow the internal representation.
    fn try_from(kcal: f64) -> Result<Self, Self::Error> {
        if !kcal.is_finite() {
            return Err(format!("non-finite energy: {kcal}"));
        }
        let score = kcal * Self::SCALE;
        if score < f64::from(i32::MIN) || score > f64::from(i32::MAX) {
            return Err(format!("energy overflow: {kcal} kcal/mol"));
        }
        Ok(Energy(score.round() as i32))
    }
}

impl From<Energy> for f64 {
    fn from(e: Energy) -> Self {
        e.to_kcal()
    }
}

impl PartialOrd for Energy {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Energy {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

impl std::str::FromStr for Energy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let kcal = s.parse::<f64>().map_err(|e| e.to_string())?;
        Self::try_from(kcal)
    }
}

#[cfg(test)]
mod tests {
    use super::Base;

    #[test]
    fn base_discriminants_are_canonical_internal_rank() {
        assert_eq!(Base::Gap.as_u8(), 0);
        assert_eq!(Base::A.as_u8(), 1);
        assert_eq!(Base::C.as_u8(), 2);
        assert_eq!(Base::G.as_u8(), 3);
        assert_eq!(Base::N.as_u8(), 4);
        assert_eq!(Base::U.as_u8(), 5);
    }
}

impl std::ops::Add for Energy {
    type Output = Self;

    #[track_caller]
    fn add(self, rhs: Self) -> Self::Output {
        Self(
            self.0
                .checked_add(rhs.0)
                .expect("Energy addition overflowed i32"),
        )
    }
}

impl std::ops::Sub for Energy {
    type Output = Self;

    #[track_caller]
    fn sub(self, rhs: Self) -> Self::Output {
        Self(
            self.0
                .checked_sub(rhs.0)
                .expect("Energy subtraction overflowed i32"),
        )
    }
}

impl std::ops::Mul<usize> for Energy {
    type Output = Self;

    #[track_caller]
    fn mul(self, rhs: usize) -> Self::Output {
        let rhs_i32 = i32::try_from(rhs).expect("Energy multiplier overflows i32");
        Self(
            self.0
                .checked_mul(rhs_i32)
                .expect("Energy multiplication overflowed i32"),
        )
    }
}

impl std::fmt::Display for Energy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.2}", self.to_kcal())
    }
}

/// Nucleic acid type for query/target strand identity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequenceType {
    Rna,
    Dna,
}

impl TryFrom<&str> for SequenceType {
    type Error = String;
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        match s {
            "rna" => Ok(SequenceType::Rna),
            "dna" => Ok(SequenceType::Dna),
            _ => Err(format!(
                "invalid sequence type '{s}', expected 'rna' or 'dna'"
            )),
        }
    }
}

/// Identifier for a bundled dinucleotide stacking model (DSM).
/// The set of valid IDs is determined by data/dsm/manifest.toml.
#[derive(Clone, Debug, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
pub struct DsmId(pub String);

impl std::fmt::Display for DsmId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<&str> for DsmId {
    fn from(s: &str) -> Self {
        Self(s.to_string())
    }
}
