//! Core domain types used across the codebase

use clap::ValueEnum;
use serde::{Deserialize, Serialize};

/// Nucleotide/gap representation for DSM indexing and sequence operations
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
static BASE_TO_UPPER: [u8; 6] = [b'-', b'A', b'C', b'G', b'N', b'U'];

/// Base → lowercase ASCII byte (for to_byte())
static BASE_TO_BYTE: [u8; 6] = [b'-', b'a', b'c', b'g', b'n', b'u'];

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
    pub const fn idx(self) -> usize {
        self as usize
    }

    /// Convert from usize index to Base.
    ///
    /// # Safety
    /// Caller must ensure `i < 6`.
    #[inline(always)]
    pub unsafe fn from_idx(i: usize) -> Self {
        debug_assert!(i < 6, "Invalid Base index: {}", i);
        std::mem::transmute(i as u8)
    }

    /// Get standard uppercase ASCII byte (A, G, C, U, N)
    #[inline]
    pub fn to_u8_upper(self) -> u8 {
        BASE_TO_UPPER[self as usize]
    }

    /// Convert Base to lowercase ASCII byte (a, g, c, t, n, -)
    /// Used for display and output formatting.
    #[inline]
    pub fn to_byte(self) -> u8 {
        BASE_TO_BYTE[self as usize]
    }

    /// Get char representation
    #[inline]
    pub fn as_char(self) -> char {
        self.to_u8_upper() as char
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

    /// Classify pairing between two bases (in complement-transformed target space).
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

/// Classification of a base pair in complement-transformed target space.
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

/// Number of nucleotide types (Gap, A, G, C, U, N)
pub const BASE_COUNT: usize = 6;

/// Constant for Gap index used in array indexing and DSM lookups.
pub const GAP: u8 = Base::Gap as u8;

/// Strand direction for search

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
}

/// Seed pairing mode (wobble vs strict)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "snake_case")]
pub enum SeedPairingMode {
    AllowWobble,
    Strict,
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

#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct Energy(pub(crate) i32);

impl Energy {
    pub const SCALE: f64 = 10000.0;

    /// Create Energy from kcal/mol. Panics if value is non-finite or overflows i32.
    pub fn from_kcal(kcal: f64) -> Self {
        Self::try_from_kcal(kcal).expect("invalid energy value")
    }

    /// Safely create Energy from kcal/mol.
    pub fn try_from_kcal(kcal: f64) -> Result<Self, String> {
        if !kcal.is_finite() {
            return Err(format!("non-finite energy: {kcal}"));
        }
        let score = kcal * Self::SCALE;
        if score < i32::MIN as f64 || score > i32::MAX as f64 {
            return Err(format!("energy overflow: {kcal} kcal/mol"));
        }
        Ok(Energy(score.round() as i32))
    }

    pub fn to_kcal(self) -> f64 {
        f64::from(self.0) / Self::SCALE
    }
}

impl From<f64> for Energy {
    fn from(v: f64) -> Self {
        Energy::from_kcal(v)
    }
}

impl From<Energy> for f64 {
    fn from(e: Energy) -> Self {
        e.to_kcal()
    }
}

impl PartialOrd for Energy {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl std::str::FromStr for Energy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let kcal = s.parse::<f64>().map_err(|e| e.to_string())?;
        Self::try_from_kcal(kcal)
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
