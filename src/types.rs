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

// =============================================================================
// BASE CONVERSION LUTS - Constant-time lookups for Base enum
// =============================================================================

/// Base → uppercase ASCII byte
static BASE_TO_UPPER: [u8; 6] = [b'-', b'A', b'C', b'G', b'N', b'U'];

/// Base → lowercase ASCII byte (for to_byte())
static BASE_TO_BYTE: [u8; 6] = [b'-', b'a', b'c', b'g', b'n', b'u'];

/// Index → Base (for from_idx)
static IDX_TO_BASE: [Base; 6] = [Base::Gap, Base::A, Base::C, Base::G, Base::N, Base::U];

impl Base {
    /// Convert to array index
    #[inline]
    pub const fn idx(self) -> usize {
        self as usize
    }

    /// Convert from usize index to Base
    #[inline]
    pub fn from_idx(i: usize) -> Self {
        if i < 6 {
            IDX_TO_BASE[i]
        } else {
            panic!("Invalid Base index: {}", i)
        }
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
pub const GAP: usize = Base::Gap as usize;

// =============================================================================
// IDENTIFIERS / UNITS
// =============================================================================

/// Index into query registry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct QueryId(pub u32);

/// Index into target registry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct TargetId(pub u32);

/// Seed length (always positive by construction).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SeedLen(u16);

impl SeedLen {
    pub fn new(len: usize) -> Option<Self> {
        let len_u16 = u16::try_from(len).ok()?;
        if len_u16 == 0 {
            return None;
        }
        Some(Self(len_u16))
    }

    pub const fn get(self) -> usize {
        self.0 as usize
    }
}

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

// =============================================================================
// STRAND - Display and conversion impls
// =============================================================================

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

// =============================================================================
// ENERGY NEWTYPE
// =============================================================================

/// Energy value in kcal/mol.
///
/// Thin wrapper for type safety only.
#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub struct Energy(pub f64);

impl Energy {
    /// Create from raw kcal/mol value.
    #[inline]
    pub fn new(value: f64) -> Self {
        Energy(value)
    }

    /// Get the raw kcal/mol value.
    #[inline]
    pub fn as_f64(&self) -> f64 {
        self.0
    }
}

impl From<f64> for Energy {
    fn from(v: f64) -> Self {
        Energy(v)
    }
}

/// Convert DSM raw integer units → kcal/mol. 559 = terminal penalty offset.
impl From<i32> for Energy {
    fn from(raw: i32) -> Self {
        Energy((raw as f64 - 559.0) / -100.0)
    }
}
