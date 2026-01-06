//! Core domain types used across the codebase

use clap::ValueEnum;

/// Nucleotide/gap representation for DSM indexing and sequence operations
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Base {
    #[default]
    Gap = 0,
    A = 1,
    G = 2,
    C = 3,
    U = 4,
    N = 5,
}

impl Base {
    /// Convert ASCII nucleotide byte to Base enum
    #[inline]
    pub fn from_byte(b: u8) -> Self {
        match b {
            b'A' | b'a' => Base::A,
            b'G' | b'g' => Base::G,
            b'C' | b'c' => Base::C,
            b'U' | b'u' | b'T' | b't' => Base::U,
            b'-' | b'.' => Base::Gap,
            _ => Base::N,
        }
    }

    /// Convert to array index
    #[inline]
    pub const fn idx(self) -> usize {
        self as usize
    }

    /// Convert from usize index to Base
    #[inline]
    pub fn from_idx(i: usize) -> Self {
        match i {
            0 => Base::Gap,
            1 => Base::A,
            2 => Base::G,
            3 => Base::C,
            4 => Base::U,
            5 => Base::N,
            _ => panic!("Invalid Base index: {}", i),
        }
    }

    /// Get standard uppercase ASCII byte (A, G, C, U, N)
    #[inline]
    pub fn to_u8_upper(self) -> u8 {
        match self {
            Base::A => b'A',
            Base::G => b'G',
            Base::C => b'C',
            Base::U => b'U',
            Base::Gap => b'-',
            Base::N => b'N',
        }
    }

    /// Get char representation
    #[inline]
    pub fn as_char(self) -> char {
        self.to_u8_upper() as char
    }

    /// Watson-Crick complement (A <-> U/T, G <-> C)
    #[inline]
    pub fn complement(self) -> Self {
        match self {
            Base::A => Base::U,
            Base::U => Base::A,
            Base::G => Base::C,
            Base::C => Base::G,
            _ => self,
        }
    }

    /// Get fingerprint character for this pair (P=Paired, W=Wobble, U=Unpaired)
    #[inline]
    pub fn pairing_class(self, other: Base) -> char {
        match (self, other) {
            (Base::A, Base::U) | (Base::U, Base::A) | (Base::G, Base::C) | (Base::C, Base::G) => {
                'P'
            }
            (Base::G, Base::U) | (Base::U, Base::G) => 'W',
            _ => 'U',
        }
    }
}

/// Number of nucleotide types (Gap, A, G, C, U, N)
pub const BASE_COUNT: usize = 6;

/// Strand direction for search
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

/// Seed pairing mode (wobble vs strict)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "snake_case")]
pub enum SeedPairing {
    AllowWobble,
    Strict,
}

impl From<bool> for SeedPairing {
    fn from(wobble_arg: bool) -> Self {
        if wobble_arg {
            SeedPairing::Strict
        } else {
            SeedPairing::AllowWobble
        }
    }
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
// ID NEWTYPES
// =============================================================================

/// Query sequence identifier (newtype for type safety).
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub struct QueryId(pub String);

impl QueryId {
    /// Create from string.
    pub fn new(s: impl Into<String>) -> Self {
        QueryId(s.into())
    }

    /// Get the ID string.
    pub fn as_str(&self) -> &str {
        &self.0
    }

    /// Truncate at first whitespace (for C output compatibility).
    pub fn truncated(&self) -> &str {
        self.0.split_whitespace().next().unwrap_or(&self.0)
    }
}

impl std::fmt::Display for QueryId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<&str> for QueryId {
    fn from(s: &str) -> Self {
        QueryId(s.to_string())
    }
}

impl From<String> for QueryId {
    fn from(s: String) -> Self {
        QueryId(s)
    }
}

/// Target sequence identifier (newtype for type safety).
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub struct TargetId(pub String);

impl TargetId {
    /// Create from string.
    pub fn new(s: impl Into<String>) -> Self {
        TargetId(s.into())
    }

    /// Get the ID string.
    pub fn as_str(&self) -> &str {
        &self.0
    }

    /// Truncate at first whitespace (for C output compatibility).
    pub fn truncated(&self) -> &str {
        self.0.split_whitespace().next().unwrap_or(&self.0)
    }
}

impl std::fmt::Display for TargetId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<&str> for TargetId {
    fn from(s: &str) -> Self {
        TargetId(s.to_string())
    }
}

impl From<String> for TargetId {
    fn from(s: String) -> Self {
        TargetId(s)
    }
}

// =============================================================================
// ENERGY NEWTYPE
// =============================================================================

/// Energy value in kcal/mol (newtype for type safety and formatting).
#[derive(Clone, Copy, PartialEq, PartialOrd, Debug, Default)]
pub struct Energy(pub f64);

impl Energy {
    /// Create from f64.
    pub fn new(value: f64) -> Self {
        Energy(value)
    }

    /// Get the raw f64 value.
    pub fn as_f64(&self) -> f64 {
        self.0
    }

    /// Format as string with 2 decimal places (matches C output).
    pub fn format(&self) -> String {
        format!("{:.2}", self.0)
    }

    /// Parse from string.
    pub fn parse(s: &str) -> Option<Self> {
        s.parse::<f64>().ok().map(Energy)
    }
}

impl std::fmt::Display for Energy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.2}", self.0)
    }
}

impl Eq for Energy {}

impl Ord for Energy {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap_or(std::cmp::Ordering::Equal)
    }
}

impl From<f64> for Energy {
    fn from(v: f64) -> Self {
        Energy(v)
    }
}
