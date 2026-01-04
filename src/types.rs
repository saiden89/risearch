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
            1 => Base::A,
            2 => Base::G,
            3 => Base::C,
            4 => Base::U,
            0 => Base::Gap,
            _ => Base::N,
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
