//! Alignment representation for RNA interactions.
//!
//! This module contains the `Pairing` and `Alignment` types that represent
//! the result of RNA-RNA interaction prediction. These types are constructed
//! from DP traceback operations and used for output formatting.

use smallvec::SmallVec;
use std::ops::Range;

use crate::types::Base;

// =============================================================================
// PAIRING - Single position in an alignment
// =============================================================================

/// Classification of a query-target base pair in an alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pairing {
    Match(Base, Base),    // Watson-Crick pair (A-U, G-C)
    Wobble(Base, Base),   // G-U wobble pair
    Mismatch(Base, Base), // Non-complementary bases
    TargetBulge(Base),    // Target has extra base (gap in query alignment)
    QueryBulge(Base),     // Query has extra base (gap in target alignment)
}

impl Pairing {
    /// Construct from Base enums, automatically classifying the pair type
    #[inline]
    pub fn from_bases(q: Base, t: Base) -> Self {
        match (q, t) {
            (Base::G, Base::C) | (Base::C, Base::G) | (Base::A, Base::U) | (Base::U, Base::A) => {
                Pairing::Match(q, t)
            }
            (Base::G, Base::U) | (Base::U, Base::G) => Pairing::Wobble(q, t),
            _ => Pairing::Mismatch(q, t),
        }
    }

    /// Parse from fingerprint character (for C output compatibility).
    /// Returns None for unrecognized characters.
    #[inline]
    pub fn from_fingerprint_char(c: char, target_base: Base) -> Option<Self> {
        match c {
            'P' => Some(Pairing::Match(Base::N, target_base)),
            'W' => Some(Pairing::Wobble(Base::N, target_base)),
            'U' => Some(Pairing::Mismatch(Base::N, target_base)),
            'T' => Some(Pairing::TargetBulge(target_base)),
            'Q' => Some(Pairing::QueryBulge(Base::N)),
            _ => None,
        }
    }

    /// Get the fingerprint character for this pairing
    #[inline]
    pub fn to_char(&self) -> char {
        match self {
            Pairing::Match(_, _) => 'P',
            Pairing::Wobble(_, _) => 'W',
            Pairing::Mismatch(_, _) => 'U',
            Pairing::TargetBulge(_) => 'T',
            Pairing::QueryBulge(_) => 'Q',
        }
    }

    /// Get the query base character for alignment display ('-' for gaps)
    #[inline]
    pub fn query_char(&self) -> char {
        match self {
            Pairing::Match(q, _)
            | Pairing::Wobble(q, _)
            | Pairing::Mismatch(q, _)
            | Pairing::QueryBulge(q) => q.as_char().to_ascii_lowercase(),
            Pairing::TargetBulge(_) => '-',
        }
    }

    /// Get the target base character for alignment display ('-' for gaps)
    #[inline]
    pub fn target_char(&self) -> char {
        match self {
            Pairing::Match(_, t)
            | Pairing::Wobble(_, t)
            | Pairing::Mismatch(_, t)
            | Pairing::TargetBulge(t) => t.as_char().to_ascii_lowercase(),
            Pairing::QueryBulge(_) => '-',
        }
    }
}

// =============================================================================
// ALIGNMENT - Full query-target alignment structure
// =============================================================================

/// Represents a full biological alignment between Query and Target.
/// Stores the sequence of interactions and metadata about the seed location.
/// Uses SmallVec to avoid heap allocation for typical alignments (<128 steps).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    /// The complete sequence of pairing steps (5' -> 3' of Query).
    /// SmallVec keeps small alignments on the stack (128 Pairings ≈ 384 bytes).
    steps: SmallVec<[Pairing; 128]>,

    /// The range of indices in `steps` that corresponds to the initial Seed match.
    /// This allows easy extraction of the "core" interaction vs extensions.
    seed_range: Range<usize>,
}

impl Alignment {
    /// Constructor from the three phases of extension.
    /// This fits naturally into `extend_seed` which generates these 3 parts.
    /// Accepts slices to avoid forcing callers to allocate Vecs.
    pub fn new(left: &[Pairing], seed: &[Pairing], right: &[Pairing]) -> Self {
        let left_len = left.len();
        let seed_len = seed.len();

        let mut steps = SmallVec::with_capacity(left_len + seed_len + right.len());
        steps.extend_from_slice(left);
        steps.extend_from_slice(seed);
        steps.extend_from_slice(right);

        Self {
            steps,
            seed_range: left_len..(left_len + seed_len),
        }
    }

    /// Returns the full alignment steps
    pub fn steps(&self) -> &[Pairing] {
        &self.steps
    }

    /// Write query sequence (RNA-normalized) into a buffer.
    pub fn write_query_seq(&self, buf: &mut Vec<u8>) {
        for p in self.steps.iter() {
            buf.push(p.query_char() as u8);
        }
    }

    /// Write target sequence into a buffer.
    pub fn write_target_seq(&self, buf: &mut Vec<u8>) {
        for p in self.steps.iter() {
            buf.push(p.target_char() as u8);
        }
    }

    /// Write alignment line (|, :, or space) into a buffer.
    pub fn write_alignment_line(&self, buf: &mut Vec<u8>) {
        for p in self.steps.iter() {
            let c = match p {
                Pairing::Match(_, _) => b'|',
                Pairing::Wobble(_, _) => b':',
                _ => b' ',
            };
            buf.push(c);
        }
    }

    /// Write pairing fingerprint (P/W/U/...) into a buffer.
    pub fn write_pairing_string(&self, buf: &mut Vec<u8>) {
        for p in self.steps.iter() {
            buf.push(p.to_char() as u8);
        }
    }

    /// Returns only the seed region steps
    pub fn seed(&self) -> &[Pairing] {
        &self.steps[self.seed_range.start..self.seed_range.end]
    }

    /// Returns the 5' extension (Left of seed)
    pub fn left_extension(&self) -> &[Pairing] {
        &self.steps[..self.seed_range.start]
    }

    /// Returns the 3' extension (Right of seed)
    pub fn right_extension(&self) -> &[Pairing] {
        &self.steps[self.seed_range.end..]
    }

    /// Generates the interaction string (e.g. "PPPWUUU")
    pub fn fingerprint(&self) -> String {
        self.steps.iter().map(|p| p.to_char()).collect()
    }

    /// Compute a fast hash of the fingerprint for O(1) equality pre-check.
    /// Uses FNV-1a hash for speed (no allocation, just iterates over Pairings).
    #[inline]
    pub fn fingerprint_hash(&self) -> u64 {
        // FNV-1a hash constants for 64-bit
        const FNV_OFFSET: u64 = 0xcbf29ce484222325;
        const FNV_PRIME: u64 = 0x100000001b3;

        self.steps.iter().fold(FNV_OFFSET, |hash, p| {
            (hash ^ (p.to_char() as u64)).wrapping_mul(FNV_PRIME)
        })
    }

    /// Generates the target sequence string (e.g. "accu--cg")
    pub fn target_sequence(&self) -> String {
        self.steps.iter().map(|p| p.target_char()).collect()
    }

    /// Generates the query sequence string (e.g. "gc--uuca")
    pub fn query_sequence(&self) -> String {
        self.steps.iter().map(|p| p.query_char()).collect()
    }

    /// Generates the alignment string (e.g. "|| :  ")
    pub fn alignment_string(&self) -> String {
        self.steps
            .iter()
            .map(|p| match p {
                Pairing::Match(_, _) => '|',
                Pairing::Wobble(_, _) => ':',
                _ => ' ',
            })
            .collect()
    }

    /// Create from C output (fingerprint + target sequence + seed markers).
    /// C output uses 'y' and 'x' markers to delimit seed region.
    pub fn from_c_output(
        fingerprint: &str,
        target_seq: &str,
        seed_start: Option<usize>,
        seed_end: Option<usize>,
    ) -> Self {
        let fp_chars: Vec<char> = fingerprint.chars().collect();
        let tgt_chars: Vec<char> = target_seq.chars().collect();

        let mut steps = SmallVec::with_capacity(fp_chars.len());

        for (i, &fp_char) in fp_chars.iter().enumerate() {
            let t_base = tgt_chars.get(i).copied().unwrap_or('-');
            let t = Base::from_byte(t_base as u8);

            let pairing =
                Pairing::from_fingerprint_char(fp_char, t).unwrap_or(Pairing::Mismatch(Base::N, t));
            steps.push(pairing);
        }

        let seed_range = match (seed_start, seed_end) {
            (Some(s), Some(e)) => s..e,
            _ => 0..steps.len(), // If no markers, treat entire thing as seed
        };

        Self { steps, seed_range }
    }

    /// Get seed start index
    pub fn seed_start(&self) -> Option<usize> {
        if self.seed_range.start < self.steps.len() {
            Some(self.seed_range.start)
        } else {
            None
        }
    }

    /// Get seed end index
    pub fn seed_end(&self) -> Option<usize> {
        if self.seed_range.end <= self.steps.len() {
            Some(self.seed_range.end)
        } else {
            None
        }
    }
}
