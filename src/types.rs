//! Core domain types used across the codebase

use clap::ValueEnum;
use smallvec::SmallVec;

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

/// Static lookup table for byte-to-Base conversion (256 entries, O(1) access).
/// Avoids branch prediction misses from match chain in hot paths.
static BYTE_TO_BASE: [Base; 256] = {
    let mut table = [Base::N; 256];
    table[b'A' as usize] = Base::A;
    table[b'a' as usize] = Base::A;
    table[b'G' as usize] = Base::G;
    table[b'g' as usize] = Base::G;
    table[b'C' as usize] = Base::C;
    table[b'c' as usize] = Base::C;
    table[b'U' as usize] = Base::U;
    table[b'u' as usize] = Base::U;
    table[b'T' as usize] = Base::U;
    table[b't' as usize] = Base::U;
    table[b'-' as usize] = Base::Gap;
    table[b'.' as usize] = Base::Gap;
    table
};

/// Static lookup table for RNA reverse complement (byte → complemented byte).
/// Single lookup, no enum conversion, SIMD-vectorizable.
/// A↔U, G↔C, unknown→N, outputs uppercase.
pub static RC_RNA_TABLE: [u8; 256] = {
    let mut t = [b'N'; 256];
    t[b'A' as usize] = b'U';
    t[b'a' as usize] = b'U';
    t[b'U' as usize] = b'A';
    t[b'u' as usize] = b'A';
    t[b'T' as usize] = b'A';
    t[b't' as usize] = b'A';
    t[b'G' as usize] = b'C';
    t[b'g' as usize] = b'C';
    t[b'C' as usize] = b'G';
    t[b'c' as usize] = b'G';
    t[b'-' as usize] = b'-';
    t[b'.' as usize] = b'-';
    t
};

/// DNA reverse complement LUT: byte → complemented lowercase byte (T not U).
pub static RC_DNA_TABLE: [u8; 256] = {
    let mut t = [b'n'; 256];
    t[b'A' as usize] = b't';
    t[b'a' as usize] = b't';
    t[b'U' as usize] = b'a';
    t[b'u' as usize] = b'a';
    t[b'T' as usize] = b'a';
    t[b't' as usize] = b'a';
    t[b'G' as usize] = b'c';
    t[b'g' as usize] = b'c';
    t[b'C' as usize] = b'g';
    t[b'c' as usize] = b'g';
    t[b'-' as usize] = b'-';
    t[b'.' as usize] = b'-';
    t
};

// =============================================================================
// BASE CONVERSION LUTS - Constant-time lookups for Base enum
// =============================================================================

/// Base → uppercase ASCII byte
static BASE_TO_UPPER: [u8; 6] = [b'-', b'A', b'G', b'C', b'U', b'N'];

/// Base → complement Base (indexed by Base as usize)
static BASE_COMPLEMENT: [Base; 6] = [Base::Gap, Base::U, Base::C, Base::G, Base::A, Base::N];

/// Byte → complement byte (256-entry LUT for direct ASCII lookup)
/// A<->U/T, C<->G, N->N, others->N (all lowercase output)
pub static COMPLEMENT: [u8; 256] = {
    let mut lut = [b'n'; 256];
    lut[b'a' as usize] = b't';
    lut[b'A' as usize] = b't';
    lut[b't' as usize] = b'a';
    lut[b'T' as usize] = b'a';
    lut[b'u' as usize] = b'a';
    lut[b'U' as usize] = b'a';
    lut[b'c' as usize] = b'g';
    lut[b'C' as usize] = b'g';
    lut[b'g' as usize] = b'c';
    lut[b'G' as usize] = b'c';
    lut[b'n' as usize] = b'n';
    lut[b'N' as usize] = b'n';
    lut
};

/// Index → Base (for from_idx)
static IDX_TO_BASE: [Base; 6] = [Base::Gap, Base::A, Base::G, Base::C, Base::U, Base::N];

impl Base {
    /// Convert ASCII nucleotide byte to Base enum via lookup table.
    #[inline(always)]
    pub fn from_byte(b: u8) -> Self {
        BYTE_TO_BASE[b as usize]
    }

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

    /// Get char representation
    #[inline]
    pub fn as_char(self) -> char {
        self.to_u8_upper() as char
    }

    /// Watson-Crick complement (A <-> U/T, G <-> C)
    #[inline]
    pub fn complement(self) -> Self {
        BASE_COMPLEMENT[self as usize]
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

// =============================================================================
// SPAN - Lightweight region reference
// =============================================================================

/// Lightweight, Copy-able reference to a sequence region.
///
/// 16 bytes total: no lifetime, cache-optimal for batch processing.
/// Use `seq_id` to look up actual sequence data from an index.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct Span {
    /// Start position in sequence (0-based, inclusive)
    pub start: u32,
    /// End position in sequence (0-based, exclusive)
    pub end: u32,
    /// Index into sequence database (e.g., target index)
    pub seq_id: u32,
    /// Strand direction
    pub strand: Strand,
}

impl Span {
    /// Create a new Span
    #[inline]
    pub const fn new(seq_id: u32, start: u32, end: u32, strand: Strand) -> Self {
        Self {
            start,
            end,
            seq_id,
            strand,
        }
    }

    /// Length of the region
    #[inline]
    pub const fn len(&self) -> u32 {
        self.end.saturating_sub(self.start)
    }

    /// Check if span is empty
    #[inline]
    pub const fn is_empty(&self) -> bool {
        self.start >= self.end
    }
}

// =============================================================================
// PAIRING - Base pair classification for alignments
// =============================================================================

/// Classification of a query-target base pair in an alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pairing {
    Match(Base, Base),    // Watson-Crick pair (A-U, G-C)
    Wobble(Base, Base),   // G-U wobble pair
    Mismatch(Base, Base), // Non-complementary bases
    GapQuery(Base),       // Gap in Query, Base in Target
    GapTarget(Base),      // Gap in Target, Base in Query
}

impl Pairing {
    /// Construct from raw bytes (convenience for callsites that have bytes)
    #[inline]
    pub fn from_bytes(q_byte: u8, t_byte: u8) -> Self {
        Self::from_bases(Base::from_byte(q_byte), Base::from_byte(t_byte))
    }

    /// Construct from Base enums (preferred)
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

    /// Get the fingerprint character for this pairing
    #[inline]
    pub fn to_char(&self) -> char {
        match self {
            Pairing::Match(_, _) => 'P',
            Pairing::Wobble(_, _) => 'W',
            Pairing::Mismatch(_, _) => 'U',
            Pairing::GapQuery(_) => 'T',  // Gap in Query = Target Bulge
            Pairing::GapTarget(_) => 'Q', // Gap in Target = Query Bulge
        }
    }

    /// Get the target base character for alignment display
    #[inline]
    pub fn target_char(&self) -> char {
        match self {
            Pairing::Match(_, t)
            | Pairing::Wobble(_, t)
            | Pairing::Mismatch(_, t)
            | Pairing::GapQuery(t) => t.as_char().to_ascii_lowercase(),
            Pairing::GapTarget(_) => '-',
        }
    }

    /// Get the query base character for alignment display
    #[inline]
    pub fn query_char(&self) -> char {
        match self {
            Pairing::Match(q, _)
            | Pairing::Wobble(q, _)
            | Pairing::Mismatch(q, _)
            | Pairing::GapTarget(q) => q.as_char().to_ascii_lowercase(),
            Pairing::GapQuery(_) => '-',
        }
    }

    /// Get query base (if present)
    #[inline]
    pub fn query_base(&self) -> Option<Base> {
        match self {
            Pairing::Match(q, _)
            | Pairing::Wobble(q, _)
            | Pairing::Mismatch(q, _)
            | Pairing::GapTarget(q) => Some(*q),
            Pairing::GapQuery(_) => None,
        }
    }

    /// Get target base (if present)
    #[inline]
    pub fn target_base(&self) -> Option<Base> {
        match self {
            Pairing::Match(_, t)
            | Pairing::Wobble(_, t)
            | Pairing::Mismatch(_, t)
            | Pairing::GapQuery(t) => Some(*t),
            Pairing::GapTarget(_) => None,
        }
    }
}

// =============================================================================
// ALIGNMENT - Full query-target alignment structure
// =============================================================================

use std::ops::Range;

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

    /// Returns only the seed region steps
    pub fn seed(&self) -> &[Pairing] {
        &self.steps[self.seed_range.clone()]
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

        for (i, fp_char) in fp_chars.iter().enumerate() {
            let t_base = tgt_chars.get(i).copied().unwrap_or('-');
            let t = Base::from_byte(t_base as u8);

            // We don't have query bases from C output, use N as placeholder
            let pairing = match fp_char {
                'P' => Pairing::Match(Base::N, t),
                'W' => Pairing::Wobble(Base::N, t),
                'U' => Pairing::Mismatch(Base::N, t),
                'T' => Pairing::GapQuery(t),        // Gap in query
                'Q' => Pairing::GapTarget(Base::N), // Gap in target
                _ => Pairing::Mismatch(Base::N, t),
            };
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

/// Query with ID and sequence (borrowed for processing).
#[derive(Debug, Clone, Copy)]
pub struct Query<'a> {
    pub id: &'a str,
    pub seq: &'a [u8],
}

impl<'a> Query<'a> {
    pub fn new(id: &'a str, seq: &'a [u8]) -> Self {
        Self { id, seq }
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
#[derive(Clone, Copy, PartialEq, Debug, Default)]
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

impl PartialOrd for Energy {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Energy {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0
            .partial_cmp(&other.0)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

impl From<f64> for Energy {
    fn from(v: f64) -> Self {
        Energy(v)
    }
}
