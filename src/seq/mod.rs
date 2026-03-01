//! Sequence abstraction for clean base access.
//!
//! Normalized sequences are stored as Base enums (A/G/C/U/N/Gap).
//! This module centralizes normalization, reverse-complement, and RNA-formatting helpers so the
//! rest of the repository relies on a single canonical representation.

pub mod normalize;
pub mod sequence;
pub mod utils;
pub mod view;

pub use sequence::Sequence;
pub use utils::{bases_to_rna_string, bytes_to_rna_string, push_bases_as_rna, push_bytes_as_rna};
pub use view::SeqView;

// ALIGNED SEQ - Sentinel-padded buffer for safe pointer arithmetic
// =============================================================================

/// Sentinel padding sizes
const LOOKBEHIND_PAD: usize = 16; // Safe offset(-1) for 16 positions
const SIMD_PAD: usize = 64; // AVX-512 over-read safety

/// A sequence buffer with sentinel padding for safe lookbehind and SIMD over-read.
///
/// Guarantees:
/// - `ptr.offset(-1)` is valid for the first LOOKBEHIND_PAD bytes
/// - Reading past the end by up to SIMD_PAD bytes is safe (returns 'N')
///
/// Use this when you need pointer arithmetic in hot loops.
#[derive(Debug, Clone)]
pub struct AlignedSeq {
    /// Allocation: [sentinel(16)] + [data] + [sentinel(64)]
    buffer: Vec<u8>,
    /// Start index of actual data within buffer
    data_start: usize,
    /// Length of actual sequence data
    data_len: usize,
}

impl AlignedSeq {
    /// Create a new AlignedSeq from raw sequence bytes.
    pub fn new(data: &[u8]) -> Self {
        let mut buffer = Vec::with_capacity(LOOKBEHIND_PAD + data.len() + SIMD_PAD);

        // Leading sentinels (N for unknown base)
        buffer.extend(std::iter::repeat_n(b'N', LOOKBEHIND_PAD));
        let data_start = buffer.len();

        // Actual sequence
        buffer.extend_from_slice(data);
        let data_len = data.len();

        // Trailing sentinels
        buffer.extend(std::iter::repeat_n(b'N', SIMD_PAD));

        Self {
            buffer,
            data_start,
            data_len,
        }
    }

    /// Get a pointer to the start of actual data.
    /// Safe to offset(-1) up to LOOKBEHIND_PAD times.
    #[inline]
    pub fn as_ptr(&self) -> *const u8 {
        unsafe { self.buffer.as_ptr().add(self.data_start) }
    }

    /// Get a slice of the actual data (no padding).
    #[inline]
    pub fn as_slice(&self) -> &[u8] {
        &self.buffer[self.data_start..self.data_start + self.data_len]
    }

    /// Length of actual sequence data.
    #[inline]
    pub fn len(&self) -> usize {
        self.data_len
    }

    /// Check if empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data_len == 0
    }
}
