//! Sequence abstraction for clean base access.
//!
//! Normalized sequences are stored as Base enums (A/G/C/U/N/Gap).
//! This module centralizes normalization, reverse-complement, and RNA-formatting helpers so the
//! rest of the repository relies on a single canonical representation.

pub mod normalize;
pub mod sequence;
pub mod utils;

pub use sequence::Sequence;
pub use utils::{bases_to_rna_string, bytes_to_rna_string, push_bases_as_rna, push_bytes_as_rna};
