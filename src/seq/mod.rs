//! Sequence abstraction for clean base access.
//!
//! Normalized sequences are stored as Base enums (A/G/C/U/N/Gap).
//! This module centralizes normalization and reverse-complement so the
//! rest of the repository relies on a single canonical representation.

pub mod normalize;
pub mod sequence;

pub use sequence::Sequence;
