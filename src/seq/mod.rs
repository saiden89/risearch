//! Sequence abstraction for clean base access.
//!
//! Normalized sequences are stored as Base enums (A/G/C/U/N/Gap).
//! This module centralizes normalization so the rest of the repository relies
//! on a single canonical representation.

pub(crate) mod normalize;
pub mod sequence;

pub use normalize::NormalizationStats;
pub use sequence::Sequence;
