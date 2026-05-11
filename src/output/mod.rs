//! Output formatting and compression.

pub mod format;
pub use format::format_hit_into;

pub mod writer;
pub use writer::{HitFormatter, OutputChunk, OutputWriter};
