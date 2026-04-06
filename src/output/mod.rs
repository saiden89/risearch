//! Output formatting and compression.

use crate::config::{OutputCompression, OutputConfig};

pub mod format;
pub use format::format_hit_into;

pub mod writer;
pub use writer::{HitFormatter, OutputChunk, OutputWriter};

/// Return a suitable file extension (including leading dot) for the given output
/// configuration, e.g. `".tsv"`, `".tsv.gz"`, `".tsv.zst"`.
pub fn output_extension(cfg: &OutputConfig) -> &'static str {
    match cfg.compress {
        OutputCompression::None => ".tsv",
        OutputCompression::Gzip(_) => ".tsv.gz",
        OutputCompression::Zstd(_) => ".tsv.zst",
    }
}
