//! Output formatting and compression.

use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use zstd::stream;

use crate::config::{OutputCompression, OutputConfig};

pub mod format;
pub use format::{append_hit_minimal_names_vec, write_hit, write_hit_names, OutputBuffers};

/// Open output writer from parsed search output config.
pub fn open_output(path: Option<impl AsRef<Path>>, cfg: &OutputConfig) -> Result<Box<dyn Write>> {
    let path_ref = path.as_ref().map(|p| p.as_ref());

    // Create underlying writer
    let inner: Box<dyn Write> = match path_ref {
        Some(p) if p != Path::new("-") => {
            Box::new(std::fs::File::create(p).context("Failed to create output file")?)
        }
        _ => Box::new(std::io::stdout()),
    };

    let compression = cfg
        .compress
        .unwrap_or_else(|| path_ref.map(infer_compression).unwrap_or_default());

    // Wrap with compression encoder if needed
    match compression {
        OutputCompression::None => {
            if cfg.level.is_some() {
                bail!("--output-level requires compressed output");
            }
            Ok(Box::new(BufWriter::with_capacity(256 * 1024, inner)))
        }
        OutputCompression::Gzip => {
            let level = match cfg.level {
                None => Compression::default(),
                Some(v) if (0..=9).contains(&v) => Compression::new(v as u32),
                Some(v) => bail!("gzip level must be 0-9 (got {})", v),
            };
            Ok(Box::new(BufWriter::with_capacity(
                256 * 1024,
                GzEncoder::new(inner, level),
            )))
        }
        OutputCompression::Zstd => {
            let level = cfg.level.unwrap_or(0);
            if !(-7..=22).contains(&level) {
                bail!("zstd level must be -7..22 (got {})", level);
            }
            let encoder = stream::write::Encoder::new(inner, level)
                .context("zstd encoder init failed")?
                .auto_finish();
            Ok(Box::new(BufWriter::with_capacity(256 * 1024, encoder)))
        }
    }
}

fn infer_compression(path: &Path) -> OutputCompression {
    match path.extension().and_then(|e| e.to_str()) {
        Some(ext) => match ext.to_ascii_lowercase().as_str() {
            "gz" | "gzip" => OutputCompression::Gzip,
            "zst" | "zstd" => OutputCompression::Zstd,
            _ => OutputCompression::None,
        },
        None => OutputCompression::None,
    }
}
