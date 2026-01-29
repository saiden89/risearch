use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use zstd::stream;

use crate::config::OutputCompression;

pub mod format;

#[derive(Debug, Clone, Copy)]
pub enum CompressionConfig {
    None,
    Gzip(Compression),
    Zstd(i32),
}

pub fn resolve_compression(
    codec: OutputCompression,
    level: Option<i32>,
) -> Result<CompressionConfig> {
    match codec {
        OutputCompression::None => {
            if level.is_some() {
                anyhow::bail!("--output-level requires compressed output");
            }
            Ok(CompressionConfig::None)
        }
        OutputCompression::Gzip => {
            let lvl = match level {
                None => Compression::default(),
                Some(v) => {
                    if !(0..=9).contains(&v) {
                        anyhow::bail!("gzip level must be in 0..=9 (got {})", v);
                    }
                    Compression::new(v as u32)
                }
            };
            Ok(CompressionConfig::Gzip(lvl))
        }
        OutputCompression::Zstd => {
            let lvl = level.unwrap_or(0);
            if !(-7..=22).contains(&lvl) {
                anyhow::bail!("zstd level must be in -7..=22 (got {})", lvl);
            }
            Ok(CompressionConfig::Zstd(lvl))
        }
    }
}

pub fn compress_bytes(config: CompressionConfig, input: Vec<u8>) -> Result<Vec<u8>> {
    match config {
        CompressionConfig::None => Ok(input),
        CompressionConfig::Gzip(level) => {
            let mut encoder = GzEncoder::new(Vec::new(), level);
            encoder.write_all(&input)?;
            let out = encoder.finish().context("gzip finish failed")?;
            Ok(out)
        }
        CompressionConfig::Zstd(level) => {
            let mut encoder = stream::write::Encoder::new(Vec::new(), level)
                .context("zstd encoder init failed")?;
            encoder.write_all(&input)?;
            let out = encoder.finish().context("zstd finish failed")?;
            Ok(out)
        }
    }
}

pub fn infer_compression(path: &Path) -> OutputCompression {
    let Some(ext) = path.extension().and_then(|e| e.to_str()) else {
        return OutputCompression::None;
    };
    match ext.to_ascii_lowercase().as_str() {
        "gz" | "gzip" => OutputCompression::Gzip,
        "zst" | "zstd" => OutputCompression::Zstd,
        _ => OutputCompression::None,
    }
}

pub fn open_output(
    path: &Path,
    override_compress: Option<OutputCompression>,
) -> Result<(BufWriter<Box<dyn Write>>, OutputCompression)> {
    let inner: Box<dyn Write> = if path == Path::new("-") {
        Box::new(std::io::stdout())
    } else {
        Box::new(std::fs::File::create(path).context("Failed to create output file")?)
    };

    let compression = override_compress.unwrap_or_else(|| {
        if path == Path::new("-") {
            OutputCompression::None
        } else {
            infer_compression(path)
        }
    });

    let writer = BufWriter::with_capacity(256 * 1024, inner);
    Ok((writer, compression))
}

/// Open output file with streaming compression applied.
/// Returns a boxed writer that compresses on-the-fly.
pub fn open_compressed_output(
    path: Option<impl AsRef<Path>>,
    compress: Option<OutputCompression>,
    level: Option<i32>,
) -> Result<Box<dyn Write>> {
    let path_ref = path.as_ref().map(|p| p.as_ref());
    let inner: Box<dyn Write> = match path_ref {
        Some(p) if p != Path::new("-") => {
            Box::new(std::fs::File::create(p).context("Failed to create output file")?)
        }
        _ => Box::new(std::io::stdout()),
    };

    let compression = compress.unwrap_or_else(|| {
        path_ref
            .map(infer_compression)
            .unwrap_or(OutputCompression::None)
    });

    let config = resolve_compression(compression, level)?;

    match config {
        CompressionConfig::None => Ok(Box::new(BufWriter::with_capacity(256 * 1024, inner))),
        CompressionConfig::Gzip(lvl) => Ok(Box::new(BufWriter::with_capacity(
            256 * 1024,
            GzEncoder::new(inner, lvl),
        ))),
        CompressionConfig::Zstd(lvl) => {
            let encoder = stream::write::Encoder::new(inner, lvl)
                .context("zstd encoder init failed")?
                .auto_finish();
            Ok(Box::new(BufWriter::with_capacity(256 * 1024, encoder)))
        }
    }
}
