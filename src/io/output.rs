use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use zstd::stream;

use crate::config::{OutputCompression, OutputFormat};
use crate::search::output::{write_results_to, write_results_with_format_to};
use crate::search::SearchHit;

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

pub fn write_results(hits: &[SearchHit], output: impl AsRef<Path>) -> Result<()> {
    let (mut writer, _compression) = open_output(output.as_ref(), None)?;
    write_results_to(hits, &mut writer)?;
    writer.flush().context("Failed to flush output")?;
    Ok(())
}

pub fn write_results_with_format(
    hits: &[SearchHit],
    output: impl AsRef<Path>,
    format: OutputFormat,
) -> Result<()> {
    let (mut writer, _compression) = open_output(output.as_ref(), None)?;
    write_results_with_format_to(hits, &mut writer, format)?;
    writer.flush().context("Failed to flush output")?;
    Ok(())
}
