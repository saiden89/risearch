use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use risearch::args::OutputCompression;


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
