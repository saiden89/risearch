use std::collections::HashSet;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use zstd::stream;

use crate::config::{OutputCompression, OutputConfig, OutputFormat};
use crate::registry::QueryRegistry;
use crate::search::SearchHit;

use super::format::format_hit_into;

const CHUNK_SIZE_THRESHOLD: usize = 64 * 1024;

/// A batch of pre-formatted hit lines ready for writing.
pub struct OutputChunk {
    pub data: Vec<u8>,
    pub hits: usize,
}

/// Formats `SearchHit`s into `OutputChunk`s, buffering to avoid excessive channel traffic.
pub struct HitFormatter {
    itoa: itoa::Buffer,
    chunk_data: Vec<u8>,
    chunk_hits: usize,
    output_format: OutputFormat,
}

impl HitFormatter {
    pub fn new(output_format: OutputFormat) -> Self {
        Self {
            itoa: itoa::Buffer::new(),
            chunk_data: Vec::with_capacity(CHUNK_SIZE_THRESHOLD + 1024),
            chunk_hits: 0,
            output_format,
        }
    }

    #[inline]
    #[allow(clippy::too_many_arguments)]
    pub fn add_hit(
        &mut self,
        hit: &SearchHit,
        q_name: &str,
        q_seq: &[crate::types::Base],
        t_name: &str,
        t_fwd: &[crate::types::Base],
        t_rc: &[crate::types::Base],
    ) -> Option<OutputChunk> {
        format_hit_into(
            &mut self.chunk_data,
            &mut self.itoa,
            hit,
            q_name,
            q_seq,
            t_name,
            t_fwd,
            t_rc,
            self.output_format,
        );
        self.chunk_hits += 1;

        if self.chunk_data.len() >= CHUNK_SIZE_THRESHOLD {
            Some(self.take_chunk())
        } else {
            None
        }
    }

    pub fn take_chunk(&mut self) -> OutputChunk {
        let data = std::mem::replace(
            &mut self.chunk_data,
            Vec::with_capacity(CHUNK_SIZE_THRESHOLD + 1024),
        );
        let hits = std::mem::replace(&mut self.chunk_hits, 0);
        OutputChunk { data, hits }
    }

    pub fn flush(&mut self) -> Option<OutputChunk> {
        if self.chunk_hits > 0 {
            Some(self.take_chunk())
        } else {
            None
        }
    }
}

/// Owns the full lifecycle of a single output destination: path resolution,
/// compression wrapping, buffering, writing, and flushing.
///
/// Pass `-` as `output_path` to write to stdout.
pub struct OutputWriter {
    writer: Box<dyn Write + Send>,
}

impl OutputWriter {
    pub fn new(config: &OutputConfig, output_path: &Path) -> Result<Self> {
        let inner: Box<dyn Write + Send> = if output_path == Path::new("-") {
            Box::new(std::io::stdout())
        } else {
            Box::new(
                std::fs::File::create(output_path)
                    .with_context(|| format!("Failed to create output file {:?}", output_path))?,
            )
        };

        let writer: Box<dyn Write + Send> = match config.compress {
            OutputCompression::None => Box::new(BufWriter::with_capacity(256 * 1024, inner)),
            OutputCompression::Gzip(level) => Box::new(BufWriter::with_capacity(
                256 * 1024,
                GzEncoder::new(inner, Compression::new(level as u32)),
            )),
            OutputCompression::Zstd(level) => {
                let encoder = stream::write::Encoder::new(inner, level)
                    .context("zstd encoder init failed")?
                    .auto_finish();
                Box::new(BufWriter::with_capacity(256 * 1024, encoder))
            }
        };

        Ok(Self { writer })
    }

    pub fn write_chunk(&mut self, chunk: &OutputChunk) -> Result<()> {
        self.writer
            .write_all(&chunk.data)
            .context("Failed to write output chunk")
    }

    pub fn flush_all(&mut self) -> Result<()> {
        self.writer.flush().context("Failed to flush output")
    }
}

/// Replace characters that are unsafe in filenames and handle reserved names.
fn sanitize_filename(name: &str) -> String {
    let mut s: String = name
        .chars()
        .map(|c| match c {
            '/' | '\\' | '\0' | ':' | '*' | '?' | '"' | '<' | '>' | '|' => '_',
            _ => c,
        })
        .collect();

    if s.is_empty() {
        return "query".to_string();
    }

    // Handle Windows reserved names (case-insensitive)
    let upper = s.to_uppercase();
    match upper.as_str() {
        "CON" | "PRN" | "AUX" | "NUL" | "COM1" | "COM2" | "COM3" | "COM4" | "COM5" | "COM6"
        | "COM7" | "COM8" | "COM9" | "LPT1" | "LPT2" | "LPT3" | "LPT4" | "LPT5" | "LPT6"
        | "LPT7" | "LPT8" | "LPT9" => {
            s.push('_');
        }
        _ => {}
    }
    s
}

fn unique_filename_stem(stem: &str, used: &mut HashSet<String>) -> String {
    if used.insert(stem.to_string()) {
        return stem.to_string();
    }

    let mut suffix = 1usize;
    loop {
        let candidate = format!("{stem}_{suffix}");
        if used.insert(candidate.clone()) {
            return candidate;
        }
        suffix += 1;
    }
}

pub(crate) fn build_multifile_paths(
    queries: &QueryRegistry,
    output_dir: &Path,
    ext: &str,
) -> Vec<PathBuf> {
    let mut used_stems: HashSet<String> = HashSet::with_capacity(queries.len());
    let mut out: Vec<PathBuf> = Vec::with_capacity(queries.len());

    for query_idx in 0..queries.len() {
        let name = queries.get_name(query_idx);
        let stem_raw = sanitize_filename(name);
        let stem = unique_filename_stem(&stem_raw, &mut used_stems);
        out.push(output_dir.join(format!("{stem}{ext}")));
    }

    out
}

#[cfg(test)]
mod tests {
    use super::{sanitize_filename, unique_filename_stem};
    use std::collections::HashSet;

    #[test]
    fn sanitize_replaces_unsafe_chars() {
        assert_eq!(sanitize_filename("a/b:c*?"), "a_b_c__");
    }

    #[test]
    fn unique_stem_avoids_collisions() {
        let mut used = HashSet::new();
        assert_eq!(unique_filename_stem("a_b", &mut used), "a_b");
        assert_eq!(unique_filename_stem("a_b", &mut used), "a_b_1");
        assert_eq!(unique_filename_stem("a_b", &mut used), "a_b_2");
    }
}
