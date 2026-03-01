use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::config::{OutputConfig, OutputFormat};
use crate::registry::QueryRegistry;
use crate::search::SearchHit;
use crate::seq::SeqView;
use crate::types::Base;

use super::format::{append_hit_names_vec, OutputBuffers};
use super::{open_output, output_extension};

const CHUNK_SIZE_THRESHOLD: usize = 64 * 1024;

/// A batch of pre-formatted hit lines ready for writing.
pub struct OutputChunk {
    pub data: Vec<u8>,
    pub hits: usize,
    /// When multi-file output is active, identifies which query produced this chunk.
    pub query_idx: Option<u32>,
}

/// Formats `SearchHit`s into `OutputChunk`s, buffering to avoid excessive channel traffic.
pub struct HitFormatter {
    fmt_bufs: OutputBuffers,
    chunk_data: Vec<u8>,
    chunk_hits: usize,
    output_format: OutputFormat,
}

impl HitFormatter {
    pub fn new(output_format: OutputFormat) -> Self {
        Self {
            fmt_bufs: OutputBuffers::new(),
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
        query_name: &str,
        target_name: &str,
        query_seq: &[Base],
        t_fwd: &[Base],
        t_rc: &[Base],
        query_idx: Option<u32>,
    ) -> Option<OutputChunk> {
        append_hit_names_vec(
            &mut self.fmt_bufs,
            hit,
            self.output_format,
            &mut self.chunk_data,
            query_name,
            target_name,
            SeqView::from(query_seq),
            SeqView::from(t_fwd),
            SeqView::from(t_rc),
        );
        self.chunk_hits += 1;

        if self.chunk_data.len() >= CHUNK_SIZE_THRESHOLD {
            Some(self.take_chunk(query_idx))
        } else {
            None
        }
    }

    pub fn take_chunk(&mut self, query_idx: Option<u32>) -> OutputChunk {
        let data = std::mem::replace(
            &mut self.chunk_data,
            Vec::with_capacity(CHUNK_SIZE_THRESHOLD + 1024),
        );
        let hits = std::mem::replace(&mut self.chunk_hits, 0);
        OutputChunk {
            data,
            hits,
            query_idx,
        }
    }

    pub fn flush(&mut self, query_idx: Option<u32>) -> Option<OutputChunk> {
        if self.chunk_hits > 0 {
            Some(self.take_chunk(query_idx))
        } else {
            None
        }
    }
}

/// Orchestrates actual output writing, handling single-file and multi-file logic.
pub struct OutputWriter {
    config: OutputConfig,
    single_writer: Option<Box<dyn Write>>,
    multi_writers: HashMap<u32, Box<dyn Write>>,
    multi_paths: Vec<PathBuf>,
}

impl OutputWriter {
    pub fn new(config: &OutputConfig, output_path: &Path, queries: &QueryRegistry) -> Result<Self> {
        let config = config.clone();
        if config.multifile {
            std::fs::create_dir_all(output_path).with_context(|| {
                format!(
                    "Failed to create output directory {:?} for --multifile",
                    output_path
                )
            })?;
            let ext = output_extension(&config);
            let multi_paths = build_multifile_paths(queries, output_path, ext);
            Ok(Self {
                config,
                single_writer: None,
                multi_writers: HashMap::new(),
                multi_paths,
            })
        } else {
            let single_writer = open_output(Some(output_path), &config)?;
            Ok(Self {
                config,
                single_writer: Some(single_writer),
                multi_writers: HashMap::new(),
                multi_paths: Vec::new(),
            })
        }
    }

    pub fn write_chunk(&mut self, chunk: &OutputChunk) -> Result<()> {
        if let Some(w) = self.single_writer.as_mut() {
            w.write_all(&chunk.data)
                .context("Failed to write output chunk")?;
        } else {
            let qi = chunk
                .query_idx
                .expect("BUG: multifile mode active but OutputChunk has no query_idx");
            let file_path = self.multi_paths.get(qi as usize)
                .ok_or_else(|| anyhow::anyhow!("Invalid query index {} for multi-file output", qi))?;
            let config = &self.config;
            let writer = match self.multi_writers.get_mut(&qi) {
                Some(w) => w,
                None => {
                    let w = open_output(Some(file_path), config)?;
                    self.multi_writers.entry(qi).or_insert(w)
                }
            };
            writer
                .write_all(&chunk.data)
                .context("Failed to write output chunk to per-query file")?;
        }
        Ok(())
    }

    pub fn flush_all(&mut self) -> Result<()> {
        if let Some(w) = self.single_writer.as_mut() {
            w.flush().context("Failed to flush output")?;
        }
        for w in self.multi_writers.values_mut() {
            w.flush().context("Failed to flush per-query output")?;
        }
        Ok(())
    }
}

/// Replace characters that are unsafe in filenames and handle reserved names.
fn sanitize_filename(name: &str) -> String {
    let mut s: String = name.chars()
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
        "CON" | "PRN" | "AUX" | "NUL" | "COM1" | "COM2" | "COM3" | "COM4" | "COM5" | "COM6" | "COM7" | "COM8" | "COM9"
        | "LPT1" | "LPT2" | "LPT3" | "LPT4" | "LPT5" | "LPT6" | "LPT7" | "LPT8" | "LPT9" => {
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

fn build_multifile_paths(queries: &QueryRegistry, output_dir: &Path, ext: &str) -> Vec<PathBuf> {
    let mut used_stems: HashSet<String> = HashSet::with_capacity(queries.len());
    let mut out: Vec<PathBuf> = Vec::with_capacity(queries.len());

    for query_idx in 0..queries.len() {
        let name = queries.get_name(query_idx as u32);
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
