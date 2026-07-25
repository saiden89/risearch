use std::collections::HashSet;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use zstd::stream;

use crate::config::OutputCompression;

/// Where a run's bytes go: one destination for everything, or one per key.
///
/// Routing is on an opaque `key`; this layer never learns what a key means. Both
/// variants open their files on first write, so constructing an `OutputWriter`
/// touches nothing — a caller whose config is later rejected leaves an existing
/// output file untruncated and creates no directory.
///
/// gzip and zstd write their trailer when the underlying stream drops rather than
/// on flush, so drop this only after the last write.
pub enum OutputWriter {
    /// Every key writes to the same stream, serialized. `-` means stdout.
    Single {
        path: PathBuf,
        stream: Mutex<Option<Box<dyn Write + Send>>>,
        compress: OutputCompression,
    },
    /// One file per key, opened and closed inside each write. Paths are supplied
    /// up front so naming never depends on the order keys arrive in.
    PerKey {
        dir: PathBuf,
        paths: Vec<PathBuf>,
        compress: OutputCompression,
    },
}

impl OutputWriter {
    /// One stream for every key, at `path` (`-` for stdout).
    pub fn single(path: &Path, compress: OutputCompression) -> Self {
        Self::Single {
            path: path.to_path_buf(),
            stream: Mutex::new(None),
            compress,
        }
    }

    /// One file per key. `paths` is indexed by key, so it must cover every key
    /// the caller will write.
    pub fn per_key(dir: &Path, paths: Vec<PathBuf>, compress: OutputCompression) -> Self {
        Self::PerKey {
            dir: dir.to_path_buf(),
            paths,
            compress,
        }
    }

    pub fn write_block(&self, key: usize, block: &[u8]) -> Result<()> {
        match self {
            Self::Single {
                path,
                stream,
                compress,
            } => {
                let mut slot = stream.lock().unwrap();
                open_stream(path, *compress, &mut slot)?
                    .write_all(block)
                    .context("Failed to write output")
            }
            Self::PerKey {
                dir,
                paths,
                compress,
            } => {
                ensure_dir(dir)?;
                let mut stream = open_path(&paths[key], *compress)?;
                stream.write_all(block).context("Failed to write output")?;
                stream.flush().context("Failed to flush output")
            }
        }
    }

    /// Finalize. `Single` opens its stream even if nothing was ever written, which
    /// is what truncates a stale output file after a run with no hits; `PerKey`
    /// flushed each file as it closed it, but still owes the caller the directory.
    pub fn finish(&self) -> Result<()> {
        match self {
            Self::Single {
                path,
                stream,
                compress,
            } => {
                let mut slot = stream.lock().unwrap();
                open_stream(path, *compress, &mut slot)?
                    .flush()
                    .context("Failed to flush output")
            }
            Self::PerKey { dir, .. } => ensure_dir(dir),
        }
    }
}

/// Opens `path` on first use. Idempotent by design: reopening would `File::create`
/// over bytes this run already wrote.
fn open_stream<'s>(
    path: &Path,
    compress: OutputCompression,
    slot: &'s mut Option<Box<dyn Write + Send>>,
) -> Result<&'s mut Box<dyn Write + Send>> {
    if slot.is_none() {
        *slot = Some(open_path(path, compress)?);
    }
    Ok(slot.as_mut().expect("just populated"))
}

fn open_path(path: &Path, compress: OutputCompression) -> Result<Box<dyn Write + Send>> {
    let inner: Box<dyn Write + Send> = if path == Path::new("-") {
        Box::new(std::io::stdout())
    } else {
        Box::new(
            fs_err::File::create(path)
                .with_context(|| format!("Failed to create output file {:?}", path))?,
        )
    };

    Ok(match compress {
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
    })
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

/// Creates the multifile output directory. Callers do this on first use rather
/// than up front, so a config the search rejects leaves the filesystem alone.
/// Idempotent, and races benignly between workers.
fn ensure_dir(dir: &Path) -> Result<()> {
    fs_err::create_dir_all(dir)
        .with_context(|| format!("Failed to create output directory {dir:?} for --multifile"))
}

pub(crate) fn build_multifile_paths<'a>(
    names: impl ExactSizeIterator<Item = &'a str>,
    output_dir: &Path,
    ext: &str,
) -> Vec<PathBuf> {
    let mut used_stems: HashSet<String> = HashSet::with_capacity(names.len());
    let mut out: Vec<PathBuf> = Vec::with_capacity(names.len());

    for name in names {
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
