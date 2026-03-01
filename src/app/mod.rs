//! Application entry point - handles CLI dispatch and orchestration.

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::CommandFactory;
use log::{debug, info, trace};

use risearch::{output, search, QueryRegistry, TargetStore};

use crate::cli::legacy::emit_legacy_warnings;
use crate::cli::{Cli, Commands};

pub(crate) fn run(cli: Cli) -> Result<()> {
    init_logging(cli.verbose);
    init_thread_pool(cli.jobs)?;

    match &cli.command {
        Some(Commands::Index(cmd)) => cmd_index(&cmd.input, &cmd.output),
        Some(Commands::Search(cmd)) => {
            cmd_search(&cmd.query, &cmd.target, &cmd.output.path, &cmd.opts)
        }
        None => {
            Cli::command().print_help()?;
            println!();
            Ok(())
        }
    }
}

// =============================================================================
// COMMAND HANDLERS
// =============================================================================

fn cmd_index(input: &Path, output: &Path) -> Result<()> {
    info!("Creating index: {:?} -> {:?}", input, output);
    TargetStore::build_from_fasta(input, output).context("Failed to write index file")?;
    info!("Index saved to {:?}", output);
    Ok(())
}

fn cmd_search(
    query_path: &Path,
    target_path: &Path,
    output_path: &Path,
    cli_opts: &risearch::cli::args::SearchArgs,
) -> Result<()> {
    let raw_args: Vec<String> = std::env::args().collect();

    // Convert CLI args to config (handles deprecated flag translation)
    let mut opts: risearch::config::SearchArgs = cli_opts.clone().into();
    emit_legacy_warnings(&raw_args, &mut opts)?;

    debug!("Loading queries from {:?}", query_path);
    let queries =
        QueryRegistry::from_fasta(query_path, &opts.seed).context("Failed to load queries")?;
    info!("Loaded {} queries", queries.len());

    debug!("Loading target index from {:?}", target_path);
    let targets = TargetStore::open(target_path).context("Failed to load index")?;
    trace!("Index loaded: {} targets", targets.len());

    debug!("Starting search...");

    let hits = if opts.output.multifile {
        // -- Multi-file mode: one output file per query ----------------------
        std::fs::create_dir_all(output_path).with_context(|| {
            format!(
                "Failed to create output directory {:?} for --output-multifile",
                output_path
            )
        })?;

        let ext = output::output_extension(&opts.output);
        let query_output_paths = build_multifile_paths(&queries, output_path, ext);
        let mut writers: HashMap<u32, Box<dyn Write>> = HashMap::new();

        let hits = search::run_search(&queries, &targets, &opts, |chunk| -> Result<()> {
            let qi = chunk
                .query_idx
                .expect("BUG: multifile mode active but OutputChunk has no query_idx");
            let writer = match writers.get_mut(&qi) {
                Some(w) => w,
                None => {
                    let file_path = &query_output_paths[qi as usize];
                    let w = output::open_output(Some(&file_path), &opts.output)?;
                    writers.entry(qi).or_insert(w)
                }
            };
            writer
                .write_all(&chunk.data)
                .context("Failed to write output chunk")?;
            Ok(())
        })?;

        for (_, mut w) in writers.drain() {
            w.flush().context("Failed to flush per-query output")?;
        }

        hits
    } else {
        // -- Single-file mode (default) --------------------------------------
        let mut writer = output::open_output(Some(output_path), &opts.output)?;

        let hits = search::run_search(&queries, &targets, &opts, |chunk| -> Result<()> {
            writer
                .write_all(&chunk.data)
                .context("Failed to write output chunk")?;
            Ok(())
        })?;

        writer.flush().context("Failed to flush output")?;
        hits
    };

    info!("Done: {} hits", hits);
    Ok(())
}

/// Replace characters that are unsafe in filenames.
fn sanitize_filename(name: &str) -> String {
    name.chars()
        .map(|c| match c {
            '/' | '\\' | '\0' | ':' | '*' | '?' | '"' | '<' | '>' | '|' => '_',
            _ => c,
        })
        .collect()
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
        let stem_base = if stem_raw.is_empty() {
            "query"
        } else {
            stem_raw.as_str()
        };
        let stem = unique_filename_stem(stem_base, &mut used_stems);
        out.push(output_dir.join(format!("{stem}{ext}")));
    }

    out
}

// =============================================================================
// INITIALIZATION
// =============================================================================

fn init_thread_pool(threads: usize) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .context("Failed to initialize thread pool")?;
    debug!("Thread pool: {} threads", threads);
    Ok(())
}

pub(crate) fn init_logging(verbosity: u8) {
    use std::io::Write;

    let level = match verbosity {
        0 => log::LevelFilter::Warn,
        1 => log::LevelFilter::Info,
        2 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };

    env_logger::Builder::new()
        .filter_level(level)
        .format_timestamp(None)
        .format(|buf, record| {
            let color = match record.level() {
                log::Level::Error => "\x1b[31m", // red
                log::Level::Warn => "\x1b[33m",  // yellow
                log::Level::Info => "\x1b[32m",  // green
                log::Level::Debug => "\x1b[34m", // blue
                log::Level::Trace => "\x1b[35m", // magenta
            };
            writeln!(
                buf,
                "[\x1b[1m{}{:<5}\x1b[0m] {}",
                color,
                record.level(),
                record.args()
            )
        })
        .init();
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
