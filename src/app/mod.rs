//! Application entry point - handles CLI dispatch and orchestration.

use std::collections::HashMap;
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
        let mut writers: HashMap<u32, Box<dyn Write>> = HashMap::new();

        let hits = search::run_search(&queries, &targets, &opts, |chunk| -> Result<()> {
            let qi = chunk
                .query_idx
                .expect("BUG: multifile mode active but OutputChunk has no query_idx");
            let writer = match writers.get_mut(&qi) {
                Some(w) => w,
                None => {
                    let qname = sanitize_filename(queries.get_name(qi));
                    let file_path: PathBuf = [output_path, Path::new(&format!("{qname}{ext}"))]
                        .iter()
                        .collect();
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
