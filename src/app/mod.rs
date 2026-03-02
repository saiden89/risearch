//! Application entry point - handles CLI dispatch and orchestration.

use std::path::Path;

use anyhow::{Context, Result};
use clap::CommandFactory;
use log::{debug, info, trace};

use risearch::{search, QueryRegistry, TargetStore};

use crate::cli::legacy::emit_legacy_warnings;
use crate::cli::{Cli, Commands};

pub(crate) fn run(cli: Cli) -> Result<()> {
    init_logging(cli.verbose);
    init_thread_pool(cli.jobs)?;

    match &cli.command {
        Some(Commands::Index(cmd)) => cmd_index(&cmd.input, &cmd.output),
        Some(Commands::Search(cmd)) => cmd_search(&cmd.query, &cmd.target, &cmd.opts),
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
    cli_opts: &risearch::cli::args::SearchArgs,
) -> Result<()> {
    let output_path = &cli_opts.output.path;
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
        search::run_search_multifile(&queries, &targets, &opts, output_path)?
    } else {
        let mut writer =
            risearch::output::writer::OutputWriter::new(&opts.output, output_path, &queries)?;

        let hits = search::run_search(&queries, &targets, &opts, |chunk| {
            writer.write_chunk(&chunk)
        })?;

        writer.flush_all().context("Failed to flush output")?;
        hits
    };

    info!("Done: {} hits", hits);
    Ok(())
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
