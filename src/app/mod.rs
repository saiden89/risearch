//! Application entry point - handles CLI dispatch and orchestration.

use std::path::Path;

use anyhow::{Context, Result};
use log::{debug, info, trace};

use risearch::{search, QueryRegistry, SearchConfig, TargetRegistry};

use crate::cli::legacy::emit_legacy_warnings;
use crate::cli::{Cli, Commands, SearchArgs};

pub(crate) fn run(cli: Cli) -> Result<()> {
    match &cli.command {
        Commands::Index(cmd) => {
            init_runtime(cli.verbose, cli.jobs)?;
            cmd_index(&cmd.input, &cmd.output, cli.jobs)
        }
        Commands::Search(cmd) => {
            init_runtime(cli.verbose, cli.jobs)?;
            cmd_search(cmd)
        }
    }
}

// =============================================================================
// COMMAND HANDLERS
// =============================================================================

fn cmd_index(input: &Path, output: &Path, threads: Option<usize>) -> Result<()> {
    info!("Creating index: {:?} -> {:?}", input, output);
    TargetRegistry::build(input, output, threads).context("Failed to write index file")?;
    info!("Index saved to {:?}", output);
    Ok(())
}

fn cmd_search(cmd: &SearchArgs) -> Result<()> {
    let query_path = &cmd.input.query;
    let target_path = cmd.input.target_path();
    let output_path = &cmd.output.path;

    // Convert CLI args to config (handles deprecated flag translation)
    let opts: SearchConfig = cmd.clone().try_into()?;
    emit_legacy_warnings(cmd, cmd.input.uses_legacy_target());

    debug!("Loading queries from {:?}", query_path);
    let queries =
        QueryRegistry::from_fasta(query_path, &opts.seed).context("Failed to load queries")?;
    info!("Loaded {} queries", queries.len());

    debug!("Loading target index from {:?}", target_path);
    let targets = TargetRegistry::open(target_path).context("Failed to load index")?;
    trace!("Index loaded: {} targets", targets.len());

    debug!("Starting search...");

    search::run_search(&queries, &targets, &opts, output_path)?;
    info!("Done");
    Ok(())
}

// =============================================================================
// INITIALIZATION
// =============================================================================

fn init_runtime(verbosity: u8, threads: Option<usize>) -> Result<()> {
    init_logging(verbosity);
    init_thread_pool(threads)
}

fn init_thread_pool(threads: Option<usize>) -> Result<()> {
    let mut builder = rayon::ThreadPoolBuilder::new();
    // Only pin the count when the user passed one; otherwise let rayon read
    // RAYON_NUM_THREADS and fall back to all cores on its own.
    if let Some(n) = threads {
        builder = builder.num_threads(n);
    }
    builder
        .build_global()
        .context("Failed to initialize thread pool")?;
    debug!("Thread pool: {} threads", rayon::current_num_threads());
    Ok(())
}

pub(crate) fn init_logging(verbosity: u8) {
    let level = match verbosity {
        0 => log::LevelFilter::Warn,
        1 => log::LevelFilter::Info,
        2 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };

    env_logger::Builder::new()
        .filter_level(level)
        .format_timestamp(None)
        .format_target(false)
        .format_module_path(false)
        .init();
}
