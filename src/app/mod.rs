//! Application entry point - handles CLI dispatch and orchestration.

use std::path::Path;

use anyhow::{Context, Result};
use log::{debug, info, trace};

use risearch::{search, QueryRegistry, SearchConfig, TargetStore};

use crate::cli::legacy::emit_legacy_warnings;
use crate::cli::{Cli, Commands, SearchArgs};

pub(crate) fn run(cli: Cli) -> Result<()> {
    match &cli.command {
        Commands::Index(cmd) => {
            init_runtime(cli.verbose, cli.jobs)?;
            cmd_index(&cmd.input, &cmd.output)
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

fn cmd_index(input: &Path, output: &Path) -> Result<()> {
    info!("Creating index: {:?} -> {:?}", input, output);
    TargetStore::build_from_fasta(input, output).context("Failed to write index file")?;
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
    let targets = TargetStore::open(target_path).context("Failed to load index")?;
    trace!("Index loaded: {} targets", targets.len());

    debug!("Starting search...");

    search::run_search(&queries, &targets, &opts, output_path)?;
    info!("Done");
    Ok(())
}

// =============================================================================
// INITIALIZATION
// =============================================================================

fn init_runtime(verbosity: u8, threads: usize) -> Result<()> {
    init_logging(verbosity);
    init_thread_pool(threads)
}

fn init_thread_pool(threads: usize) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .context("Failed to initialize thread pool")?;
    debug!("Thread pool: {} threads", threads);
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
