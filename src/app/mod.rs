//! Application entry point - handles CLI dispatch and orchestration.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use clap::CommandFactory;
use log::{debug, info, trace};

use risearch::{output, sa, search, QueryRegistry};

use crate::cli::warnings::emit_legacy_warnings;
use crate::cli::{Cli, Commands};

pub(crate) fn run(cli: Cli) -> Result<()> {
    init_logging(cli.verbose);
    init_thread_pool(cli.threads)?;

    match &cli.command {
        Some(Commands::Index { input, output }) => cmd_index(input, output),
        Some(Commands::Search {
            query,
            target,
            output,
            opts,
        }) => cmd_search(query, target, output, opts),
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
    sa::create_suffix_array(input, output)?;
    info!("Index saved to {:?}", output);
    Ok(())
}

fn cmd_search(
    query_path: &Path,
    index_path: &Path,
    output_path: &Path,
    cli_opts: &risearch::cli_args::SearchArgs,
) -> Result<()> {
    let raw_args: Vec<String> = std::env::args().collect();

    debug!("Loading queries from {:?}", query_path);
    let queries = sa::process_sequences(query_path).context("Failed to load queries")?;

    // Convert CLI args to config (handles deprecated flag translation)
    let output_compress = cli_opts.output_compress;
    let output_level = cli_opts.output_level;
    let mut opts: risearch::config::SearchArgs = cli_opts.clone().into();
    emit_legacy_warnings(&raw_args, &mut opts)?;

    let queries = QueryRegistry::from_indices(queries.into_entries(), &opts.seed);
    info!("Loaded {} queries", queries.len());

    debug!("Loading index from {:?}", index_path);
    let targets = sa::load_index_file(index_path).context("Failed to load index")?;
    trace!("Index loaded: {} targets", targets.len());

    // Open output with compression
    let mut writer =
        output::open_compressed_output(Some(output_path), output_compress, output_level)?;
    opts.output.compress = Some(risearch::config::OutputCompression::None);

    debug!("Starting search...");
    let hits = search::run_search_streaming(&queries, &targets, &opts, &mut writer)?;
    writer.flush().context("Failed to flush output")?;

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
