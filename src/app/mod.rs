use std::io::Write;

use anyhow::{Context, Result};
use clap::CommandFactory;
use log::{debug, info, trace};

use risearch::{io, sa, search};

use crate::cli::{Cli, Commands};
use crate::cli::warnings::emit_legacy_warnings;

/// Initialize logging based on verbosity level with colored output.
pub fn init_logging(verbosity: u8) {
    let level = match verbosity {
        0 => log::LevelFilter::Warn,
        1 => log::LevelFilter::Info,
        2 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };

    // ANSI color codes (matching tracing palette)
    const RESET: &str = "\x1b[0m";
    const BOLD: &str = "\x1b[1m";
    const RED: &str = "\x1b[31m";
    const GREEN: &str = "\x1b[32m";
    const YELLOW: &str = "\x1b[33m";
    const BLUE: &str = "\x1b[34m";
    const MAGENTA: &str = "\x1b[35m";
    const CYAN: &str = "\x1b[36m";

    env_logger::Builder::new()
        .filter_level(level)
        .format_timestamp(None)
        .format(move |buf, record| {
            let msg = record.args().to_string();

            // Color code for log level (matching tracing: trace=purple)
            let level_color = match record.level() {
                log::Level::Error => RED,
                log::Level::Warn => YELLOW,
                log::Level::Info => GREEN,
                log::Level::Debug => BLUE,
                log::Level::Trace => MAGENTA, // Purple like tracing
            };

            // Color code for message based on component prefix
            // Messages without recognized prefix use level color for consistency
            let msg_color = if msg.starts_with("DP_LEFT") || msg.starts_with("DP_RIGHT") {
                CYAN
            } else if msg.starts_with("SA_SEARCH") || msg.starts_with("[FIND_CAND]") {
                YELLOW
            } else if msg.starts_with("SEED") {
                GREEN
            } else if msg.starts_with("EXTEND") || msg.starts_with("[MAXIMALITY]") {
                MAGENTA
            } else if msg.starts_with("DEDUP")
                || msg.starts_with("[HIT")
                || msg.starts_with("[PROC_CAND]")
            {
                BLUE
            } else if msg.starts_with("[QUERY]")
                || msg.starts_with("Starting")
                || msg.starts_with("Search")
            {
                GREEN
            } else {
                // Fallback: use level color for consistent appearance
                level_color
            };

            writeln!(
                buf,
                "[{}{}{:<5}{}] {}{}{}",
                BOLD,
                level_color,
                record.level(),
                RESET,
                msg_color,
                msg,
                RESET
            )
        })
        .init();
}

pub fn run(cli: Cli) -> Result<()> {
    let raw_args: Vec<String> = std::env::args().collect();

    // Initialize logging based on verbosity
    init_logging(cli.verbose);

    // Configure rayon thread pool globally
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;

    debug!("Thread pool initialized with {} threads", cli.threads);
    trace!("CLI arguments: {:?}", cli);

    match &cli.command {
        Some(Commands::Index { input, output }) => {
            info!("Creating index from {:?} -> {:?}", input, output);
            sa::create_suffix_array(input, output)?;
            info!("Saved index to {:?}", output);
        }
        Some(Commands::Search {
            query,
            index,
            output,
            opts,
        }) => {
            debug!(
                "Search command: query={:?} index={:?} output={:?}",
                query, index, output
            );
            trace!("Search options: {:?}", opts);

            let queries = sa::process_sequences(query)
                .context("Failed to process query sequences")?
                .sequences
                .into_iter()
                .map(|s| (s.name, s.sequence))
                .collect::<Vec<_>>();

            info!("Loaded {} query sequences", queries.len());

            debug!("Loading suffix array index...");
            let idx = sa::load_index_file(index).context("Failed to load index file")?;
            trace!("Index loaded successfully");

            let wrapper = search::SaIndex { index: &idx };
            debug!("Starting search with streaming output...");

            let (mut writer, compression) =
                io::output::open_output(output.as_ref(), opts.output_compress.map(Into::into))?;
            let mut opts: risearch::config::SearchArgs = opts.clone().into();
            emit_legacy_warnings(&raw_args, &mut opts);
            opts.output.compress = Some(compression);

            let hit_count =
                search::run_search_streaming(&queries, &wrapper, &opts, &mut writer)?;
            writer.flush().context("Failed to flush output")?;
            info!("Search completed: {} hits written", hit_count);
        }
        None => {
            // No subcommand: show help
            Cli::command().print_help()?;
            println!();
        }
    }

    Ok(())
}
