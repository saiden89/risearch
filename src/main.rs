#[cfg(feature = "fm-index")]
use risearch::fm;
use risearch::{sa, search};

use anyhow::{Context, Result, bail};
use clap::{CommandFactory, Parser, Subcommand};
use log::{debug, info, trace};
use std::path::PathBuf;

/// Initialize logging based on verbosity level with colored output
fn init_logging(verbosity: u8) {
    use std::io::Write;

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

use risearch::args::*;

#[derive(Parser, Debug)]
#[command(name = "risearch")]
#[command(author, version = "3.alpha.1", about = "Energy based RNA-RNA interaction predictions", long_about = None)]
struct Cli {
    /// Increase verbosity (-v = info, -vv = debug, -vvv = trace)
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, global = true)]
    verbose: u8,

    /// Set threads for parallel processing (global)
    #[arg(short = 't', long = "threads", value_name = "N", default_value_t = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1), global = true)]
    threads: usize,

    /// Choose backend (fm or sa)
    #[arg(short = 'b', long = "backend", value_enum, default_value_t = Backend::Sa, global = true)]
    backend: Backend,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
#[allow(clippy::large_enum_variant)] // CLI parsing - allocation overhead is negligible
enum Commands {
    /// Create index for target sequence(s)
    Index {
        /// Input file in FASTA format.
        #[arg(value_name = "INPUT")]
        input: PathBuf,

        /// Save index to given index file path
        #[arg(value_name = "OUTPUT")]
        output: PathBuf,
    },

    /// Search for interactions in the given sequence(s)
    Search {
        /// FASTA file for query sequence(s) (.fa or .fa.gz) -- use '-' for stdin
        #[arg(short = 'q', long = "query", value_name = "FILE")]
        query: PathBuf,

        /// Index file (created by `index` command)
        #[arg(short = 'i', long = "index", value_name = "INDEX")]
        index: PathBuf,

        /// Output file for search results (use '-' for stdout)
        #[arg(short = 'o', long = "output", value_name = "FILE")]
        output: PathBuf,

        /// Search-related options (seed, extension, energy, matrix, penalty, threads, format)
        #[command(flatten)]
        opts: SearchArgs,
    },
}

fn main() -> Result<()> {
    let cli: Cli = Cli::parse();

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
            info!(
                "Creating index ({:?}) from {:?} -> {:?}",
                cli.backend, input, output
            );
            match cli.backend {
                Backend::Fm => {
                    #[cfg(feature = "fm-index")]
                    {
                        fm::create_fm_index(input, output)?;
                    }
                    #[cfg(not(feature = "fm-index"))]
                    {
                        bail!(
                            "FM-Index backend is not available. \
                            Recompile with `--features fm-index` to enable it."
                        );
                    }
                }
                Backend::Sa => sa::create_suffix_array(input, output)?,
            }
            info!("Saved index to {:?}", output);
        }
        Some(Commands::Search {
            query,
            index,
            output,
            opts,
        }) => {
            debug!(
                "Search command: backend={:?} query={:?} index={:?} output={:?}",
                cli.backend, query, index, output
            );
            trace!("Search options: {:?}", opts);

            let queries = sa::process_sequences(query)
                .context("Failed to process query sequences")?
                .sequences
                .into_iter()
                .map(|s| (s.name, s.sequence))
                .collect::<Vec<_>>();

            info!("Loaded {} query sequences", queries.len());

            match cli.backend {
                Backend::Fm => {
                    #[cfg(feature = "fm-index")]
                    {
                        warn!("FM-Index search is experimental");
                        // TODO: Implement FM-Index search
                        bail!("FM-Index search not yet implemented");
                    }
                    #[cfg(not(feature = "fm-index"))]
                    {
                        bail!(
                            "FM-Index backend is not available. \
                            Recompile with `--features fm-index` to enable it."
                        );
                    }
                }
                Backend::Sa => {
                    debug!("Loading suffix array index...");
                    let idx = sa::load_index_file(index).context("Failed to load index file")?;
                    trace!("Index loaded successfully");

                    let wrapper = search::SaIndex { index: &idx };
                    debug!("Starting search...");
                    let hits = search::run_search(&queries, &wrapper, opts)?;
                    search::write_results(&hits, output)?;
                    info!("Search completed successfully");
                }
            }
        }
        None => {
            // No subcommand: show help
            Cli::command().print_help()?;
            println!();
        }
    }

    Ok(())
}
