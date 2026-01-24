#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[cfg(feature = "fm-index")]
use risearch::fm;
use risearch::{sa, search};
use risearch::seed::MismatchSpec;

use anyhow::{Context, Result, bail};
use clap::{CommandFactory, Parser, Subcommand};
use log::{debug, info, trace, warn};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

/// Initialize logging based on verbosity level with colored output
fn init_logging(verbosity: u8) {
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

fn infer_compression(path: &std::path::Path) -> OutputCompression {
    let Some(ext) = path.extension().and_then(|e| e.to_str()) else {
        return OutputCompression::None;
    };
    match ext.to_ascii_lowercase().as_str() {
        "gz" | "gzip" => OutputCompression::Gzip,
        "zst" | "zstd" => OutputCompression::Zstd,
        _ => OutputCompression::None,
    }
}

fn extract_legacy_mismatch_arg(args: &[String]) -> Option<String> {
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        if arg == "-m" || arg == "--mismatch" {
            if let Some(val) = iter.next() {
                return Some(val.clone());
            }
        } else if let Some(val) = arg.strip_prefix("--mismatch=") {
            return Some(val.to_string());
        } else if arg.starts_with("-m") && arg.len() > 2 {
            return Some(arg[2..].to_string());
        }
    }
    None
}

fn has_seed_pairing_arg(args: &[String]) -> bool {
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        if arg == "--seed-pairing" {
            if iter.peek().is_some() {
                return true;
            }
        } else if arg.starts_with("--seed-pairing=") {
            return true;
        }
    }
    false
}

fn has_wobble_legacy_arg(args: &[String]) -> bool {
    args.iter().skip(1).any(|arg| arg == "-w" || arg == "--wobble")
}

#[derive(Parser, Debug)]
#[command(name = "risearch")]
#[command(
    author,
    version = "3.alpha.1",
    about = "Energy based RNA-RNA interaction predictions",
    long_about = None,
    after_help = "Subcommand help:\n  risearch index --help\n  risearch search --help"
)]
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
                    debug!("Starting search with streaming output...");

                    let output_path: &std::path::Path = output.as_ref();
                    let inner: Box<dyn std::io::Write> = if output_path == std::path::Path::new("-")
                    {
                        Box::new(std::io::stdout())
                    } else {
                        Box::new(
                            std::fs::File::create(output_path)
                                .context("Failed to create output file")?,
                        )
                    };
                    let compression = opts.output_compress.unwrap_or_else(|| {
                        if output_path == std::path::Path::new("-") {
                            OutputCompression::None
                        } else {
                            infer_compression(output_path)
                        }
                    });
                    let mut writer = BufWriter::with_capacity(256 * 1024, inner);
                    let mut opts = opts.clone();
                    let legacy_mismatch = extract_legacy_mismatch_arg(&raw_args);
                    let legacy_wobble = has_wobble_legacy_arg(&raw_args);
                    let explicit_pairing = has_seed_pairing_arg(&raw_args);
                    if let Some(raw) = legacy_mismatch.as_deref() {
                        if opts.seed.has_named_mismatch() {
                            warn!(
                                "Both legacy -m/--mismatch ({}) and named --mismatch-* flags were provided; named flags take precedence.",
                                raw
                            );
                        }
                        let suggestion = raw
                            .parse::<MismatchSpec>()
                            .ok()
                            .map(|spec| {
                                format!(
                                    "--mismatch-max {} --mismatch-prefix {} --mismatch-suffix {}",
                                    spec.max_mismatches, spec.min_position, spec.min_matches_after
                                )
                            });
                        if let Some(s) = suggestion {
                            warn!(
                                "Legacy mismatch syntax '-m c[:ps[:pe]]' is deprecated; use {}.",
                                s
                            );
                        } else {
                            warn!(
                                "Legacy mismatch syntax '-m c[:ps[:pe]]' is deprecated; use --mismatch-max/--mismatch-prefix/--mismatch-suffix."
                            );
                        }
                    }
                    if legacy_wobble {
                        warn!(
                            "Legacy -w/--wobble is deprecated and redundant (default is allow_wobble). Prefer --seed-pairing allow_wobble."
                        );
                        if opts.seed.no_guseed {
                            warn!(
                                "Both -w/--wobble and --no-guseed were provided; --no-guseed takes precedence."
                            );
                        } else if explicit_pairing {
                            warn!(
                                "Both -w/--wobble and --seed-pairing were provided; --seed-pairing takes precedence."
                            );
                        }
                    }
                    opts.seed.apply_mismatch_overrides();
                    opts.seed.apply_pairing_overrides(explicit_pairing);
                    opts.output_compress = Some(compression);
                    let hit_count =
                        search::run_search_streaming(&queries, &wrapper, &opts, &mut writer)?;
                    writer.flush().context("Failed to flush output")?;
                    info!("Search completed: {} hits written", hit_count);
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
