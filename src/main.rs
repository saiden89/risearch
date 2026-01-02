mod dsm;
#[cfg(feature = "fm-index")]
mod fm;
mod io;
mod sa;
mod search;
mod seed;

use anyhow::{Context, Result, bail};
use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use log::{debug, info, trace};
use std::path::PathBuf;

/// Initialize logging based on verbosity level
fn init_logging(verbosity: u8) {
    let level = match verbosity {
        0 => log::LevelFilter::Warn,
        1 => log::LevelFilter::Info,
        2 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };
    env_logger::Builder::new()
        .filter_level(level)
        .format_timestamp(None)
        .init();
}

#[derive(ValueEnum, Clone, Debug)]
#[clap(rename_all = "snake_case")]
enum Matrix {
    T99,
    T04,
}

#[derive(clap::ValueEnum, Clone, Debug)]
#[value(rename_all = "lowercase")]
enum OutputFormat {
    /// Report predictions in detailed format
    Detailed,
    /// Report predictions in a simple format together with CIGAR-like string for interaction structure
    Cigar,
    /// Report predictions in a simple format together with binding site (3'->5'), flanking 5'end (3'->5') and flanking 3'end (5'->3') sequences of the target x(required for post-processing of CRISPR off-target predictions)
    BindingSite,
}

#[derive(ValueEnum, Clone, Debug)]
#[clap(rename_all = "lowercase")]
enum Backend {
    /// FM-Index backend (requires 'fm-index' feature)
    Fm,
    /// Suffix Array backend
    Sa,
}

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

/// Arguments for seed generation
#[derive(clap::Args, Debug, Clone)]
pub struct SeedArgs {
    /// Set seed length (-s l = length only; -s n:m = full interval)
    #[arg(
        short = 's',
        long = "seed",
        value_name = "start:end/length",
        default_value = "6"
    )]
    seed: Option<String>,

    /// Consider G-U wobble pairs as mismatch within the seed (only for locating seeds, energy model is not affected)
    #[arg(short = 'w', long = "wobble", action = clap::ArgAction::SetTrue)]
    wobble: bool,

    /// Introduce mismatched seeds
    /// Set the max num of mismatches (c) allowed in the seed and min num of consecutive matches required at seed start/end (p)
    ///These seeds will not overlap with perfect complementary seeds.
    #[arg(
        short = 'm',
        long = "mismatch",
        value_name = "c:p",
        default_value = "0:0"
    )]
    mismatch_seed: Option<String>,
}

/// Arguments for seed extension and scoring
#[derive(clap::Args, Debug, Clone)]
pub struct ExtendArgs {
    /// Max extension length on the seed (do DP for max this length up- and downstream of seed)
    #[arg(
        short = 'l',
        long = "extension",
        value_name = "LENGTH",
        default_value_t = 20
    )]
    max_extension: u8,

    /// Set deltaG energy threshold (in kcal/mol) to filter predictions
    #[arg(
            short = 'e',
            long = "energy",
            value_name = "dG",
            default_value_t = -20.0,
            allow_hyphen_values = true
        )]
    delta_g: f64,

    /// Energy matrix for RNA-RNA duplexes
    #[arg(short = 'z', long = "matrix", value_name = "MATRIX", default_value_t = Matrix::T04, value_enum)]
    matrix: Matrix,

    /// Per-nucleotide extension penalty (in kcal/mol)
    #[arg(
        short = 'p',
        long = "penalty",
        value_name = "PENALTY",
        default_value_t = 0.0
    )]
    penalty: f64,

    /// Energy per length threshold that filters seeds
    #[arg(
        short = 'x',
        long = "seed-energy",
        value_name = "THRESHOLD",
        default_value_t = 0.0
    )]
    seed_energy: f64,
}

/// Options that apply to the `search` subcommand
#[derive(clap::Args, Debug, Clone)]
pub struct SearchArgs {
    #[command(flatten)]
    seed: SeedArgs,

    #[command(flatten)]
    extend: ExtendArgs,

    /// Output format
    #[arg(
        short = 'f',
        long = "format",
        value_name = "LEVEL",
        default_missing_value = "detailed",
        default_value = "detailed",
        require_equals = true,
        value_enum
    )]
    report_format: Option<OutputFormat>,
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
                    search::run_search(&queries, &wrapper, output, opts)?;
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
