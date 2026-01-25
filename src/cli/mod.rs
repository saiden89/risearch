pub mod warnings;

use clap::{Parser, Subcommand};
use std::path::PathBuf;

use risearch::args::SearchArgs;

#[derive(Parser, Debug)]
#[command(name = "risearch")]
#[command(
    author,
    version = "3.alpha.1",
    about = "Energy based RNA-RNA interaction predictions",
    long_about = None,
    after_help = "Subcommand help:\n  risearch index --help\n  risearch search --help"
)]
pub struct Cli {
    /// Increase verbosity (-v = info, -vv = debug, -vvv = trace)
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, global = true)]
    pub verbose: u8,

    /// Set threads for parallel processing (global)
    #[arg(
        short = 't',
        long = "threads",
        value_name = "N",
        default_value_t = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1),
        global = true
    )]
    pub threads: usize,

    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
#[allow(clippy::large_enum_variant)] // CLI parsing - allocation overhead is negligible
pub enum Commands {
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
