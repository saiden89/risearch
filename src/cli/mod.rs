pub(crate) mod legacy;

use clap::{Parser, Subcommand};
use std::path::PathBuf;

use risearch::cli::args::SearchArgs;

#[derive(Parser, Debug)]
#[command(name = "risearch")]
#[command(
    author,
    version = "3.alpha.1",
    about = "Energy based RNA-RNA interaction predictions",
    long_about = None,
    subcommand_required = true,
    arg_required_else_help = true,
    after_help = "Subcommand help:\n  risearch index --help\n  risearch search --help"
)]
pub(crate) struct Cli {
    /// Increase verbosity (-v = info, -vv = debug, -vvv = trace)
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, global = true)]
    pub(crate) verbose: u8,

    /// Set number of parallel jobs (global)
    #[arg(
        short = 'j',
        long = "jobs",
        alias = "threads",
        value_name = "N",
        default_value_t = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1),
        global = true
    )]
    pub(crate) jobs: usize,

    #[command(subcommand)]
    pub(crate) command: Commands,
}

#[derive(clap::Args, Debug)]
pub(crate) struct IndexCommand {
    /// Input file in FASTA format.
    #[arg(value_name = "INPUT")]
    pub(crate) input: PathBuf,

    /// Save index to given index file path
    #[arg(value_name = "OUTPUT")]
    pub(crate) output: PathBuf,
}

#[derive(Subcommand, Debug)]
#[allow(clippy::large_enum_variant)] // CLI parsing - allocation overhead is negligible
pub(crate) enum Commands {
    /// Create index for target sequence(s)
    Index(IndexCommand),

    /// Search for interactions in the given sequence(s)
    Search(SearchArgs),
}
