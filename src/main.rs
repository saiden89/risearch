#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

mod app;
mod cli;

use anyhow::Result;
use clap::Parser;
use crate::cli::Cli;

fn main() -> Result<()> {
    let cli = Cli::parse();
    app::run(cli)
}
