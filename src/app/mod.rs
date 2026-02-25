//! Application entry point - handles CLI dispatch and orchestration.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use clap::CommandFactory;
use log::{debug, info, trace};

use risearch::{output, search, QueryRegistry, TargetStore};

use crate::cli::legacy::emit_legacy_warnings;
use crate::cli::{Cli, Commands};

pub(crate) fn run(cli: Cli) -> Result<()> {
    init_logging(cli.verbose);
    init_thread_pool(cli.jobs)?;

    match &cli.command {
        Some(Commands::Index(cmd)) => cmd_index(&cmd.input, &cmd.output),
        Some(Commands::Search(cmd)) => {
            cmd_search(&cmd.query, &cmd.target, &cmd.output.path, &cmd.opts)
        }
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
    TargetStore::build_from_fasta(input, output).context("Failed to write index file")?;
    info!("Index saved to {:?}", output);
    Ok(())
}

fn cmd_search(
    query_path: &Path,
    target_path: &Path,
    output_path: &Path,
    cli_opts: &risearch::cli::args::SearchArgs,
) -> Result<()> {
    let raw_args: Vec<String> = std::env::args().collect();

    // Convert CLI args to config (handles deprecated flag translation)
    let mut opts: risearch::config::SearchArgs = cli_opts.clone().into();
    emit_legacy_warnings(&raw_args, &mut opts)?;

    debug!("Loading queries from {:?}", query_path);
    let queries =
        QueryRegistry::from_fasta(query_path, &opts.seed).context("Failed to load queries")?;
    info!("Loaded {} queries", queries.len());

    debug!("Loading target index from {:?}", target_path);
    let targets = TargetStore::open(target_path).context("Failed to load index")?;
    trace!("Index loaded: {} targets", targets.len());

    // Open output with compression
    let mut writer = output::open_output(Some(output_path), &opts.output)?;
    let mut out_buf = Vec::with_capacity(64 * 1024);
    let mut fmt_bufs = output::OutputBuffers::new();
    let format = opts.output.format;
    let mut target_cache: Option<(u32, risearch::index::store::TargetView<'_>)> = None;

    debug!("Starting search...");
    let hits = search::run_search(&queries, &targets, &opts, |hit| -> Result<()> {
        if target_cache.as_ref().map(|(idx, _)| *idx) != Some(hit.target_idx) {
            let view = targets
                .target_view(hit.target_idx as usize)
                .with_context(|| format!("Failed to load target #{}", hit.target_idx))?;
            target_cache = Some((hit.target_idx, view));
        }

        let target = &target_cache
            .as_ref()
            .expect("target cache must be populated")
            .1;
        let q_name = queries.get_name(hit.query_idx);
        let q_seq = queries.get(hit.query_idx).sequence();
        let t_fwd = &target.combined_seq[..target.seq_len];
        let t_rc = &target.combined_seq[target.seq_len + 1..2 * target.seq_len + 1];

        if format == risearch::config::OutputFormat::Minimal {
            output::append_hit_minimal_names_vec(
                &mut fmt_bufs,
                &hit,
                &mut out_buf,
                q_name,
                target.name,
            );
        } else {
            output::write_hit_names(
                &mut fmt_bufs,
                &hit,
                format,
                &mut out_buf,
                q_name,
                target.name,
                q_seq,
                t_fwd,
                t_rc,
            )
            .context("Failed to format output hit")?;
        }

        if out_buf.len() >= 64 * 1024 {
            writer
                .write_all(&out_buf)
                .context("Failed to write output chunk")?;
            out_buf.clear();
        }

        Ok(())
    })?;

    if !out_buf.is_empty() {
        writer
            .write_all(&out_buf)
            .context("Failed to write output chunk")?;
        out_buf.clear();
    }
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
