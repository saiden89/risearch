use anyhow::{Result, bail};
use log::info;
use rayon::prelude::*;
use std::io::Write;

use crate::config::{OutputFormat, SearchArgs};
use crate::dsm::DsmModel;
use crate::registry::QueryRegistry;
use crate::seed::SeedHit;
use crate::seed::search::find_seeds;
use crate::types::Strand;

use super::core::extend_seed;
use super::{FilterReason, SearchContext, THREAD_MATCHES, THREAD_SEEDS};
use crate::output::format::fill_line_buf;
use crate::sa::TargetRegistry;

pub(super) fn run_search_streaming_impl<M: DsmModel, W: std::io::Write>(
    queries: &QueryRegistry,
    index: &TargetRegistry,
    opts: &SearchArgs,
    writer: &mut W,
    extender_tls: &'static std::thread::LocalKey<std::cell::RefCell<crate::dp::DpExtender<M>>>,
) -> Result<usize> {

    for (q_idx, q) in queries.iter() {
        if let Err(err) = opts.seed.seed.normalize(q.sequence().len()) {
            bail!(
                "Invalid seed spec for query '{}': {}",
                queries.get_name(q_idx),
                err
            );
        }
    }
    info!(
        "Starting streaming search: {} queries, seed={:?}, max_ext={}, delta_g={}",
        queries.len(),
        opts.seed.seed,
        opts.extend.max_extension,
        opts.extend.delta_g
    );

    let mut hit_count = 0;
    let format = opts.output.format.unwrap_or(OutputFormat::Detailed);
    let chunk_size = std::cmp::max(1, rayon::current_num_threads() * 4);

    let process_query = |q_idx: u32,
                         q: &crate::registry::QueryEntry,
                         writer: &mut dyn Write,
                         line_buf: &mut Vec<u8>|
     -> usize {
        let mut local_hits = 0usize;
        extender_tls.with(|ext| {
            THREAD_SEEDS.with(|seeds_cell| {
                THREAD_MATCHES.with(|matches_cell| {
                    let mut extender = ext.borrow_mut();
                    let mut seeds = seeds_cell.borrow_mut();
                    let mut matches = matches_cell.borrow_mut();
                    let mut ctx = SearchContext::<M>::with_extender(index, opts, &mut extender);

                    // Reuse thread-local Vecs instead of allocating new ones
                    find_seeds(q, ctx.index, &ctx.args.seed, &mut seeds, &mut matches);

                    for candidate in seeds.iter() {
                        if process_candidate_streaming::<M, _>(
                            q_idx, queries, candidate, &mut ctx, writer, format, line_buf,
                        ) {
                            local_hits += 1;
                        }
                    }
                });
            });
        });
        local_hits
    };

    let mut chunk_start = 0usize;
    for chunk in queries.entries().chunks(chunk_size) {
        // Parallel processing: each thread produces raw bytes (no compression here)
        // Compression is handled by the streaming writer wrapper
        let outputs: Vec<(Vec<u8>, usize)> = chunk
            .par_iter()
            .enumerate()
            .map(|(i, q)| {
                let mut line_buf: Vec<u8> = Vec::new();
                let mut out: Vec<u8> = Vec::new();
                let q_idx = (chunk_start + i) as u32;
                let hits = process_query(q_idx, q, &mut out, &mut line_buf);
                (out, hits)
            })
            .collect();

        for (buf, count) in outputs {
            if !buf.is_empty() {
                writer.write_all(&buf)?;
            }
            hit_count += count;
        }
        chunk_start += chunk.len();
    }

    info!("Streaming search complete: {} hits", hit_count);
    Ok(hit_count)
}

/// Process candidate and write directly to output - zero String allocation variant.
/// Returns true if a hit was written.
fn process_candidate_streaming<M: DsmModel, W: Write + ?Sized>(
    query_idx: u32,
    queries: &QueryRegistry,
    candidate: &SeedHit,
    ctx: &mut SearchContext<'_, '_, M>,
    writer: &mut W,
    format: OutputFormat,
    line_buf: &mut Vec<u8>,
) -> bool {
    let query_seq = queries.get(query_idx).sequence();
    let t_idx = candidate.target_idx;
    let t_seq = match candidate.strand {
        Strand::Reverse => ctx.index.get_sequence_rc(t_idx),
        Strand::Forward => ctx.index.get_sequence(t_idx),
    };

    let t_start_idx = candidate.target_start;
    let seed_len = candidate.len;
    let q_pos = candidate.query_pos;

    if t_start_idx + seed_len > t_seq.len() {
        ctx.stats.record_filter(FilterReason::SeedOutOfBounds);
        return false;
    }

    // Call extend_seed
    let Some(ext) = extend_seed(ctx, query_seq, t_seq, candidate) else {
        return false;
    };

    let score = ext.score;
    if score > ctx.args.extend.delta_g {
        ctx.stats.record_filter(FilterReason::EnergyAboveThreshold);
        return false;
    }

    // Final coordinates (0-based)
    let final_q_start = q_pos - ext.l_q;
    let final_q_end = (q_pos + seed_len - 1) + ext.r_q;
    let final_t_start = t_start_idx - ext.r_t;
    let final_t_end = (t_start_idx + seed_len - 1) + ext.l_t;

    // Output coordinates (1-based, strand-aware for target)
    let original_len = ctx.index.get_sequence_len(candidate.target_idx);
    let (out_t_start, out_t_end, strand_char) = match candidate.strand {
        Strand::Reverse => {
            let fwd_start = original_len - 1 - final_t_end;
            let fwd_end = original_len - 1 - final_t_start;
            (fwd_start + 1, fwd_end + 1, '-')
        }
        Strand::Forward => (final_t_start + 1, final_t_end + 1, '+'),
    };

    let target_id = ctx.index.get_name(t_idx as u32);
    let ctx_len = 20;

    // Flanking Base sequences for binding-site output
    let (flank_5_range, flank_5_rev) = if format == OutputFormat::BindingSite {
        let start = final_t_start.saturating_sub(ctx_len);
        (start..final_t_start, true)
    } else {
        (0..0, false)
    };

    let (flank_3_range, flank_3_rev) = if format == OutputFormat::BindingSite {
        let t_3_end = (final_t_end + 1 + ctx_len).min(t_seq.len());
        if final_t_end + 1 < t_seq.len() {
            ((final_t_end + 1)..t_3_end, false)
        } else {
            (0..0, false)
        }
    } else {
        (0..0, false)
    };

    let mut itoa_buf = itoa::Buffer::new();
    let mut zmij_buf = zmij::Buffer::new();
    fill_line_buf(
        line_buf,
        &mut itoa_buf,
        &mut zmij_buf,
        format,
        queries.get_name(query_idx),
        final_q_start + 1,
        final_q_end + 1,
        target_id,
        out_t_start,
        out_t_end,
        strand_char,
        score,
        &ext.alignment,
        (t_seq, flank_5_range, flank_5_rev),
        (t_seq, flank_3_range, flank_3_rev),
        None,
    );

    if writer.write_all(line_buf).is_err() {
        return false;
    }

    true
}
