use anyhow::{Result, bail};
use log::info;
use rayon::prelude::*;
use std::io::Write;

use crate::config::{OutputCompression, OutputFormat, SearchArgs};
use crate::io::output::{compress_bytes, resolve_compression};
use crate::seed::{SeedCandidate, find_seeds_into};
use crate::types::Strand;

use super::core::extend_seed;
use super::output::fill_line_buf;
use super::{
    FilterReason, SaIndex, SearchContext, THREAD_EXTENDER, THREAD_MATCHES, THREAD_SEEDS,
};

/// Run search with streaming output - writes hits directly instead of collecting.
/// This avoids memory overhead for large result sets.
pub fn run_search_streaming<W: std::io::Write>(
    queries: &[(String, Vec<u8>)],
    index: &SaIndex<'_>,
    opts: &SearchArgs,
    writer: &mut W,
) -> Result<usize> {
    for (q_id, q_seq) in queries {
        if let Err(err) = opts.seed.seed.normalize(q_seq.len()) {
            bail!("Invalid seed spec for query '{}': {}", q_id, err);
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
    let format = opts
        .output
        .format
        .clone()
        .unwrap_or(OutputFormat::Detailed);
    let compression = opts.output.compress.unwrap_or(OutputCompression::None);
    let compression_cfg = resolve_compression(compression, opts.output.level)?;
    let chunk_size = std::cmp::max(1, rayon::current_num_threads() * 4);

    let process_query = |q_id: &str,
                         q_seq: &[u8],
                         writer: &mut dyn Write,
                         line_buf: &mut Vec<u8>|
     -> usize {
        let mut local_hits = 0usize;
        THREAD_EXTENDER.with(|ext| {
            THREAD_SEEDS.with(|seeds_cell| {
                THREAD_MATCHES.with(|matches_cell| {
                    let mut extender = ext.borrow_mut();
                    let mut seeds = seeds_cell.borrow_mut();
                    let mut matches = matches_cell.borrow_mut();
                    let mut ctx = SearchContext::with_extender(index, opts, &mut *extender);

                    // Reuse thread-local Vecs instead of allocating new ones
                    find_seeds_into(q_seq, ctx.index.index, &ctx.args.seed, &mut seeds, &mut matches);

                    for candidate in seeds.iter() {
                        if process_candidate_streaming(
                            crate::types::Query::new(q_id, q_seq),
                            candidate,
                            &mut ctx,
                            writer,
                            format,
                            line_buf,
                        ) {
                            local_hits += 1;
                        }
                    }
                });
            });
        });
        local_hits
    };

    for chunk in queries.chunks(chunk_size) {
        let outputs: Vec<Result<(Vec<u8>, usize)>> = chunk
            .par_iter()
            .map(|(q_id, q_seq)| {
                let mut line_buf: Vec<u8> = Vec::new();
                let mut out: Vec<u8> = Vec::new();
                let hits = process_query(q_id, q_seq, &mut out, &mut line_buf);
                let out = compress_bytes(compression_cfg, out)?;
                Ok((out, hits))
            })
            .collect();

        for item in outputs {
            let (buf, count) = item?;
            if !buf.is_empty() {
                writer.write_all(&buf)?;
            }
            hit_count += count;
        }
    }

    info!("Streaming search complete: {} hits", hit_count);
    Ok(hit_count)
}

/// Process candidate and write directly to output - zero String allocation variant.
/// Returns true if a hit was written.
fn process_candidate_streaming<W: Write + ?Sized>(
    query: crate::types::Query<'_>,
    candidate: &SeedCandidate,
    ctx: &mut SearchContext<'_, '_>,
    writer: &mut W,
    format: OutputFormat,
    line_buf: &mut Vec<u8>,
) -> bool {
    let t_idx = candidate.target_idx;
    let t_seq: &[u8] = match candidate.strand {
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
    let Some(ext) = extend_seed(ctx, query.seq, t_seq, candidate) else {
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

    let target_id = ctx.index.get_id(t_idx);
    let ctx_len = 20;
    let (flank_5, flank_5_rev, flank_3, flank_3_rev) = if format == OutputFormat::BindingSite {
        let flank_5_slice = &t_seq[final_t_start.saturating_sub(ctx_len)..final_t_start];
        let t_3_end = (final_t_end + 1 + ctx_len).min(t_seq.len());
        let flank_3_slice = if final_t_end + 1 < t_seq.len() {
            &t_seq[final_t_end + 1..t_3_end]
        } else {
            &t_seq[0..0]
        };
        (flank_5_slice, true, flank_3_slice, false)
    } else {
        (&t_seq[0..0], false, &t_seq[0..0], false)
    };

    let mut itoa_buf = itoa::Buffer::new();
    let mut zmij_buf = zmij::Buffer::new();
    fill_line_buf(
        line_buf,
        &mut itoa_buf,
        &mut zmij_buf,
        format,
        query.id,
        final_q_start + 1,
        final_q_end + 1,
        target_id,
        out_t_start,
        out_t_end,
        strand_char,
        score,
        &ext.alignment,
        (flank_5, flank_5_rev),
        (flank_3, flank_3_rev),
        Some(50),
    );

    if writer.write_all(line_buf).is_err() {
        return false;
    }

    true
}
