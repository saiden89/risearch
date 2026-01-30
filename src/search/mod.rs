//! Search module for miRNA target finding.
//!
//! ## Module Structure
//! - `mod.rs` - Entry points and shared types
//! - `pipeline.rs` - Unified query processing logic  
//! - `extend.rs` - DP extension algorithm

use anyhow::{Result, bail};
use log::info;
use rayon::prelude::*;
use std::collections::HashMap;

use crate::alignment::Alignment;
use crate::config::{Matrix, SearchArgs};
use crate::dp;
use crate::dsm::{DsmModel, T04, T99};
use crate::registry::QueryRegistry;
use crate::sa::TargetRegistry;
use crate::seq::Sequence;
use crate::types::{Energy, Strand};

mod extend;
mod pipeline;

use pipeline::{QueryProcessingState, process_query_hits};

const MAX_DP_EXT: usize = 50;

// =============================================================================
// PUBLIC ENTRY POINTS
// =============================================================================

/// Run search and return all hits (in-memory mode).
pub fn run_search(
    queries: &QueryRegistry,
    index: &TargetRegistry,
    opts: &SearchArgs,
) -> Result<Vec<SearchHit>> {
    match opts.extend.matrix {
        Matrix::T04 => pipeline::run_search_impl::<T04>(queries, index, opts),
        Matrix::T99 => pipeline::run_search_impl::<T99>(queries, index, opts),
    }
}

/// Run search and stream results directly to writer (streaming mode).
pub fn run_search_streaming<W: std::io::Write>(
    queries: &QueryRegistry,
    index: &TargetRegistry,
    opts: &SearchArgs,
    writer: &mut W,
) -> Result<usize> {
    match opts.extend.matrix {
        Matrix::T04 => run_search_streaming_impl::<T04, W>(queries, index, opts, writer),
        Matrix::T99 => run_search_streaming_impl::<T99, W>(queries, index, opts, writer),
    }
}

/// Streaming search implementation using shared pipeline.
fn run_search_streaming_impl<M: DsmModel, W: std::io::Write>(
    queries: &QueryRegistry,
    index: &TargetRegistry,
    opts: &SearchArgs,
    writer: &mut W,
) -> Result<usize> {
    // Validate seed specs upfront
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

    let format = opts.output.format;
    let chunk_size = std::cmp::max(1, rayon::current_num_threads() * 4);

    let mut hit_count = 0;
    let mut chunk_start = 0usize;

    for chunk in queries.entries().chunks(chunk_size) {
        let outputs: Vec<(Vec<u8>, usize)> = chunk
            .par_iter()
            .enumerate()
            .map_init(QueryProcessingState::<M>::new, |state, (i, q)| {
                let q_idx = (chunk_start + i) as u32;
                let hits = process_query_hits(q_idx, q, state, queries, index, opts);

                let mut out_buf = Vec::with_capacity(4096);
                for hit in &hits {
                    hit.write_with_format(&mut out_buf, format, queries, index)
                        .ok();
                }
                (out_buf, hits.len())
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

// =============================================================================
// TYPES
// =============================================================================

/// High-level algorithm stages for structured logging
#[derive(Debug, Clone, Copy)]
enum SearchStage {
    Extend,
}

impl std::fmt::Display for SearchStage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Extend => write!(f, "[EXTEND]"),
        }
    }
}

/// Reasons why a seed, candidate, or hit was filtered out
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FilterReason {
    SeedContainsN,
    SeedOutOfBounds,
    MaximalityLeft,
    MaximalityRight,
    EnergyAboveThreshold,
}

/// Statistics for search filtering
#[derive(Debug, Default)]
pub struct SearchStats {
    pub seeds_tried: usize,
    pub candidates_processed: usize,
    pub hits_final: usize,
    pub filtered: HashMap<FilterReason, usize>,
}

impl SearchStats {
    pub fn record_filter(&mut self, reason: FilterReason) {
        *self.filtered.entry(reason).or_insert(0) += 1;
    }
}

/// A search hit representing a miRNA-target interaction.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SearchHit {
    pub query_idx: u32,
    pub target_idx: u32,
    pub q_start: usize,
    pub q_end: usize,
    pub t_start: usize,
    pub t_end: usize,
    pub output_q_start: usize,
    pub output_q_end: usize,
    pub output_t_start: usize,
    pub output_t_end: usize,
    pub strand: Strand,
    pub energy: Energy,
    pub alignment: Option<Alignment>,
    pub flank_5: Sequence,
    pub flank_3: Sequence,
}

/// Search context holding index, args, and per-query state.
pub struct SearchContext<'a, 'e, M: DsmModel> {
    pub index: &'a TargetRegistry,
    pub args: &'a SearchArgs,
    pub extender: &'e mut dp::DpExtender<M>,
    pub stats: SearchStats,
}

impl<'a, 'e, M: DsmModel> SearchContext<'a, 'e, M> {
    pub fn with_extender(
        index: &'a TargetRegistry,
        args: &'a SearchArgs,
        extender: &'e mut dp::DpExtender<M>,
    ) -> Self {
        Self {
            index,
            args,
            extender,
            stats: SearchStats::default(),
        }
    }
}

// =============================================================================
// SHARED FILTER FUNCTIONS
// =============================================================================

/// Validate that the seed hit is within target sequence bounds.
#[inline]
fn validate_seed_bounds<M: DsmModel>(
    candidate: &crate::seed::SeedHit,
    t_seq: &Sequence,
    ctx: &mut SearchContext<'_, '_, M>,
) -> bool {
    let t_start_idx = candidate.target_start;
    let seed_len = candidate.len;

    if t_start_idx + seed_len > t_seq.len() {
        ctx.stats.record_filter(FilterReason::SeedOutOfBounds);
        return false;
    }
    true
}

/// Filter by energy threshold.
#[inline]
fn filter_by_energy<M: DsmModel>(score: f64, ctx: &mut SearchContext<'_, '_, M>) -> bool {
    if score > ctx.args.extend.delta_g {
        ctx.stats.record_filter(FilterReason::EnergyAboveThreshold);
        return false;
    }
    true
}
