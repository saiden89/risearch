//! Unified search pipeline for processing queries.
//!
//! This module contains the shared query processing logic used by both
//! in-memory and streaming search paths.

use anyhow::Result;
use rayon::prelude::*;

use crate::config::OutputFormat;
use crate::dp;
use crate::dsm::DsmModel;
use crate::registry::QueryRegistry;
use crate::sa::TargetRegistry;
use crate::seed::SeedHit;
use crate::seq::Sequence;
use crate::types::Strand;

use super::extend::{ExtensionResult, build_alignment_from_extension, extend_seed};
use super::{SearchContext, SearchHit, filter_by_energy, validate_seed_bounds};

// =============================================================================
// PER-QUERY PROCESSING STATE
// =============================================================================

/// Per-query processing state (reusable across queries for efficiency).
pub(super) struct QueryProcessingState<M: DsmModel> {
    pub extender: dp::DpExtender<M>,
    pub seeds: Vec<SeedHit>,
    pub matches: Vec<crate::seed::SeedMatch>,
}

impl<M: DsmModel> QueryProcessingState<M> {
    pub fn new() -> Self {
        Self {
            extender: dp::DpExtender::<M>::new(),
            seeds: Vec::with_capacity(128_000),
            matches: Vec::with_capacity(1024),
        }
    }

    pub fn reset(&mut self) {
        self.seeds.clear();
        self.matches.clear();
    }
}

impl<M: DsmModel> Default for QueryProcessingState<M> {
    fn default() -> Self {
        Self::new()
    }
}

// =============================================================================
// SHARED QUERY PROCESSING
// =============================================================================

/// Process a single query and return all valid hits.
///
/// This is the shared core logic used by both in-memory and streaming paths.
/// All seeding, extension, and filtering happens here.
pub(super) fn process_query_hits<M: DsmModel>(
    q_idx: u32,
    q: &crate::registry::QueryData,
    state: &mut QueryProcessingState<M>,
    queries: &QueryRegistry,
    index: &TargetRegistry,
    opts: &crate::config::SearchArgs,
) -> Vec<SearchHit> {
    state.reset();

    let mut ctx = SearchContext::<M>::with_extender(index, opts, &mut state.extender);
    crate::seed::find_seeds(
        q,
        ctx.index,
        &ctx.args.seed,
        &mut state.seeds,
        &mut state.matches,
    );

    state
        .seeds
        .iter()
        .filter_map(|candidate| {
            // Retrieve target sequence based on strand
            let t_seq = match candidate.strand {
                Strand::Reverse => ctx.index.get_sequence_rc(candidate.target_idx),
                Strand::Forward => ctx.index.get_sequence(candidate.target_idx),
            };

            // PRE-EXTENSION VALIDATION: Validate seed bounds
            if !validate_seed_bounds(candidate, t_seq, &mut ctx) {
                return None;
            }

            // EXTENSION: Extend seed and build hit
            let hit = extend_and_build_hit(q_idx, queries, candidate, t_seq, &mut ctx)?;

            // POST-EXTENSION FILTERING: Energy threshold
            if !filter_by_energy(hit.energy.as_f64(), &mut ctx) {
                return None;
            }

            Some(hit)
        })
        .collect()
}

// =============================================================================
// HIT BUILDING
// =============================================================================

/// Extend a seed candidate and build a `SearchHit` structure.
fn extend_and_build_hit<M: DsmModel>(
    query_idx: u32,
    queries: &QueryRegistry,
    candidate: &SeedHit,
    t_seq: &Sequence,
    ctx: &mut SearchContext<'_, '_, M>,
) -> Option<SearchHit> {
    let query_seq = queries.get(query_idx).sequence();
    let extension = extend_seed(ctx, query_seq, t_seq, candidate)?;

    Some(build_hit::<M>(
        query_idx, ctx, candidate, extension, query_seq, t_seq,
    ))
}

/// Finalize coordinates and build SearchHit from extension result.
fn build_hit<M: DsmModel>(
    query_idx: u32,
    ctx: &SearchContext<'_, '_, M>,
    candidate: &SeedHit,
    ext: ExtensionResult,
    query_seq: &Sequence,
    t_seq: &Sequence,
) -> SearchHit {
    let seed_len = candidate.len;
    let q_pos = candidate.query_pos;
    let t_start_idx = candidate.target_start;

    let final_q_start = q_pos.saturating_sub(ext.l_q);
    let final_q_end = (q_pos + seed_len - 1) + ext.r_q;
    let final_t_start = t_start_idx.saturating_sub(ext.r_t);
    let final_t_end = (t_start_idx + seed_len - 1) + ext.l_t;

    let original_len = ctx.index.get_sequence_len(candidate.target_idx);
    let (out_t_start, out_t_end, strand_char) = match candidate.strand {
        Strand::Reverse => {
            let fwd_start = original_len - 1 - final_t_end;
            let fwd_end = original_len - 1 - final_t_start;
            (fwd_start + 1, fwd_end + 1, '-')
        }
        Strand::Forward => (final_t_start + 1, final_t_end + 1, '+'),
    };

    // Build alignment only when output format requires it
    let alignment = match ctx.args.output.format {
        OutputFormat::Minimal => None,
        _ => Some(build_alignment_from_extension(
            &ext,
            query_seq,
            t_seq,
            t_start_idx,
        )),
    };

    SearchHit {
        query_idx,
        target_idx: candidate.target_idx as u32,
        q_start: final_q_start,
        q_end: final_q_end,
        t_start: final_t_start,
        t_end: final_t_end,
        output_q_start: final_q_start + 1,
        output_q_end: final_q_end + 1,
        output_t_start: out_t_start,
        output_t_end: out_t_end,
        strand: strand_char.into(),
        energy: ext.score.into(),
        alignment,
        flank_5: Sequence::from(Vec::new()),
        flank_3: Sequence::from(Vec::new()),
    }
}

// =============================================================================
// IN-MEMORY SEARCH ENTRY POINT
// =============================================================================

/// Run search and collect all hits into a Vec (in-memory mode).
pub(super) fn run_search_impl<M: DsmModel>(
    queries: &QueryRegistry,
    index: &TargetRegistry,
    opts: &crate::config::SearchArgs,
) -> Result<Vec<SearchHit>> {
    let per_query_hits: Vec<Vec<SearchHit>> = queries
        .entries()
        .par_iter()
        .enumerate()
        .map_init(QueryProcessingState::<M>::new, |state, (i, q)| {
            process_query_hits(i as u32, q, state, queries, index, opts)
        })
        .collect();

    Ok(per_query_hits.into_iter().flatten().collect())
}
