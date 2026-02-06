//! Search module - finds miRNA-target interactions.
//!
//! Flow: queries → seed finding → DP extension → hits

use anyhow::{bail, Result};
use log::{info, trace};
use rayon::prelude::*;
use smallvec::SmallVec;

use crate::alignment::{Alignment, Pairing};
use crate::config::{Matrix, OutputFormat, SearchArgs};
use crate::dp::{DpExtender, DpView};
use crate::dsm::{
    pair_mat, seed_energy_with_penalty, terminal_3p_with_penalty, terminal_5p_with_penalty,
    DsmModel, T04, T99,
};
use crate::registry::QueryRegistry;
use crate::sa::TargetRegistry;
use crate::seed::SeedHit;
use crate::seq::Sequence;
use crate::types::{Energy, Strand};

const MAX_DP_EXT: usize = 50;

// =============================================================================
// PUBLIC API
// =============================================================================

/// A search hit representing a miRNA-target interaction.
#[derive(Debug, Clone)]
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

/// Run search and return all hits (in-memory, parallel).
pub fn run_search(
    queries: &QueryRegistry,
    index: &TargetRegistry,
    opts: &SearchArgs,
) -> Result<Vec<SearchHit>> {
    match opts.extend.matrix {
        Matrix::T04 => run_search_parallel::<T04>(queries, index, opts),
        Matrix::T99 => run_search_parallel::<T99>(queries, index, opts),
    }
}

/// Run search and stream results directly to writer.
pub fn run_search_streaming<W: std::io::Write>(
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
        "Starting search: {} queries, seed={:?}, max_ext={}, delta_g={}",
        queries.len(),
        opts.seed.seed,
        opts.extend.max_extension,
        opts.extend.delta_g
    );

    let count = match opts.extend.matrix {
        Matrix::T04 => run_search_streaming_impl::<T04, W>(queries, index, opts, writer)?,
        Matrix::T99 => run_search_streaming_impl::<T99, W>(queries, index, opts, writer)?,
    };

    info!("Search complete: {} hits", count);
    Ok(count)
}

// =============================================================================
// IMPLEMENTATION
// =============================================================================

/// Reusable state for processing queries.
struct SearchState<M: DsmModel> {
    extender: DpExtender<M>,
    penalty: i32,
    seeds: Vec<SeedHit>,
    matches: Vec<crate::seed::SeedMatch>,
    out_buf: Vec<u8>,
    fmt_bufs: crate::output::OutputBuffers,
}

impl<M: DsmModel> SearchState<M> {
    fn new(penalty: i32) -> Self {
        Self {
            extender: DpExtender::<M>::with_penalty(penalty),
            penalty,
            seeds: Vec::with_capacity(128_000),
            matches: Vec::with_capacity(1024),
            out_buf: Vec::with_capacity(64 * 1024),
            fmt_bufs: crate::output::OutputBuffers::new(),
        }
    }

    fn clear(&mut self) {
        self.seeds.clear();
        self.matches.clear();
        self.out_buf.clear();
    }
}

fn run_search_parallel<M: DsmModel>(
    queries: &QueryRegistry,
    index: &TargetRegistry,
    opts: &SearchArgs,
) -> Result<Vec<SearchHit>> {
    let penalty = (opts.extend.penalty * 100.0).round() as i32;
    let hits: Vec<Vec<SearchHit>> = queries
        .entries()
        .par_iter()
        .enumerate()
        .map_init(|| SearchState::<M>::new(penalty), |state, (i, q)| {
            let mut hits = Vec::new();
            process_query::<M, _>(i as u32, q, state, queries, index, opts, |hit| {
                hits.push(hit);
            });
            hits
        })
        .collect();

    Ok(hits.into_iter().flatten().collect())
}

fn run_search_streaming_impl<M: DsmModel, W: std::io::Write>(
    queries: &QueryRegistry,
    index: &TargetRegistry,
    opts: &SearchArgs,
    writer: &mut W,
) -> Result<usize> {
    let format = opts.output.format;
    let mut hit_count = 0;
    let penalty = (opts.extend.penalty * 100.0).round() as i32;
    let mut state = SearchState::<M>::new(penalty);

    for (q_idx, q) in queries.entries().iter().enumerate() {
        // Move buffers out to avoid borrow conflicts with closure
        let mut out_buf = std::mem::take(&mut state.out_buf);
        let mut fmt_bufs = std::mem::take(&mut state.fmt_bufs);

        process_query::<M, _>(q_idx as u32, q, &mut state, queries, index, opts, |hit| {
            let _ = crate::output::write_hit_with_format(
                &mut fmt_bufs,
                &hit,
                format,
                &mut out_buf,
                queries,
                index,
            );
            hit_count += 1;

            // Flush periodically
            if out_buf.len() >= 64 * 1024 {
                let _ = writer.write_all(&out_buf);
                out_buf.clear();
            }
        });

        // Flush remaining and restore buffers
        if !out_buf.is_empty() {
            writer.write_all(&out_buf)?;
            out_buf.clear();
        }
        state.out_buf = out_buf;
        state.fmt_bufs = fmt_bufs;
    }

    Ok(hit_count)
}

/// Core query processing: find seeds → extend → filter → emit hits.
fn process_query<M: DsmModel, F: FnMut(SearchHit)>(
    q_idx: u32,
    q: &crate::registry::QueryData,
    state: &mut SearchState<M>,
    queries: &QueryRegistry,
    index: &TargetRegistry,
    opts: &SearchArgs,
    mut on_hit: F,
) {
    state.clear();

    // Find seeds
    crate::seed::find_seeds(q, index, &opts.seed, &mut state.seeds, &mut state.matches);

    let query = queries.get(q_idx);
    let q_seq = query.sequence();
    let interval = query.seed_interval();

    for seed in state.seeds.iter() {
        let t_seq = match seed.strand {
            Strand::Reverse => index.get_sequence_rc(seed.target_id.0 as usize),
            Strand::Forward => index.get_sequence(seed.target_id.0 as usize),
        };

        // Bounds check
        if seed.target_start + seed.seed_len.get() > t_seq.len() {
            continue;
        }

        // Extend seed
        let Some(ext) = extend_seed::<M>(
            &mut state.extender,
            state.penalty,
            q_seq,
            t_seq,
            seed,
            interval,
            opts,
        )
        else {
            continue;
        };

        // Energy filter
        if ext.score > opts.extend.delta_g {
            continue;
        }

        on_hit(build_hit(q_idx, index, seed, ext));
    }
}

// =============================================================================
// DP EXTENSION
// =============================================================================

struct Extension {
    score: f64,
    l_q: usize,
    l_t: usize,
    r_q: usize,
    r_t: usize,
    alignment: Option<Alignment>,
}

fn extend_seed<M: DsmModel>(
    extender: &mut DpExtender<M>,
    penalty: i32,
    q_seq: &Sequence,
    t_seq: &Sequence,
    seed: &SeedHit,
    interval: crate::types::Interval,
    opts: &SearchArgs,
) -> Option<Extension> {
    let q_pos = seed.query_pos;
    let t_pos = seed.target_start;
    let len = seed.seed_len.get();
    let pair_mat = pair_mat(opts.seed.allows_wobble());

    // Maximality check: left
    if q_pos > interval.start && t_pos + len < t_seq.len() {
        let p_class = pair_mat[q_seq[q_pos - 1].idx()][t_seq[t_pos + len].idx()];
        if p_class != 0 && !opts.extend.no_max_prune {
            trace!("Filtered: non-maximal left");
            return None;
        }
    }

    // Maximality check: right
    if q_pos + len < interval.end && t_pos > 0 {
        let p_class = pair_mat[q_seq[q_pos + len].idx()][t_seq[t_pos - 1].idx()];
        if p_class != 0 && !opts.extend.no_max_prune {
            trace!("Filtered: non-maximal right");
            return None;
        }
    }

    let t_match_end = t_pos + len - 1;
    let max_ext = (opts.extend.max_extension as usize).min(MAX_DP_EXT);
    let needs_alignment = opts.output.format != OutputFormat::Minimal;

    // Seed energy
    let seed_energy = seed_energy_with_penalty::<M>(q_seq, t_seq, q_pos, t_match_end, len, penalty);

    // Check if extension is possible
    let can_extend_left = q_pos > 0 && t_pos + len < t_seq.len();
    let can_extend_right = q_pos + len < q_seq.len() && t_pos > 0;

    if max_ext == 0 || (!can_extend_left && !can_extend_right) {
        // Seed only - add terminal penalties
        let term_5p = terminal_5p_with_penalty::<M>(q_seq[q_pos], t_seq[t_match_end], penalty);
        let term_3p = terminal_3p_with_penalty::<M>(q_seq[q_pos + len - 1], t_seq[t_pos], penalty);
        let nt_count = (2 * len) as i32;
        let alignment = if needs_alignment {
            let seed_pairs = build_seed_pairs(q_seq, t_seq, q_pos, t_match_end, len);
            Some(Alignment::new(&[], &seed_pairs, &[]))
        } else {
            None
        };
        return Some(Extension {
            score: M::to_kcal(seed_energy + term_5p + term_3p + nt_count * penalty),
            l_q: 0,
            l_t: 0,
            r_q: 0,
            r_t: 0,
            alignment,
        });
    }

    // DP extension: left (scope block releases borrow on extender)
    let (l_score, l_q, l_t, left_pairs) = {
        let view = DpView::<M>::left(q_seq, t_seq, q_pos, t_match_end, max_ext, penalty);
        let result = extender.extend(&view);
        let pairs = if needs_alignment { result.traceback(&view) } else { SmallVec::new() };
        (result.score, result.q_len, result.t_len, pairs)
    };

    // DP extension: right
    let (r_score, r_q, r_t, right_pairs) = {
        let view = DpView::<M>::right(q_seq, t_seq, q_pos + len - 1, t_pos, max_ext, penalty);
        let result = extender.extend(&view);
        let pairs = if needs_alignment { result.traceback(&view) } else { SmallVec::new() };
        (result.score, result.q_len, result.t_len, pairs)
    };

    let nt_count = (l_q + l_t + r_q + r_t + 2 * len) as i32;
    let alignment = if needs_alignment {
        let seed_pairs = build_seed_pairs(q_seq, t_seq, q_pos, t_match_end, len);
        Some(Alignment::new(&left_pairs, &seed_pairs, &right_pairs))
    } else {
        None
    };

    Some(Extension {
        score: M::to_kcal(seed_energy + l_score + r_score + nt_count * penalty),
        l_q,
        l_t,
        r_q,
        r_t,
        alignment,
    })
}

// =============================================================================
// HIT BUILDING
// =============================================================================

fn build_seed_pairs(
    q_seq: &Sequence,
    t_seq: &Sequence,
    q_pos: usize,
    t_match_end: usize,
    len: usize,
) -> SmallVec<[Pairing; 64]> {
    let mut pairs = SmallVec::with_capacity(len);
    for i in 0..len {
        pairs.push(Pairing::from_bases(q_seq[q_pos + i], t_seq[t_match_end - i]));
    }
    pairs
}

fn build_hit(
    query_idx: u32,
    index: &TargetRegistry,
    seed: &SeedHit,
    ext: Extension,
) -> SearchHit {
    let q_pos = seed.query_pos;
    let t_start = seed.target_start;
    let len = seed.seed_len.get();

    let final_q_start = q_pos.saturating_sub(ext.l_q);
    let final_q_end = (q_pos + len - 1) + ext.r_q;
    let final_t_start = t_start.saturating_sub(ext.r_t);
    let final_t_end = (t_start + len - 1) + ext.l_t;

    let original_len = index.get_sequence_len(seed.target_id.0 as usize);
    let (out_t_start, out_t_end, strand_char) = match seed.strand {
        Strand::Reverse => {
            let fwd_start = original_len - 1 - final_t_end;
            let fwd_end = original_len - 1 - final_t_start;
            (fwd_start + 1, fwd_end + 1, '-')
        }
        Strand::Forward => (final_t_start + 1, final_t_end + 1, '+'),
    };

    SearchHit {
        query_idx,
        target_idx: seed.target_id.0,
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
        alignment: ext.alignment,
        flank_5: Sequence::from(Vec::new()),
        flank_3: Sequence::from(Vec::new()),
    }
}
