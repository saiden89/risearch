//! Search module - finds miRNA-target interactions.
//!
//! Flow: queries → seed finding → DP extension → hits

use anyhow::{Context, Result};
use log::{info, trace};
use rayon::prelude::*;
use smallvec::SmallVec;

use crate::alignment::{Alignment, Pairing};
use crate::config::{ExtendConfig, Matrix, OutputFormat, SearchArgs};
use crate::dp::{DpExtender, DpView};
use crate::dsm::{pair_mat, seed_energy, terminal_3p, terminal_5p, DsmModel, T04, T99};
use crate::index::store::{TargetStore, TargetView};
use crate::registry::{QueryRegistry, TargetRegistry};
use crate::seed::{for_each_seed_one_target, SeedHit, TargetSeedView};
use crate::seq::Sequence;
use crate::types::{Base, Energy, Strand};

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
    info!(
        "Starting search: {} queries, seed={:?}, max_ext={}, delta_g={}",
        queries.len(),
        opts.seed.seed,
        opts.extend.max_extension,
        opts.extend.delta_g
    );

    let count = match opts.extend.matrix {
        Matrix::T04 => stream_hits::<T04, W>(queries, index, opts, writer)?,
        Matrix::T99 => stream_hits::<T99, W>(queries, index, opts, writer)?,
    };

    info!("Search complete: {} hits", count);
    Ok(count)
}

/// Run search against mmap-backed target store and stream results directly to writer.
pub fn run_search_streaming_store<W: std::io::Write>(
    queries: &QueryRegistry,
    store: &TargetStore,
    opts: &SearchArgs,
    writer: &mut W,
) -> Result<usize> {
    info!(
        "Starting search: {} queries x {} targets, seed={:?}, max_ext={}, delta_g={}",
        queries.len(),
        store.len(),
        opts.seed.seed,
        opts.extend.max_extension,
        opts.extend.delta_g
    );

    let count = match opts.extend.matrix {
        Matrix::T04 => stream_hits_store::<T04, W>(queries, store, opts, writer)?,
        Matrix::T99 => stream_hits_store::<T99, W>(queries, store, opts, writer)?,
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
    out_buf: Vec<u8>,
    fmt_bufs: crate::output::OutputBuffers,
}

impl<M: DsmModel> SearchState<M> {
    fn new(penalty: i32) -> Self {
        Self {
            extender: DpExtender::<M>::with_penalty(penalty),
            penalty,
            seeds: Vec::with_capacity(128_000),
            out_buf: Vec::with_capacity(64 * 1024),
            fmt_bufs: crate::output::OutputBuffers::new(),
        }
    }

    fn clear(&mut self) {
        self.seeds.clear();
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
        .map_init(
            || SearchState::<M>::new(penalty),
            |state, (i, q)| {
                let mut hits = Vec::new();
                process_query::<M, _>(i as u32, q, state, index, opts, |hit| {
                    hits.push(hit);
                });
                hits
            },
        )
        .collect();

    Ok(hits.into_iter().flatten().collect())
}

fn stream_hits<M: DsmModel, W: std::io::Write>(
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
        let mut write_err: Option<std::io::Error> = None;

        process_query::<M, _>(q_idx as u32, q, &mut state, index, opts, |hit| {
            if write_err.is_some() {
                return;
            }

            if let Err(err) =
                crate::output::write_hit(&mut fmt_bufs, &hit, format, &mut out_buf, queries, index)
            {
                write_err = Some(err);
                return;
            }

            hit_count += 1;

            // Flush periodically
            if out_buf.len() >= 64 * 1024 {
                if let Err(err) = writer.write_all(&out_buf) {
                    write_err = Some(err);
                    return;
                }
                out_buf.clear();
            }
        });

        if let Some(err) = write_err {
            return Err(err.into());
        }

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

fn stream_hits_store<M: DsmModel, W: std::io::Write>(
    queries: &QueryRegistry,
    store: &TargetStore,
    opts: &SearchArgs,
    writer: &mut W,
) -> Result<usize> {
    let format = opts.output.format;
    let mut hit_count = 0;
    let penalty = (opts.extend.penalty * 100.0).round() as i32;
    let mut state = SearchState::<M>::new(penalty);

    for (target_idx, _) in store.iter_meta() {
        let target = store
            .target_view(target_idx as usize)
            .with_context(|| format!("Failed to load target #{}", target_idx))?;

        for (q_idx, q) in queries.entries().iter().enumerate() {
            // Move buffers out to avoid borrow conflicts with closure.
            let mut out_buf = std::mem::take(&mut state.out_buf);
            let mut fmt_bufs = std::mem::take(&mut state.fmt_bufs);
            let mut write_err: Option<std::io::Error> = None;

            process_query_one_target::<M, _>(
                q_idx as u32,
                q,
                &target,
                target_idx,
                &mut state,
                opts,
                |hit| {
                    if write_err.is_some() {
                        return;
                    }

                    let q_name = queries.get_name(hit.query_idx);
                    let t_name = target.name;
                    if let Err(err) = crate::output::write_hit_names(
                        &mut fmt_bufs,
                        &hit,
                        format,
                        &mut out_buf,
                        q_name,
                        t_name,
                    ) {
                        write_err = Some(err);
                        return;
                    }

                    hit_count += 1;

                    if out_buf.len() >= 64 * 1024 {
                        if let Err(err) = writer.write_all(&out_buf) {
                            write_err = Some(err);
                            return;
                        }
                        out_buf.clear();
                    }
                },
            );

            if let Some(err) = write_err {
                return Err(err.into());
            }

            if !out_buf.is_empty() {
                writer.write_all(&out_buf)?;
                out_buf.clear();
            }
            state.out_buf = out_buf;
            state.fmt_bufs = fmt_bufs;
        }
    }

    Ok(hit_count)
}

/// Core query processing: find seeds → extend → filter → emit hits.
struct QueryCtx<'a> {
    q_idx: u32,
    q_seq: &'a [Base],
    interval: crate::types::Interval,
    include_alignment: bool,
    pair_matrix: &'static [[u8; 6]; 6],
    delta_g: f64,
    extend_cfg: &'a ExtendConfig,
}

fn emit_seed_hit<M: DsmModel, F: FnMut(SearchHit)>(
    extender: &mut DpExtender<M>,
    penalty: i32,
    ctx: &QueryCtx<'_>,
    seed: &SeedHit,
    t_seq: &[Base],
    original_len: usize,
    on_hit: &mut F,
) {
    if seed.target_start + seed.seed_len.get() > t_seq.len() {
        return;
    }

    let Some(ext) = extend_seed::<M>(
        extender,
        penalty,
        ctx.q_seq,
        t_seq,
        seed,
        ctx.interval,
        ctx.extend_cfg,
        ctx.pair_matrix,
        ctx.include_alignment,
    ) else {
        return;
    };

    if ext.score > ctx.delta_g {
        return;
    }

    on_hit(SearchHit::new(
        ctx.q_idx,
        ctx.q_seq,
        t_seq,
        seed,
        ext,
        ctx.include_alignment,
        original_len,
    ));
}

fn process_query<M: DsmModel, F: FnMut(SearchHit)>(
    q_idx: u32,
    q: &crate::registry::QueryData,
    state: &mut SearchState<M>,
    index: &TargetRegistry,
    opts: &SearchArgs,
    mut on_hit: F,
) {
    state.clear();
    crate::seed::find_seeds(q, index, &opts.seed, &mut state.seeds);
    let ctx = QueryCtx {
        q_idx,
        q_seq: q.sequence(),
        interval: q.seed_interval(),
        include_alignment: opts.output.format != OutputFormat::Minimal,
        pair_matrix: pair_mat(opts.seed.allows_wobble()),
        delta_g: opts.extend.delta_g,
        extend_cfg: &opts.extend,
    };

    for seed in state.seeds.iter() {
        let t_seq = match seed.strand {
            Strand::Reverse => index.get_sequence_rc(seed.target_id.0 as usize),
            Strand::Forward => index.get_sequence(seed.target_id.0 as usize),
        };
        emit_seed_hit::<M, _>(
            &mut state.extender,
            state.penalty,
            &ctx,
            seed,
            t_seq,
            index.get_sequence_len(seed.target_id.0 as usize),
            &mut on_hit,
        );
    }
}

fn process_query_one_target<M: DsmModel, F: FnMut(SearchHit)>(
    q_idx: u32,
    q: &crate::registry::QueryData,
    target: &TargetView<'_>,
    target_idx: u32,
    state: &mut SearchState<M>,
    opts: &SearchArgs,
    mut on_hit: F,
) {
    let target_seed_view = TargetSeedView {
        combined_seq: target.combined_seq,
        combined_sa: target.combined_sa,
        seq_len: target.seq_len,
    };
    let ctx = QueryCtx {
        q_idx,
        q_seq: q.sequence(),
        interval: q.seed_interval(),
        include_alignment: opts.output.format != OutputFormat::Minimal,
        pair_matrix: pair_mat(opts.seed.allows_wobble()),
        delta_g: opts.extend.delta_g,
        extend_cfg: &opts.extend,
    };
    let seq_len = target.seq_len;
    let t_fwd = &target.combined_seq[..seq_len];
    let t_rc = &target.combined_seq[seq_len + 1..2 * seq_len + 1];
    for_each_seed_one_target(q, target_idx, &target_seed_view, &opts.seed, |seed| {
        let t_seq = match seed.strand {
            Strand::Reverse => t_rc,
            Strand::Forward => t_fwd,
        };
        emit_seed_hit::<M, _>(
            &mut state.extender,
            state.penalty,
            &ctx,
            &seed,
            t_seq,
            seq_len,
            &mut on_hit,
        );
    });
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
    left_pairs: SmallVec<[Pairing; 64]>,
    right_pairs: SmallVec<[Pairing; 64]>,
}

fn extend_seed<M: DsmModel>(
    extender: &mut DpExtender<M>,
    penalty: i32,
    q_seq: &[Base],
    t_seq: &[Base],
    seed: &SeedHit,
    interval: crate::types::Interval,
    extend_cfg: &ExtendConfig,
    pair_matrix: &'static [[u8; 6]; 6],
    with_traceback: bool,
) -> Option<Extension> {
    let q_pos = seed.query_pos;
    let t_pos = seed.target_start;
    let len = seed.seed_len.get();

    // Maximality check: left
    if q_pos > interval.start && t_pos + len < t_seq.len() {
        let p_class = pair_matrix[q_seq[q_pos - 1].idx()][t_seq[t_pos + len].idx()];
        if p_class != 0 && !extend_cfg.no_max_prune {
            trace!("Filtered: non-maximal left");
            return None;
        }
    }

    // Maximality check: right
    if q_pos + len < interval.end && t_pos > 0 {
        let p_class = pair_matrix[q_seq[q_pos + len].idx()][t_seq[t_pos - 1].idx()];
        if p_class != 0 && !extend_cfg.no_max_prune {
            trace!("Filtered: non-maximal right");
            return None;
        }
    }

    let t_match_end = t_pos + len - 1;
    let max_ext = (extend_cfg.max_extension as usize).min(MAX_DP_EXT);
    // Seed energy
    let seed_energy = seed_energy::<M>(q_seq, t_seq, q_pos, t_match_end, len, penalty);

    // Check if extension is possible
    let can_extend_left = q_pos > 0 && t_pos + len < t_seq.len();
    let can_extend_right = q_pos + len < q_seq.len() && t_pos > 0;

    if max_ext == 0 || (!can_extend_left && !can_extend_right) {
        // Seed only - add terminal penalties
        let term_5p = terminal_5p::<M>(q_seq[q_pos], t_seq[t_match_end], penalty);
        let term_3p = terminal_3p::<M>(q_seq[q_pos + len - 1], t_seq[t_pos], penalty);
        let nt_count = (2 * len) as i32;
        return Some(Extension {
            score: M::to_kcal(seed_energy + term_5p + term_3p + nt_count * penalty),
            l_q: 0,
            l_t: 0,
            r_q: 0,
            r_t: 0,
            left_pairs: SmallVec::new(),
            right_pairs: SmallVec::new(),
        });
    }

    // DP extension: left (scope block releases borrow on extender)
    let (l_score, l_q, l_t, left_pairs) = {
        let view = DpView::<M>::left(q_seq, t_seq, q_pos, t_match_end, max_ext);
        let result = extender.extend(&view);
        let pairs = if with_traceback {
            result.traceback(&view)
        } else {
            SmallVec::new()
        };
        (result.score, result.q_len, result.t_len, pairs)
    };

    // DP extension: right
    let (r_score, r_q, r_t, right_pairs) = {
        let view = DpView::<M>::right(q_seq, t_seq, q_pos + len - 1, t_pos, max_ext);
        let result = extender.extend(&view);
        let pairs = if with_traceback {
            result.traceback(&view)
        } else {
            SmallVec::new()
        };
        (result.score, result.q_len, result.t_len, pairs)
    };

    let nt_count = (l_q + l_t + r_q + r_t + 2 * len) as i32;

    Some(Extension {
        score: M::to_kcal(seed_energy + l_score + r_score + nt_count * penalty),
        l_q,
        l_t,
        r_q,
        r_t,
        left_pairs,
        right_pairs,
    })
}

// =============================================================================
// HIT BUILDING
// =============================================================================

fn build_seed_pairs(
    q_seq: &[Base],
    t_seq: &[Base],
    q_pos: usize,
    t_match_end: usize,
    len: usize,
) -> SmallVec<[Pairing; 64]> {
    let mut pairs = SmallVec::with_capacity(len);
    for i in 0..len {
        pairs.push(Pairing::from_bases(
            q_seq[q_pos + i],
            t_seq[t_match_end - i],
        ));
    }
    pairs
}

impl SearchHit {
    fn new(
        query_idx: u32,
        q_seq: &[Base],
        t_seq: &[Base],
        seed: &SeedHit,
        ext: Extension,
        include_alignment: bool,
        original_len: usize,
    ) -> Self {
        let Extension {
            score,
            l_q,
            l_t,
            r_q,
            r_t,
            left_pairs,
            right_pairs,
        } = ext;

        let q_pos = seed.query_pos;
        let t_start = seed.target_start;
        let len = seed.seed_len.get();

        let final_q_start = q_pos.saturating_sub(l_q);
        let final_q_end = (q_pos + len - 1) + r_q;
        let final_t_start = t_start.saturating_sub(r_t);
        let final_t_end = (t_start + len - 1) + l_t;

        let (out_t_start, out_t_end, strand) = match seed.strand {
            Strand::Reverse => {
                let fwd_start = original_len - 1 - final_t_end;
                let fwd_end = original_len - 1 - final_t_start;
                (fwd_start + 1, fwd_end + 1, Strand::Reverse)
            }
            Strand::Forward => (final_t_start + 1, final_t_end + 1, Strand::Forward),
        };

        let alignment = if include_alignment {
            let t_match_end = seed.target_start + len - 1;
            let seed_pairs = build_seed_pairs(q_seq, t_seq, q_pos, t_match_end, len);
            Some(Alignment::new(&left_pairs, &seed_pairs, &right_pairs))
        } else {
            None
        };

        Self {
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
            strand,
            energy: score.into(),
            alignment,
            flank_5: Sequence::from(Vec::new()),
            flank_3: Sequence::from(Vec::new()),
        }
    }
}
