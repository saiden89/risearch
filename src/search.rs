use anyhow::{Context, Result};
use log::{debug, info, trace, warn};
use std::io::Write;
use std::path::Path;

use crate::args::SearchArgs;
use crate::dp;
use crate::dsm::{EnergyModel, PAIR_MAT, PAIR_MAT_NO_GU};
use crate::sa::SaIndexFile;
use crate::seed::{SeedCandidate, build_seed_alignment, find_seeds_into};
use crate::parallel_sa::ParallelSeedMatch;
use crate::seq::Seq;
use crate::types::{Alignment, Energy, Pairing, QueryId, SeedPairing, Strand, TargetId};

use std::collections::HashMap;

const MAX_DP_EXT: usize = 50;

/// Convert bytes to RNA string (T->U), optionally reversed. Single-pass, one allocation.
#[inline]
fn bytes_to_rna_string(s: &[u8], reverse: bool) -> String {
    let mut result = String::with_capacity(s.len());
    if reverse {
        for &b in s.iter().rev() {
            result.push(match b {
                b'T' => 'U',
                b't' => 'u',
                _ => b as char,
            });
        }
    } else {
        for &b in s {
            result.push(match b {
                b'T' => 'U',
                b't' => 'u',
                _ => b as char,
            });
        }
    }
    result
}

/// Write bytes as RNA (T->U) directly to writer, optionally reversed. Zero allocation.
#[inline]
fn write_bytes_as_rna<W: Write>(w: &mut W, s: &[u8], reverse: bool) -> std::io::Result<()> {
    if reverse {
        for &b in s.iter().rev() {
            let c = match b {
                b'T' => b'U',
                b't' => b'u',
                _ => b,
            };
            w.write_all(&[c])?;
        }
    } else {
        for &b in s {
            let c = match b {
                b'T' => b'U',
                b't' => b'u',
                _ => b,
            };
            w.write_all(&[c])?;
        }
    }
    Ok(())
}

/// High-level algorithm stages for structured logging
#[derive(Debug, Clone, Copy)]
enum SearchStage {
    Input,  // Query parsing and processing
    Seed,   // Suffix array search, seed generation
    Extend, // DP extension, maximality checks
    Output, // Final results
}

impl std::fmt::Display for SearchStage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Input => write!(f, "[INPUT]"),
            Self::Seed => write!(f, "[SEED]"),
            Self::Extend => write!(f, "[EXTEND]"),
            Self::Output => write!(f, "[OUTPUT]"),
        }
    }
}

/// Reasons why a seed, candidate, or hit was filtered out
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FilterReason {
    // Seed generation (find_seeds_for_query)
    SeedContainsN,

    // Candidate validation (process_candidate)
    SeedOutOfBounds,

    // Maximality check (extend_seed)
    MaximalityLeft,
    MaximalityRight,

    // Energy filtering (process_candidate)
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

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SearchHit {
    pub query_id: QueryId,
    pub target_id: TargetId,

    pub q_start: usize,        // 0-based internal
    pub q_end: usize,          // 0-based internal
    pub t_start: usize,        // 0-based internal
    pub t_end: usize,          // 0-based internal
    pub output_q_start: usize, // 1-based for output/comparison
    pub output_q_end: usize,   // 1-based for output/comparison
    pub output_t_start: usize, // 1-based, strand-aware
    pub output_t_end: usize,   // 1-based, strand-aware
    pub strand: Strand,
    pub energy: Energy,
    pub alignment: Alignment,
    pub flank_5: String,
    pub flank_3: String,
}

impl SearchHit {
    pub fn write(&self, w: &mut dyn Write) -> std::io::Result<()> {
        // Write fixed fields
        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t",
            self.query_id.truncated(),
            self.output_q_start,
            self.output_q_end,
            self.target_id.truncated(),
            self.output_t_start,
            self.output_t_end,
            self.strand,
            self.energy.as_f64(),
        )?;

        // Write fingerprint directly (no String allocation)
        for p in self.alignment.steps() {
            write!(w, "{}", p.to_char())?;
        }
        write!(w, "\t")?;

        // Write target sequence with T->U normalization (no intermediate String)
        for p in self.alignment.steps() {
            let c = p.target_char();
            let normalized = match c {
                'T' => 'U',
                't' => 'u',
                other => other,
            };
            write!(w, "{}", normalized)?;
        }

        // Write remaining fields
        writeln!(w, "\t{}\t{}", self.flank_5, self.flank_3)
    }

    /// Parse a SearchHit from C risearch output line.
    ///
    /// C output format (tab-separated):
    /// `q_id, q_start, q_end, t_id, t_start, t_end, strand, energy, interaction, target_seq, [flank_5, flank_3]`
    ///
    /// Handles C quirks:
    /// - Seed markers 'y' and 'x' in interaction/target strings
    /// - Missing optional columns (flanks)
    /// - 1-based coordinates
    pub fn from_c_output(line: &str) -> Option<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            return None;
        }

        // Parse and strip seed markers from interaction
        let (interaction, seed_start, seed_end) = Self::strip_c_markers(fields[8]);
        let target_seq = Self::strip_markers_simple(fields[9]);

        // Parse coordinates (C uses 1-based)
        let q_start: usize = fields[1].parse().ok()?;
        let q_end: usize = fields[2].parse().ok()?;
        let t_start: usize = fields[4].parse().ok()?;
        let t_end: usize = fields[5].parse().ok()?;

        // Parse strand
        let strand = fields[6].chars().next().unwrap_or('+').into();

        // Parse energy
        let energy = Energy::parse(fields[7])?;

        // Create alignment from fingerprint and target sequence
        let alignment = Alignment::from_c_output(&interaction, &target_seq, seed_start, seed_end);

        // Optional flanks
        let flank_5 = fields
            .get(10)
            .map_or(String::new(), |s| Self::strip_markers_simple(s));
        let flank_3 = fields
            .get(11)
            .map_or(String::new(), |s| Self::strip_markers_simple(s));

        Some(SearchHit {
            query_id: fields[0].into(),
            target_id: fields[3].into(),
            q_start: q_start.saturating_sub(1), // 0-based internal
            q_end: q_end.saturating_sub(1),     // 0-based internal
            t_start: t_start.saturating_sub(1), // 0-based internal
            t_end: t_end.saturating_sub(1),     // 0-based internal
            output_q_start: q_start,            // 1-based for comparison
            output_q_end: q_end,                // 1-based for comparison
            output_t_start: t_start,            // 1-based
            output_t_end: t_end,                // 1-based
            strand,
            energy,
            alignment,
            flank_5,
            flank_3,
        })
    }

    /// Strip y/x seed markers and extract seed range positions.
    fn strip_c_markers(s: &str) -> (String, Option<usize>, Option<usize>) {
        let y_pos = s.find('y');
        let x_pos = s.find('x');

        let (seed_start, seed_end) = match (y_pos, x_pos) {
            (Some(y), Some(x)) if x > y => (Some(y), Some(x - 1)),
            _ => (None, None),
        };

        let clean: String = s.chars().filter(|&c| c != 'y' && c != 'x').collect();
        (clean, seed_start, seed_end)
    }

    /// Strip y/x markers without tracking positions.
    fn strip_markers_simple(s: &str) -> String {
        s.chars().filter(|&c| c != 'y' && c != 'x').collect()
    }

    // =========================================================================
    // COMPARISON HELPERS (for parity testing)
    // =========================================================================

    /// Returns true if coordinates and strand match another hit.
    /// Uses output coordinates (1-based) for comparison.
    pub fn coords_match(&self, other: &Self) -> bool {
        self.output_q_start == other.output_q_start
            && self.output_q_end == other.output_q_end
            && self.output_t_start == other.output_t_start
            && self.output_t_end == other.output_t_end
            && self.strand == other.strand
    }

    /// Group key for matching hits (query_id:target_id).
    pub fn group_key(&self) -> String {
        format!(
            "{}:{}",
            self.query_id.truncated(),
            self.target_id.truncated()
        )
    }

    /// Fingerprint string for comparison.
    pub fn fingerprint(&self) -> String {
        self.alignment.fingerprint()
    }

    /// Target sequence for comparison.
    pub fn target_seq(&self) -> String {
        self.alignment.target_sequence()
    }

    /// Query sequence for comparison (may have N placeholders if from C output).
    pub fn query_seq(&self) -> String {
        self.alignment.query_sequence()
    }

    /// Seed start position within interaction.
    pub fn seed_start(&self) -> Option<usize> {
        Some(self.alignment.left_extension().len())
    }

    /// Seed end position within interaction.
    pub fn seed_end(&self) -> Option<usize> {
        let start = self.alignment.left_extension().len();
        let seed_len = self.alignment.seed().len();
        Some(start + seed_len)
    }

    /// Format for debug output (uses 1-based output coordinates).
    pub fn fmt_coords(&self) -> String {
        format!(
            "q=[{},{}] t=[{},{}] S={} E={}",
            self.output_q_start,
            self.output_q_end,
            self.output_t_start,
            self.output_t_end,
            self.strand,
            self.energy
        )
    }
}

// Reimplementing mapping locally for safety and speed

pub struct SaIndex<'a> {
    pub index: &'a SaIndexFile,
}

impl<'a> SaIndex<'a> {
    pub fn find_candidates(&self, seed: &[u8], pairing: SeedPairing) -> Vec<SeedCandidate> {
        let mut candidates = Vec::new();
        // Normalize seed to lowercase for matching (since index uses lowercase)
        // Also normalize U -> T since the index stores DNA (T) not RNA (U)
        let seed_normalized: Vec<u8> = seed
            .iter()
            .map(|&b| {
                let lower = b.to_ascii_lowercase();
                if lower == b'u' { b't' } else { lower }
            })
            .collect();

        trace!(
            "{} seed={} pairing={:?}",
            SearchStage::Seed,
            String::from_utf8_lossy(&seed_normalized),
            pairing
        );

        for (i, seq_idx) in self.index.sequences.iter().enumerate() {
            let fwd_count_before = candidates.len();

            // 1. Search FORWARD strand
            self.search_sa_simple(
                &seq_idx.forward_sa,
                &seq_idx.sequence,
                &seed_normalized,
                pairing,
                i,
                Strand::Forward,
                &mut candidates,
            );

            let fwd_count = candidates.len() - fwd_count_before;

            // 2. Search REVERSE COMPLEMENT (use pre-computed from index)
            let rc_seq = &seq_idx.sequence_rc;
            let rc_count_before = candidates.len();

            self.search_sa_simple(
                &seq_idx.reverse_sa,
                &rc_seq,
                &seed_normalized,
                pairing,
                i,
                Strand::Reverse,
                &mut candidates,
            );

            let rc_count = candidates.len() - rc_count_before;

            trace!(
                "{} seq_idx={} name={} fwd_hits={} rc_hits={}",
                SearchStage::Seed,
                i,
                &seq_idx.name,
                fwd_count,
                rc_count
            );
        }

        trace!(
            "{} total_candidates={}",
            SearchStage::Seed,
            candidates.len()
        );
        candidates
    }

    pub fn get_sequence(&self, seq_idx: usize) -> &[u8] {
        &self.index.sequences[seq_idx].sequence
    }

    pub fn get_sequence_rc(&self, seq_idx: usize) -> &[u8] {
        &self.index.sequences[seq_idx].sequence_rc
    }

    pub fn get_id(&self, seq_idx: usize) -> &str {
        &self.index.sequences[seq_idx].name
    }

    pub fn get_sequence_len(&self, seq_idx: usize) -> usize {
        self.index.sequences[seq_idx].sequence.len()
    }

    /// Simple SA search - returns positions where seed matches
    #[allow(clippy::too_many_arguments)]
    fn search_sa_simple(
        &self,
        sa: &[u32],
        text: &[u8],
        seed: &[u8],
        pairing: SeedPairing,
        seq_idx: usize,
        strand: Strand,
        candidates: &mut Vec<SeedCandidate>,
    ) {
        trace!(
            "{} START seed={} pairing={:?} idx={} strand={:?}",
            SearchStage::Seed,
            String::from_utf8_lossy(seed),
            pairing,
            seq_idx,
            strand
        );
        let mut stack = vec![(0, sa.len(), 0)];

        while let Some((start, end, offset)) = stack.pop() {
            if offset == seed.len() {
                // Match found - add all positions in this SA range
                for &sa_pos in &sa[start..end] {
                    candidates.push(SeedCandidate {
                        query_pos: 0, // Set by caller
                        target_idx: seq_idx,
                        target_start: sa_pos as usize,
                        len: seed.len(),
                        strand,
                    });
                }
                continue;
            }

            let target_char = seed[offset];

            // Canonical Match
            let (s, e) = self.get_sa_interval(sa, text, start, end, offset, target_char);
            if s < e {
                stack.push((s, e, offset + 1));
            }

            // Wobble Match (if Allowed)
            if matches!(pairing, SeedPairing::AllowWobble) {
                let wobble_char = match target_char {
                    b'c' => Some(b't'),
                    b'a' => Some(b'g'),
                    _ => None,
                };

                trace!(
                    "{} Wobble Check offset={} char={} pairing={:?} wobble_char={:?}",
                    SearchStage::Seed,
                    offset,
                    target_char as char,
                    pairing,
                    wobble_char
                );

                if let Some(wc) = wobble_char {
                    let (ws, we) = self.get_sa_interval(sa, text, start, end, offset, wc);
                    if ws < we {
                        trace!(
                            "{} Wobble Found! offset={} range={}-{}",
                            SearchStage::Seed,
                            offset,
                            ws,
                            we
                        );
                        stack.push((ws, we, offset + 1));
                    }
                }
            }
        }
    }

    /// Helper to find range in SA for a specific character at offset
    fn get_sa_interval(
        &self,
        sa: &[u32],
        text: &[u8],
        start: usize,
        end: usize,
        offset: usize,
        target: u8,
    ) -> (usize, usize) {
        // Find lower bound
        // We want the first index where text[sa[i] + offset] >= target
        let lower = sa[start..end].partition_point(|&idx| {
            let pos = (idx as usize) + offset;
            let c = if pos < text.len() { text[pos] } else { 0 };
            c < target
        });

        // Find upper bound
        // We want the first index where text[sa[i] + offset] > target
        let upper = sa[start..end].partition_point(|&idx| {
            let pos = (idx as usize) + offset;
            let c = if pos < text.len() { text[pos] } else { 0 };
            c <= target
        });

        (start + lower, start + upper)
    }
}

pub struct SearchContext<'a, 'e> {
    pub index: &'a SaIndex<'a>,
    pub args: &'a SearchArgs,
    pub extender: &'e mut dp::DpExtender,
    pub stats: SearchStats,
    pub energy: EnergyModel,
}

impl<'a, 'e> SearchContext<'a, 'e> {
    pub fn with_extender(
        index: &'a SaIndex<'a>,
        args: &'a SearchArgs,
        extender: &'e mut dp::DpExtender,
    ) -> Self {
        Self {
            index,
            args,
            extender,
            stats: SearchStats::default(),
            energy: EnergyModel::T04,
        }
    }
}

/// Thread-local storage for reusable buffers - avoids allocation per query.
thread_local! {
    static THREAD_EXTENDER: std::cell::RefCell<dp::DpExtender> =
        std::cell::RefCell::new(dp::DpExtender::new());
    static THREAD_SEEDS: std::cell::RefCell<Vec<SeedCandidate>> =
        std::cell::RefCell::new(Vec::with_capacity(100_000));
    static THREAD_MATCHES: std::cell::RefCell<Vec<ParallelSeedMatch>> =
        std::cell::RefCell::new(Vec::with_capacity(1024));
}

/// Core search logic - finds and deduplicates hits for all queries.
/// Uses parallel iteration over queries for multi-core utilization.
fn search_core(
    queries: &[(String, Vec<u8>)],
    index: &SaIndex<'_>,
    opts: &SearchArgs,
) -> Result<Vec<SearchHit>> {
    info!(
        "Starting search: {} queries, seed={:?}, max_ext={}, delta_g={}",
        queries.len(),
        opts.seed.seed,
        opts.extend.max_extension,
        opts.extend.delta_g
    );

    // Sequential search for fair single-threaded comparison.
    // Use par_iter() and RAYON_NUM_THREADS for multi-threaded mode.
    let all_hits: Vec<SearchHit> = queries
        .iter() // Sequential - use par_iter() for parallel
        .flat_map(|(q_id, q_seq)| {
            // Borrow thread-local extender for this query
            THREAD_EXTENDER.with(|ext| {
                let mut extender = ext.borrow_mut();
                let mut ctx = SearchContext::with_extender(index, opts, &mut *extender);

                trace!("{} id={} len={}", SearchStage::Input, q_id, q_seq.len());
                trace!("{} {}", SearchStage::Input, String::from_utf8_lossy(q_seq));

                let seeds = match find_seeds_for_query(q_seq, &mut ctx) {
                    Ok(s) => s,
                    Err(_) => return Vec::new(),
                };
                trace!("{} {} candidates found", SearchStage::Input, seeds.len());

                let mut hits = Vec::new();
                for candidate in &seeds {
                    trace!(
                        "{} q_pos={} t_idx={} t_start={} len={} strand={:?}",
                        SearchStage::Seed,
                        candidate.query_pos,
                        candidate.target_idx,
                        candidate.target_start,
                        candidate.len,
                        candidate.strand
                    );
                    if let Some(hit) = process_candidate(
                        crate::types::Query::new(q_id, q_seq),
                        candidate,
                        &mut ctx,
                    ) {
                        trace!(
                            "{} q={}-{} t={}-{} E={:.2}",
                            SearchStage::Output,
                            hit.q_start,
                            hit.q_end,
                            hit.t_start,
                            hit.t_end,
                            hit.energy
                        );
                        hits.push(hit);
                    }
                }
                hits
            })
        })
        .collect();

    info!("Search complete: {} hits", all_hits.len());

    Ok(all_hits)
}

/// Run search and return hits.
pub fn run_search(
    queries: &[(String, Vec<u8>)],
    index: &SaIndex<'_>,
    opts: &SearchArgs,
) -> Result<Vec<SearchHit>> {
    search_core(queries, index, opts)
}

/// Run search with streaming output - writes hits directly instead of collecting.
/// This avoids memory overhead for large result sets.
pub fn run_search_streaming<W: std::io::Write>(
    queries: &[(String, Vec<u8>)],
    index: &SaIndex<'_>,
    opts: &SearchArgs,
    writer: &mut W,
) -> Result<usize> {
    info!(
        "Starting streaming search: {} queries, seed={:?}, max_ext={}, delta_g={}",
        queries.len(),
        opts.seed.seed,
        opts.extend.max_extension,
        opts.extend.delta_g
    );

    let mut hit_count = 0;
    let mut line_buf: Vec<u8> = Vec::with_capacity(512);

    for (q_id, q_seq) in queries {
        // Borrow thread-local buffers for this query
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
                        // Use streaming variant - zero String allocations
                        if process_candidate_streaming(
                            crate::types::Query::new(q_id, q_seq),
                            candidate,
                            &mut ctx,
                            writer,
                            &mut line_buf,
                        ) {
                            hit_count += 1;
                        }
                    }
                });
            });
        });
    }

    info!("Streaming search complete: {} hits", hit_count);
    Ok(hit_count)
}

pub fn write_results(hits: &[SearchHit], output: impl AsRef<Path>) -> Result<()> {
    use std::io::BufWriter;

    let inner: Box<dyn Write> = if output.as_ref() == Path::new("-") {
        Box::new(std::io::stdout())
    } else {
        Box::new(std::fs::File::create(output.as_ref()).context("Failed to create output file")?)
    };

    // Larger buffer reduces syscall overhead for high-volume output.
    let mut writer = BufWriter::with_capacity(256 * 1024, inner);

    debug!("{} output={:?}", SearchStage::Output, output.as_ref());

    for hit in hits {
        hit.write(&mut writer)?;
    }

    // BufWriter flushes on drop, but explicit flush ensures errors are caught
    writer.flush().context("Failed to flush output")?;
    Ok(())
}

/// Check if hit `k` shadows hit `h`.
///
/// # Deduplication behavior
///
/// **Energy-different duplicates at same coordinates:**
/// C does NOT filter hits with identical coordinates but different energies.
/// This can occur when different extension paths (varying seeds or traceback
/// choices) produce alignments ending at the same coordinates with different
/// total energies. Both are valid biological findings that may represent
/// different sub-optimal solutions.
///
/// Whether this is "better" is use-case dependent:
/// - **Preserving all**: More complete output, shows alternative alignments
/// - **Dedup to best**: Cleaner output, but loses sub-optimal alternatives
///
/// We match C behavior (preserve all) for parity. Only true duplicates
/// (same coords AND same energy within 0.001) are filtered.
///

fn find_seeds_for_query(
    q_seq: &[u8],
    ctx: &mut SearchContext<'_, '_>,
) -> Result<Vec<SeedCandidate>> {
    use crate::seed;

    trace!(
        "{} q_len={} spec={:?} pairing={:?}",
        SearchStage::Seed,
        q_seq.len(),
        ctx.args.seed.seed,
        ctx.args.seed.pairing
    );

    Ok(seed::find_seeds(q_seq, ctx.index.index, &ctx.args.seed))
}

fn process_candidate(
    query: crate::types::Query<'_>,
    candidate: &SeedCandidate,
    ctx: &mut SearchContext<'_, '_>,
) -> Option<SearchHit> {
    let t_idx = candidate.target_idx;
    let t_seq: &[u8] = match candidate.strand {
        Strand::Reverse => ctx.index.get_sequence_rc(t_idx),
        Strand::Forward => ctx.index.get_sequence(t_idx),
    };

    let t_start_idx = candidate.target_start;
    let seed_len = candidate.len;
    let q_pos = candidate.query_pos;

    trace!(
        "{} q_id={} t_idx={} q_pos={} t_start={} seed_len={} strand={:?}",
        SearchStage::Extend,
        query.id,
        t_idx,
        q_pos,
        t_start_idx,
        seed_len,
        candidate.strand
    );

    if t_start_idx + seed_len > t_seq.len() {
        ctx.stats.record_filter(FilterReason::SeedOutOfBounds);
        warn!(
            "{} FILTERED reason={:?} t_start={} seed_len={} t_len={}",
            SearchStage::Extend,
            FilterReason::SeedOutOfBounds,
            t_start_idx,
            seed_len,
            t_seq.len()
        );
        return None;
    }

    // Call extend_seed (or essentially reproduce its valuable logic).
    let extension_result = extend_seed(ctx, query.seq, t_seq, candidate)?;

    let ext = extension_result;
    let score = ext.score;

    if score > ctx.args.extend.delta_g {
        ctx.stats.record_filter(FilterReason::EnergyAboveThreshold);
        trace!(
            "{} FILTERED reason={:?} score={:.2} > delta_g={}",
            SearchStage::Extend,
            FilterReason::EnergyAboveThreshold,
            score,
            ctx.args.extend.delta_g
        );
        return None;
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
            // Antisense hit / Reverse Strand
            let fwd_start = original_len - 1 - final_t_end;
            let fwd_end = original_len - 1 - final_t_start;
            (fwd_start + 1, fwd_end + 1, '-')
        }
        Strand::Forward => {
            // Forward hit
            (final_t_start + 1, final_t_end + 1, '+')
        }
    };

    // Flanks (20bp context, T->U normalized)
    let ctx_len = 20;
    let flank_5 = bytes_to_rna_string(
        &t_seq[final_t_start.saturating_sub(ctx_len)..final_t_start],
        true, // reverse
    );
    let t_3_end = (final_t_end + 1 + ctx_len).min(t_seq.len());
    let flank_3 = if final_t_end + 1 < t_seq.len() {
        bytes_to_rna_string(&t_seq[final_t_end + 1..t_3_end], false)
    } else {
        String::new()
    };

    Some(SearchHit {
        query_id: query.id.into(),
        target_id: ctx.index.get_id(t_idx).into(),

        q_start: final_q_start,
        q_end: final_q_end,
        t_start: final_t_start,
        t_end: final_t_end,
        output_q_start: final_q_start + 1, // 1-based for output
        output_q_end: final_q_end + 1,     // 1-based for output
        output_t_start: out_t_start,
        output_t_end: out_t_end,
        strand: strand_char.into(),
        energy: score.into(),
        alignment: ext.alignment,
        flank_5,
        flank_3,
    })
}

/// Process candidate and write directly to output - zero String allocation variant.
/// Returns true if a hit was written.
fn process_candidate_streaming<W: Write>(
    query: crate::types::Query<'_>,
    candidate: &SeedCandidate,
    ctx: &mut SearchContext<'_, '_>,
    writer: &mut W,
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
        Strand::Forward => {
            (final_t_start + 1, final_t_end + 1, '+')
        }
    };

    // Write directly - no String allocations
    let target_id = ctx.index.get_id(t_idx);
    let ctx_len = 20;

    // Write fixed fields (using truncated IDs)
    let q_id_trunc = if query.id.len() > 50 { &query.id[..50] } else { query.id };
    let t_id_trunc = if target_id.len() > 50 { &target_id[..50] } else { target_id };

    line_buf.clear();

    if write!(
        line_buf,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t",
        q_id_trunc,
        final_q_start + 1,
        final_q_end + 1,
        t_id_trunc,
        out_t_start,
        out_t_end,
        strand_char,
        score,
    )
    .is_err()
    {
        return false;
    }

    // Write fingerprint directly
    for p in ext.alignment.steps() {
        line_buf.push(p.to_char() as u8);
    }
    line_buf.push(b'\t');

    // Write target sequence with T->U normalization
    for p in ext.alignment.steps() {
        let c = p.target_char();
        let normalized = match c {
            'T' => 'U',
            't' => 'u',
            other => other,
        };
        line_buf.push(normalized as u8);
    }
    line_buf.push(b'\t');

    // Write flank_5 (reversed, T->U)
    let flank_5_slice = &t_seq[final_t_start.saturating_sub(ctx_len)..final_t_start];
    if write_bytes_as_rna(line_buf, flank_5_slice, true).is_err() {
        return false;
    }
    line_buf.push(b'\t');

    // Write flank_3 (T->U)
    let t_3_end = (final_t_end + 1 + ctx_len).min(t_seq.len());
    if final_t_end + 1 < t_seq.len() {
        if write_bytes_as_rna(line_buf, &t_seq[final_t_end + 1..t_3_end], false).is_err() {
            return false;
        }
    }
    line_buf.push(b'\n');

    if writer.write_all(line_buf).is_err() {
        return false;
    }

    true
}

#[derive(Debug, Clone)]
pub struct ExtensionResult {
    pub score: f64,
    pub alignment: Alignment,

    pub l_q: usize,
    pub l_t: usize,
    pub r_q: usize,
    pub r_t: usize,
}

fn extend_seed(
    ctx: &mut SearchContext<'_, '_>,
    q_seq: &[u8],
    t_seq: &[u8],
    candidate: &SeedCandidate,
) -> Option<ExtensionResult> {
    let q_pos = candidate.query_pos;
    let t_pos = candidate.target_start;
    let len = candidate.len;
    let opts = &ctx.args.extend;

    // Get interval bounds from SeedSpec for maximality constraint
    // For -s 8: interval is [0, q_len), maximality check uses full query
    // For -s 1:8: interval is [0, 8), seeds at edges are maximal within interval
    let (interval_start, interval_end) = ctx
        .args
        .seed
        .seed
        .normalize(q_seq.len())
        .map(|(s1, e1, _)| (s1 - 1, e1)) // Convert to 0-based start, exclusive end
        .unwrap_or((0, q_seq.len())); // Fallback: full query

    // Wrap sequences for clean base access (used throughout function)
    let query = Seq::forward(q_seq);
    let target = Seq::new(t_seq, candidate.strand);

    // Select pair matrix based on wobble mode:
    // - AllowWobble: use PAIR_MAT (G-U wobble pairs are valid)
    // - Strict: use PAIR_MAT_NO_GU (G-U wobble pairs are NOT valid)
    let pair_mat = match ctx.args.seed.pairing {
        SeedPairing::AllowWobble => &PAIR_MAT,
        SeedPairing::Strict => &PAIR_MAT_NO_GU,
    };

    // MAXIMALITY CHECK
    // Skip non-maximal seeds: if the seed can be extended by a valid base pair
    // on either end, it's a sub-seed of a longer match and will

    trace!(
        "{} ENTERING extend_seed q_pos={} t_pos={} len={} delta_g={}",
        SearchStage::Extend,
        q_pos,
        t_pos,
        len,
        opts.delta_g
    );

    // 1. Left extendable? Check if q[q_pos-1] pairs with t[t_pos+len]
    // Use interval_start to constrain: seed touching interval left edge is maximal on left
    if q_pos > interval_start && t_pos + len < t_seq.len() {
        let q_prev = query.base(q_pos - 1).idx();
        let t_next = target.base(t_pos + len).idx();
        let p_class = pair_mat[q_prev][t_next];

        trace!(
            "{} Left Check q_pos={} t_pos={} len={} q_prev={} t_next={} pair={}",
            SearchStage::Extend,
            q_pos,
            t_pos,
            len,
            q_prev,
            t_next,
            p_class
        );

        if p_class != 0 {
            if opts.no_max_prune {
                trace!(
                    "SEED: Non-maximal Left-Ext: q_pos={} t_pos={} len={} pair={} (not pruning due to --no-max-prune)",
                    q_pos, t_pos, len, p_class
                );
            } else {
                ctx.stats.record_filter(FilterReason::MaximalityLeft);
                trace!(
                    "{} FILTERED reason={:?} q_pos={} t_pos={} len={} pair={}",
                    SearchStage::Extend,
                    FilterReason::MaximalityLeft,
                    q_pos,
                    t_pos,
                    len,
                    p_class
                );
                return None;
            }
        }
    }

    // 2. Right extendable? Check if q[q_pos+len] pairs with t[t_pos-1]
    // Use interval_end to constrain: seed touching interval right edge is maximal on right
    if q_pos + len < interval_end && t_pos > 0 {
        let q_next = query.base(q_pos + len).idx();
        let t_prev = target.base(t_pos - 1).idx();
        let p_class = pair_mat[q_next][t_prev];

        trace!(
            "{} Right Check q_pos={} t_pos={} len={} q_next={} t_prev={} pair={}",
            SearchStage::Extend,
            q_pos,
            t_pos,
            len,
            q_next,
            t_prev,
            p_class
        );

        if p_class != 0 {
            if opts.no_max_prune {
                trace!(
                    "SEED: Non-maximal Right-Ext: q_pos={} t_pos={} len={} pair={} (not pruning due to --no-max-prune)",
                    q_pos, t_pos, len, p_class
                );
            } else {
                ctx.stats.record_filter(FilterReason::MaximalityRight);
                trace!(
                    "{} FILTERED reason={:?} q_pos={} t_pos={} len={} pair={}",
                    SearchStage::Extend,
                    FilterReason::MaximalityRight,
                    q_pos,
                    t_pos,
                    len,
                    p_class
                );
                return None;
            }
        }
    }

    trace!(
        "{} Accepted Seed: q_pos={} t_pos={} len={} delta_g={}",
        SearchStage::Extend,
        q_pos,
        t_pos,
        len,
        opts.delta_g
    );

    let t_match_end = t_pos + len - 1;
    let max_ext = opts.max_extension as usize;
    let safe_ext = max_ext.min(MAX_DP_EXT);

    // Check extendability: seed is extendable if there's room on at least one side
    let left_extendable = q_pos > 0 && t_pos + len < t_seq.len();
    let right_extendable = q_pos + len < q_seq.len() && t_pos > 0;
    let seed_extendable = left_extendable || right_extendable;

    // Only enter extension if max_extension is nonzero AND seed is extendable
    let do_extension = max_ext > 0 && seed_extendable;

    if !do_extension {
        trace!(
            "{} SKIPPING EXTENSION: max_ext={} extendable={} (L={} R={})",
            SearchStage::Extend,
            max_ext,
            seed_extendable,
            left_extendable,
            right_extendable
        );
    }

    // Seed energy calculation using centralized EnergyModel
    let seed_energy_raw = ctx
        .energy
        .seed_energy(&query, &target, q_pos, t_match_end, len);

    // Only build interaction string when trace logging is enabled (avoids allocation in hot path)
    if log::log_enabled!(log::Level::Trace) {
        let mut seed_int_str = String::with_capacity(len);
        for k in 0..len {
            let qc = query.base(q_pos + k);
            let tc = target.base(t_match_end - k);
            seed_int_str.push(qc.pairing_class(tc));
        }
        trace!(
            "[SEED] energy_raw={} q_pos={} t_end={} len={} interaction={}",
            seed_energy_raw, q_pos, t_match_end, len, seed_int_str
        );
    }

    // If not entering extension, return seed-only result
    // Need to add terminal penalties at both ends of the seed
    if !do_extension {
        // 5' terminal: Gap->q[0] paired with Gap->t[t_end]
        let term_5p = ctx
            .energy
            .terminal_5p(query.base(q_pos), target.base(t_match_end));
        // 3' terminal: q[q_end]->Gap paired with t[t_pos]->Gap
        let term_3p = ctx
            .energy
            .terminal_3p(query.base(q_pos + len - 1), target.base(t_pos));

        let total_raw = seed_energy_raw + term_5p + term_3p;
        let final_score = ctx.energy.to_kcal(total_raw);
        let seed_alignment = build_seed_alignment(&query, &target, q_pos, t_match_end, len);

        return Some(ExtensionResult {
            score: final_score,
            alignment: Alignment::new(&[], &seed_alignment, &[]),
            l_q: 0,
            l_t: 0,
            r_q: 0,
            r_t: 0,
        });
    }

    // DP Left: Extend Query Left (5'), Target Right (3')
    let left_res = ctx
        .extender
        .extend_left(&query, &target, q_pos, t_match_end, safe_ext);

    // DP Right: Extend Query Right (3'), Target Left (5')
    let right_res = ctx
        .extender
        .extend_right(&query, &target, q_pos + len - 1, t_pos, safe_ext);

    // Use OLD embedded DP scores for C parity (new dp module has differences)
    let final_score = ctx
        .energy
        .to_kcal(seed_energy_raw + left_res.score + right_res.score);

    trace!(
        "{} score={:.2} (seed={:.2} L={} R={}) L_len={}/{} R_len={}/{}",
        SearchStage::Extend,
        final_score,
        seed_energy_raw as f64 / -100.0,
        left_res.score,
        right_res.score,
        left_res.q_len,
        left_res.t_len,
        right_res.q_len,
        right_res.t_len
    );

    // Use SmallVec to avoid heap allocation for typical extension sizes
    use smallvec::SmallVec;
    let mut left_alignment: SmallVec<[Pairing; 64]> = SmallVec::new();

    // Left Trace: replay from extension end back toward seed
    let mut li = left_res.q_len;
    let mut lj = left_res.t_len;
    for step in &left_res.trace {
        let pairing = match step {
            dp::DpOp::Match | dp::DpOp::Stop => {
                let p = Pairing::from_bases(
                    query.left(candidate.query_pos, li),
                    target.base_or_gap(t_match_end + lj),
                );
                li = li.saturating_sub(1);
                lj = lj.saturating_sub(1);
                p
            }
            dp::DpOp::GapQ => {
                let p = Pairing::GapTarget(query.left(candidate.query_pos, li));
                li = li.saturating_sub(1);
                p
            }
            dp::DpOp::GapT => {
                let p = Pairing::GapQuery(target.base_or_gap(t_match_end + lj));
                lj = lj.saturating_sub(1);
                p
            }
        };
        left_alignment.push(pairing);
    }

    // 2. Seed itself
    let seed_alignment = build_seed_alignment(&query, &target, q_pos, t_match_end, len);

    // Right Trace: replay reversed trace (from seed toward extension end)
    let mut right_alignment: SmallVec<[Pairing; 64]> = SmallVec::new();
    let (mut ri, mut rj) = (0, 0);
    let q_base = candidate.query_pos + len - 1;
    for step in right_res.trace.iter().rev() {
        let pairing = match step {
            dp::DpOp::Match | dp::DpOp::Stop => {
                ri += 1;
                rj += 1;
                Pairing::from_bases(query.base_or_gap(q_base + ri), target.left(t_pos, rj))
            }
            dp::DpOp::GapQ => {
                ri += 1;
                Pairing::GapTarget(query.base_or_gap(q_base + ri))
            }
            dp::DpOp::GapT => {
                rj += 1;
                Pairing::GapQuery(target.left(t_pos, rj))
            }
        };
        right_alignment.push(pairing);
    }

    let alignment = Alignment::new(&left_alignment, &seed_alignment, &right_alignment);

    Some(ExtensionResult {
        score: final_score,
        alignment,

        l_q: left_res.q_len,
        l_t: left_res.t_len,
        r_q: right_res.q_len,
        r_t: right_res.t_len,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_search_hit_from_c_output() {
        // Sample C output line (C uses 1-based coords)
        let line =
            "hsa-miR-1\t1\t10\tENSG00000001\t100\t110\t+\t-15.50\tyPPPUPPPx\tacguacgu\tAA\tCC";

        let hit = SearchHit::from_c_output(line).expect("should parse");

        assert_eq!(hit.query_id.as_str(), "hsa-miR-1");
        assert_eq!(hit.target_id.as_str(), "ENSG00000001");
        // Internal coords are 0-based (converted from C's 1-based)
        assert_eq!(hit.q_start, 0);
        assert_eq!(hit.q_end, 9);
        assert_eq!(hit.t_start, 99);
        assert_eq!(hit.t_end, 109);
        // Output coords stay 1-based for display
        assert_eq!(hit.output_q_start, 1);
        assert_eq!(hit.output_q_end, 10);
        assert_eq!(hit.output_t_start, 100);
        assert_eq!(hit.output_t_end, 110);
        assert_eq!(hit.strand, Strand::Forward);
        assert!((hit.energy.as_f64() - (-15.50)).abs() < 0.01);
        assert_eq!(hit.alignment.fingerprint(), "PPPUPPP"); // markers stripped
        assert_eq!(hit.flank_5, "AA");
        assert_eq!(hit.flank_3, "CC");
    }

    #[test]
    fn test_search_hit_from_c_output_minimal() {
        // Minimal 10 columns (no flanks), C uses 1-based coords
        let line = "q1\t1\t5\tt1\t10\t15\t-\t-8.00\tPPPPP\tacgua";

        let hit = SearchHit::from_c_output(line).expect("should parse minimal");
        assert_eq!(hit.query_id.as_str(), "q1");
        // Internal coords are 0-based
        assert_eq!(hit.q_start, 0);
        assert_eq!(hit.q_end, 4);
        assert_eq!(hit.t_start, 9);
        assert_eq!(hit.t_end, 14);
        // Output coords stay 1-based
        assert_eq!(hit.output_q_start, 1);
        assert_eq!(hit.output_q_end, 5);
        assert_eq!(hit.output_t_start, 10);
        assert_eq!(hit.output_t_end, 15);
        assert_eq!(hit.strand, Strand::Reverse);
        assert_eq!(hit.flank_5, "");
        assert_eq!(hit.flank_3, "");
    }

    #[test]
    fn test_search_hit_from_c_output_invalid() {
        assert!(SearchHit::from_c_output("too\tfew\tcolumns").is_none());
        assert!(SearchHit::from_c_output("").is_none());
    }
}
