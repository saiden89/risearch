use anyhow::{Context, Result, anyhow};
use log::{debug, trace};
use std::io::Write;
use std::path::Path;
use std::str::FromStr;

use crate::dsm::{Base, DSM_T04_POS, PAIR_MAT};
use crate::sa::IndexFile;
use crate::seed::SeedSpec;
use crate::{ExtendArgs, SearchArgs, SeedArgs};

const MAX_DP_EXT: usize = 30;
const GAP_IDX: usize = Base::Gap as usize;

pub trait Sequence {
    fn reverse_complement_dna(&self) -> Vec<u8>;
    fn reverse_complement_rna(&self) -> Vec<u8>;
}

impl Sequence for [u8] {
    fn reverse_complement_dna(&self) -> Vec<u8> {
        self.iter()
            .rev()
            .map(|&b| {
                // Return lowercase T for A (mimicking old complement_dna behavior)
                // Base::from_byte(b).complement().to_u8_upper().to_ascii_lowercase()
                // But wait, Base::U -> b'U'. to_ascii_lowercase -> b'u'.
                // Index stores DNA T. So we need to map U -> T if it's U.
                let c = Base::from_byte(b).complement();
                match c {
                    Base::U => b't',
                    _ => c.to_u8_upper().to_ascii_lowercase(),
                }
            })
            .collect()
    }

    fn reverse_complement_rna(&self) -> Vec<u8> {
        self.iter()
            .rev()
            .map(|&b| Base::from_byte(b).complement().to_u8_upper())
            .collect()
    }
}

pub struct SeedCandidate {
    pub query_pos: usize,
    pub target_idx: usize,
    pub target_start: usize, // 0-based index in target
    pub len: usize,
    pub is_antisense: bool,
}

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SearchHit {
    pub query_id: String,
    pub target_id: String,
    pub query_seq: String,
    pub target_seq: String,
    pub q_start: usize,
    pub q_end: usize,
    pub t_start: usize,
    pub t_end: usize,
    pub output_t_start: usize, // 1-based, strand-aware
    pub output_t_end: usize,   // 1-based, strand-aware
    pub strand: char,
    pub energy: f64,
    pub interaction: String, // The ASCII fingerprint
    pub flank_5: String,
    pub flank_3: String,
}

// Reimplementing mapping locally for safety and speed

pub struct SaIndex<'a> {
    pub index: &'a IndexFile,
}

impl<'a> SaIndex<'a> {
    pub fn find_candidates(&self, seed: &[u8], wobble: bool) -> Vec<(usize, usize, bool)> {
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

        for (i, seq_idx) in self.index.sequences.iter().enumerate() {
            // 1. Search FORWARD strand: query binds to sense strand of target
            //    Reported as '-' strand in C convention (antisense of transcript)
            self.search_sa_simple(
                &seq_idx.forward_sa,
                &seq_idx.sequence,
                &seed_normalized,
                wobble,
                i,
                false,
                &mut candidates,
            );

            // 2. Search REVERSE COMPLEMENT: query binds to antisense strand of target
            //    Reported as '+' strand in C convention (sense transcript)
            //    Need to build RC sequence to search against
            let rc_seq = seq_idx.sequence.reverse_complement_dna();
            self.search_sa_simple(
                &seq_idx.reverse_sa,
                &rc_seq,
                &seed_normalized,
                wobble,
                i,
                true,
                &mut candidates,
            );
        }
        candidates
    }

    pub fn get_sequence(&self, seq_idx: usize) -> &[u8] {
        &self.index.sequences[seq_idx].sequence
    }

    pub fn get_sequence_rc(&self, seq_idx: usize) -> Vec<u8> {
        self.index.sequences[seq_idx]
            .sequence
            .reverse_complement_dna()
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
        sa: &[i64],
        text: &[u8],
        seed: &[u8],
        wobble: bool,
        seq_idx: usize,
        is_antisense: bool,
        candidates: &mut Vec<(usize, usize, bool)>,
    ) {
        let mut stack = vec![(0, sa.len(), 0)];

        while let Some((start, end, offset)) = stack.pop() {
            if offset == seed.len() {
                // Match found - add all positions in this SA range
                for &sa_pos in &sa[start..end] {
                    candidates.push((seq_idx, sa_pos as usize, is_antisense));
                }
                continue;
            }

            let target_char = seed[offset];

            // Canonical Match
            let (s, e) = get_sa_interval(sa, text, start, end, offset, target_char);
            if s < e {
                stack.push((s, e, offset + 1));
            }

            // Wobble Match (when wobble=false, we ALLOW G-U pairing)
            if !wobble {
                let wobble_char = match target_char {
                    b'c' => Some(b't'), // Query G binds to target U (stored as T)
                    b'a' => Some(b'g'), // Query U binds to target G
                    _ => None,
                };

                if let Some(wc) = wobble_char {
                    let (ws, we) = get_sa_interval(sa, text, start, end, offset, wc);
                    if ws < we {
                        stack.push((ws, we, offset + 1));
                    }
                }
            }
        }
    }
}

/// DNA reverse complement (for target sequences stored as DNA with T not U)

/// Helper to find range in SA for a specific character at offset
fn get_sa_interval(
    sa: &[i64],
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

// --- DP Structure ---

/// DP alignment state - replaces magic integers (0=M, 1=Bq, 2=Bt)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum DpState {
    Match, // M - Match/Mismatch state (paired bases)
    GapQ,  // Bq - Query Bulge state (gap in target, query base unpaired)
    GapT,  // Bt - Target Bulge state (gap in query, target base unpaired)
}

#[derive(Clone, Debug)]
struct Grid<T> {
    data: Vec<T>,
    width: usize,
    height: usize,
}

impl<T: Clone + Copy + Default> Grid<T> {
    fn new(width: usize, height: usize) -> Self {
        Self {
            data: vec![T::default(); width * height],
            width,
            height,
        }
    }

    #[inline(always)]
    fn get(&self, r: usize, c: usize) -> T {
        self.data[r * self.width + c]
    }

    #[inline(always)]
    fn set(&mut self, r: usize, c: usize, val: T) {
        self.data[r * self.width + c] = val;
    }

    fn resize(&mut self, width: usize, height: usize) {
        let new_len = width * height;
        self.data.clear();
        self.data.resize(new_len, T::default());
        self.width = width;
        self.height = height;
    }
}

/// Type alias for score grids - uses Option<i32> instead of sentinel value
type ScoreGrid = Grid<Option<i32>>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TraceStep {
    Stop,
    Match, // From M(i-1, j-1)
    GapQ,  // From Bq(i, j)
    GapT,  // From Bt(i, j)
}

impl Default for TraceStep {
    fn default() -> Self {
        Self::Stop
    }
}

struct DpContext {
    // Score matrices (None = not reachable, Some(score) = reachable with score)
    m: ScoreGrid,
    bq: ScoreGrid,
    bt: ScoreGrid,

    // Traceback matrices
    tb_m: Grid<TraceStep>,
    tb_bq: Grid<TraceStep>,
    tb_bt: Grid<TraceStep>,
}

impl DpContext {
    fn new(width: usize, height: usize) -> Self {
        Self {
            m: ScoreGrid::new(width, height),
            bq: ScoreGrid::new(width, height),
            bt: ScoreGrid::new(width, height),
            tb_m: Grid::new(width, height),
            tb_bq: Grid::new(width, height),
            tb_bt: Grid::new(width, height),
        }
    }

    fn resize(&mut self, width: usize, height: usize) {
        self.m.resize(width, height);
        self.bq.resize(width, height);
        self.bt.resize(width, height);
        self.tb_m.resize(width, height);
        self.tb_bq.resize(width, height);
        self.tb_bt.resize(width, height);
    }
}

pub fn run_search(
    queries: &[(String, Vec<u8>)],
    index: &SaIndex<'_>,
    output: impl AsRef<Path>,
    opts: &SearchArgs,
) -> Result<()> {
    // Create output writer
    let mut writer: Box<dyn Write> = if output.as_ref() == Path::new("-") {
        Box::new(std::io::stdout())
    } else {
        Box::new(std::fs::File::create(output.as_ref()).context("Failed to create output file")?)
    };

    debug!("run_search called with {} queries", queries.len());

    // Create DP Context (reused buffers)
    // Max extension length is roughly max_ext * 2 + seed_len?
    // Actually we resize dynamically, but initial size helps.
    // Let's assume max extension around 100-200.
    let mut ctx = DpContext::new(200, 200);

    for (q_id, q_seq) in queries {
        debug!("Processing query: {} (len={})", q_id, q_seq.len());

        // Find seeds using the new helper
        let seeds = find_seeds_for_query(q_seq, index, &opts.seed)?;
        debug!("  Found {} seed candidates", seeds.len());

        let mut hit_count = 0;
        for candidate in seeds {
            if let Some(hit) =
                process_candidate(q_id, q_seq, index, &candidate, &opts.extend, &mut ctx)
            {
                print_search_hit(&mut writer, &hit)?;
                hit_count += 1;
            }
        }
        if hit_count > 0 {
            debug!("  Reported {} hits", hit_count);
        }
    }

    Ok(())
}

//TODO: improve performance via better search strategies!

fn find_seeds_for_query(
    q_seq: &[u8],
    index: &SaIndex<'_>,
    seed_args: &SeedArgs,
) -> Result<Vec<SeedCandidate>> {
    let seed_spec_str = seed_args.seed.as_deref().unwrap_or("17");
    let seed_len_specs = SeedSpec::from_str(seed_spec_str)
        .map_err(|e| anyhow!("Failed to parse seed specification: {}", e))?
        .normalize(q_seq.len())
        .map_err(|e| anyhow!("Invalid seed spec for query length: {}", e))?;

    let mut candidates = Vec::new();
    let q_len = q_seq.len();

    let (start, end, mi_len) = seed_len_specs;
    let start0 = start - 1;
    let end0 = end - 1;

    if start0 + mi_len > q_len {
        return Ok(candidates);
    }

    let last_start = end0.saturating_sub(mi_len - 1);

    for q_pos in start0..=last_start {
        // Max seed length from this position
        let max_seed_len = (end0 + 1).saturating_sub(q_pos).min(q_len - q_pos);

        for seed_len in mi_len..=max_seed_len {
            let seed_seq = &q_seq[q_pos..q_pos + seed_len];
            if seed_seq.contains(&b'N') || seed_seq.contains(&b'n') {
                continue;
            }

            // Index search (RC of seed)
            let seed_rc = seed_seq.reverse_complement_rna();
            // find_candidates returns (target_idx, target_start, is_antisense)

            let hits = index.find_candidates(&seed_rc, seed_args.wobble);

            for (target_idx, target_start, is_antisense) in hits {
                candidates.push(SeedCandidate {
                    query_pos: q_pos,
                    target_idx,
                    target_start,
                    len: seed_len,
                    is_antisense,
                });
            }
        }
    }

    Ok(candidates)
}

fn process_candidate(
    q_id: &str,
    q_seq: &[u8],
    index: &SaIndex<'_>,
    candidate: &SeedCandidate,
    opts: &ExtendArgs,
    ctx: &mut DpContext,
) -> Option<SearchHit> {
    let t_idx = candidate.target_idx;
    let t_seq_cow = if candidate.is_antisense {
        std::borrow::Cow::Owned(index.get_sequence_rc(t_idx))
    } else {
        std::borrow::Cow::Borrowed(index.get_sequence(t_idx))
    };
    let t_seq = &t_seq_cow;
    let t_start_idx = candidate.target_start;
    let seed_len = candidate.len;
    let q_pos = candidate.query_pos;

    if t_start_idx + seed_len > t_seq.len() {
        return None;
    }

    // Call extend_seed (or essentially reproduce its valuable logic).
    // Call extend_seed (or essentially reproduce its valuable logic).
    let extension_result = extend_seed(ctx, q_seq, t_seq, candidate, opts)?;

    let ext = extension_result;
    let score = ext.score; // RESTORED

    if score > opts.delta_g {
        return None;
    }

    let l_q = ext.l_q;
    let l_t = ext.l_t;
    let r_q = ext.r_q;
    let r_t = ext.r_t;
    let seed_q = q_pos;
    let seed_t = t_start_idx;

    // Use Helper to get aligned sequence strings consistent with the trace
    let (_final_qs, final_ts) = reconstruct_seqs_from_trace(
        q_seq,
        t_seq,
        q_pos,
        t_start_idx,
        seed_len,
        &ext.l_trace,
        &ext.r_trace,
        l_q,
        l_t,
        r_q,
        r_t,
    );

    let full_fp = ext.interaction;
    // full_ts from helper is the aligned target string
    let full_ts = final_ts;

    // Normalize T -> U for output strings seems to be done in `print_search_hit`.
    // Let's store raw T strings in `SearchHit` and normalize at print time?
    // Or normalize here. `SearchHit` usually implies "ready to use".
    // I'll leave them as is (DNA T) and normalize in print or here.
    // The `SearchHit` definition has `interaction`.

    // Coordinates
    let final_q_start = seed_q - l_q;
    let final_q_end = (seed_q + seed_len - 1) + r_q;
    let final_t_start = seed_t - r_t; // 0-based index in t_seq
    let final_t_end = (seed_t + seed_len - 1) + l_t; // 0-based index in t_seq

    // Output Coordinates (Strand Aware)
    // If it's Forward search: t_start .. t_end relative to Sequence Start.
    // If it's Reverse Complement search (`is_antisense`):
    // The `t_seq` we used was the RC of the original.
    // We need to map `final_t_start` / `final_t_end` back to the original sequence coordinates.
    // Let N = original len.
    // RC index i corresponds to forward index (N - 1 - i).
    // So Range [start, end] in RC maps to [N - 1 - end, N - 1 - start] in Forward?
    // Let's verify.
    // RC: 0 1 2 ... (N-1)
    // Fwd: (N-1) ... 2 1 0
    // Yes.

    let original_len = index.get_sequence_len(candidate.target_idx);
    let (out_t_start, out_t_end, strand_char) = if candidate.is_antisense {
        // Antisense hit
        // The `final_t_start` and `final_t_end` are indices in the RC sequence.
        // Convert to forward coordinates.
        // Start in Fwd = N - 1 - End in RC.
        // End in Fwd = N - 1 - Start in RC.
        let fwd_start = original_len - 1 - final_t_end;
        let fwd_end = original_len - 1 - final_t_start;
        (fwd_start + 1, fwd_end + 1, '-')
    // Wait, C output logic:
    // If query aligns to RC of target, C reports coordinates on the... genome?
    // risearch2.x output:
    // Strand '+' usually means Sense.
    // Strand '-' usually means Antisense.
    // My `is_antisense` = true came from `search_sa_simple` on `reverse_sa`.
    // `reverse_sa` is built from RC of target.
    // If query matches RC of target, it means query binds to the "other" strand.
    // Standard conventions:
    // mRNA is the sense strand.
    // miRNA binds to mRNA.
    // So miRNA is antisense to mRNA.
    // If we find match on Forward SA (Reverse Complement of Query matching Forward Target)
    // means Query binds to Forward Target.

    // Let's check `run_search` strand logic (lines 175-187 in original file, I can't see them now).
    // Standard RIsearch output:
    // Strand is usually relative to the Target.
    // If match is on Fwd strand, strand is '+'.
    // If match is on Rev strand, strand is '-'.
    } else {
        // Forward hit
        // final_t_start is index in forward sequence.
        (final_t_start + 1, final_t_end + 1, '+')
        // Actually, let's look at `tests/c_parity.rs` output.
    };

    // Flanks
    let ctx_len = 20;
    // 5' Flank (Upstream in the sequence we searched)
    // For output, we want the flank relative to the interaction or the genome?
    // C output reports flanks from the Target Sequence.
    // If we searched RC, the `t_seq` is RC. The flanks come from RC.
    // `print_detailed_output` took `t_seq` (which was COW) and extracted.
    // So we invoke `flank` extraction on `t_seq`.

    let t_5_start = final_t_start.saturating_sub(ctx_len);
    let flank_5 = String::from_utf8_lossy(&t_seq[t_5_start..final_t_start])
        .replace('T', "U")
        .replace('t', "u")
        .chars()
        .rev()
        .collect::<String>();

    let t_3_start = final_t_end + 1;
    let t_3_end = (t_3_start + ctx_len).min(t_seq.len());
    let flank_3 = if t_3_start < t_seq.len() {
        String::from_utf8_lossy(&t_seq[t_3_start..t_3_end])
            .replace('T', "U")
            .replace('t', "u")
        // Note: `ctx_3` in `print_detailed_output` was NOT reversed?
        // Let's check line 1934 in `view_file` output.
        // "String::from_utf8_lossy... replace...". No rev().
        // Correct.
    } else {
        String::new()
    };

    Some(SearchHit {
        query_id: q_id.to_string(),
        target_id: index.get_id(t_idx).to_string(),
        query_seq: String::from_utf8_lossy(q_seq).to_string(), // Or keep raw?
        target_seq: full_ts,                                   // The alignment string
        q_start: final_q_start,
        q_end: final_q_end,
        t_start: final_t_start,
        t_end: final_t_end,
        output_t_start: out_t_start,
        output_t_end: out_t_end,
        strand: strand_char,
        energy: score,
        interaction: full_fp,
        flank_5,
        flank_3: flank_3,
    })
}

fn print_search_hit(w: &mut dyn Write, hit: &SearchHit) -> std::io::Result<()> {
    // Normalize strings for output (T->U)
    let norm_fp = &hit.interaction;
    let norm_ts = hit.target_seq.replace('T', "U").replace('t', "u");

    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}",
        hit.query_id
            .split_whitespace()
            .next()
            .unwrap_or(&hit.query_id),
        hit.q_start + 1,
        hit.q_end + 1,
        hit.target_id
            .split_whitespace()
            .next()
            .unwrap_or(&hit.target_id),
        hit.output_t_start,
        hit.output_t_end,
        hit.strand,
        hit.energy,
        norm_fp,
        norm_ts,
        hit.flank_5,
        hit.flank_3
    )
}

pub struct DpResult {
    pub score: i32,
    pub ext_q_len: usize,
    pub ext_t_len: usize,
    pub trace: String,
}

struct SeedExtension {
    score: f64,
    interaction: String,
    l_trace: String,
    r_trace: String,
    l_q: usize,
    l_t: usize,
    r_q: usize,
    r_t: usize,
}

fn reconstruct_seqs_from_trace(
    _q_seq: &[u8],
    t_seq: &[u8],
    q_seed_start: usize, // 0-based index of 5' end of seed in query
    t_seed_start: usize, // 0-based index of 5' end of seed in target (low index; target is antiparallel to query)
    seed_len: usize,
    fp_l: &str,
    fp_r: &str,
    l_q: usize,
    l_t: usize,
    _r_q: usize,
    _r_t: usize,
) -> (String, String) {
    // Reconstruct aligned target sequence from fingerprint traces.
    //
    // Coordinate system (antiparallel RNA binding):
    //   Query:  5' ----[seed]----- 3'   (indices increase left to right)
    //   Target: 3' ----[seed]----- 5'   (indices increase left to right, but binding is antiparallel)
    //
    // In the target array:
    //   - t_seed_start is the LOW index (left side of seed in array)
    //   - t_seed_start + seed_len - 1 is the HIGH index (right side of seed in array)
    //   - Due to antiparallel binding:
    //     * Query 5' (q_seed_start) pairs with Target 3' (t_seed_start + seed_len - 1)
    //     * Query 3' (q_seed_start + seed_len - 1) pairs with Target 5' (t_seed_start)
    //
    // DP extension directions:
    //   - dp_left:  Query moves 5' (index decreases), Target moves 3' (index increases)
    //   - dp_right: Query moves 3' (index increases), Target moves 5' (index decreases)
    //
    // Trace characters:
    //   - 'Q' = GapQ state: query has unpaired base, target has gap in alignment
    //   - 'T' = GapT state: target has unpaired base, query has gap in alignment
    //   - 'P'/'W'/'U' = paired/wobble/unpaired match state

    let mut aligned_ts = String::new();

    // LEFT PART (5' extension of query, 3' extension of target)
    // fp_l is ordered seed-to-far (after reversal in extend_seed).
    // We iterate REVERSE to build output in 5'->3' order (far-to-seed).
    // Start at the far-left indices and move toward the seed.
    let mut q_idx = q_seed_start - l_q;
    let mut t_idx = (t_seed_start + seed_len - 1) + l_t;

    for c in fp_l.chars().rev() {
        match c {
            'Q' => {
                // GapQ state: Query has an unpaired base (bulge in query).
                // In the aligned target string, this appears as a gap.
                aligned_ts.push('-');
                let _ = q_idx; // Suppress unused
                q_idx += 1;
            }
            'T' => {
                // GapT state: Target has an unpaired base (bulge in target).
                // In the aligned target string, we emit the target base.
                let tc = t_seq[t_idx];
                let tc_char = tc as char;
                aligned_ts.push(tc_char);
                t_idx -= 1;
            }
            _ => {
                // Match/Mismatch/Wobble
                let tc = t_seq[t_idx];
                aligned_ts.push(tc as char);
                q_idx += 1;
                t_idx -= 1;
            }
        }
    }

    // SEED PART
    // Iterate through seed bases. Due to antiparallel binding, target index decreases
    // as query index increases (from 5' to 3').
    let t_seed_end_idx = t_seed_start + seed_len - 1;
    for k in 0..seed_len {
        let t_i = t_seed_end_idx - k;
        let tc = t_seq[t_i];
        aligned_ts.push(tc as char);
    }

    // RIGHT PART (3' extension of query, 5' extension of target)
    // fp_r is ordered seed-to-far. We iterate FORWARD (seed toward 3' end).
    // Target index decreases as we move right on query.
    let mut q_idx = q_seed_start + seed_len;
    let mut t_idx = t_seed_start.wrapping_sub(1);

    for c in fp_r.chars() {
        match c {
            'Q' => {
                // GapQ: Q has base, T has gap.
                aligned_ts.push('-');
                let _ = q_idx; // Suppress unused
                q_idx += 1;
            }
            'T' => {
                // GapT: T has base, Q has gap.
                let tc = t_seq[t_idx];
                aligned_ts.push(tc as char);
                t_idx = t_idx.wrapping_sub(1);
            }
            _ => {
                let tc = t_seq[t_idx];
                aligned_ts.push(tc as char);
                q_idx += 1;
                t_idx = t_idx.wrapping_sub(1);
            }
        }
    }

    (String::new(), aligned_ts)
}

fn extend_seed(
    ctx: &mut DpContext,
    q_seq: &[u8],
    t_seq: &[u8],
    candidate: &SeedCandidate,
    opts: &ExtendArgs,
) -> Option<SeedExtension> {
    let q_pos = candidate.query_pos;
    let t_pos = candidate.target_start;
    let len = candidate.len;

    // MAXIMALITY CHECK
    // Skip non-maximal seeds: if the seed can be extended by a valid base pair
    // on either end, it's a sub-seed of a longer match and will be found later.

    // 1. Left extendable? Check if q[q_pos-1] pairs with t[t_pos+len]
    if q_pos > 0 && t_pos + len < t_seq.len() {
        let q_prev = Base::from_byte(q_seq[q_pos - 1]).idx();
        let t_next = Base::from_byte(t_seq[t_pos + len]).idx();
        if PAIR_MAT[q_prev][t_next] != 0 {
            return None;
        }
    }

    // 2. Right extendable? Check if q[q_pos+len] pairs with t[t_pos-1]
    if q_pos + len < q_seq.len() && t_pos > 0 {
        let q_next = Base::from_byte(q_seq[q_pos + len]).idx();
        let t_prev = Base::from_byte(t_seq[t_pos - 1]).idx();
        if PAIR_MAT[q_next][t_prev] != 0 {
            return None;
        }
    }

    let mut seed_energy = 0.0;

    // Seed energy calculation (antiparallel binding).
    // q[q_pos + k] pairs with t[t_pos + len - 1 - k]

    let t_match_end = t_pos + len - 1;
    let mut seed_int_str = String::with_capacity(len);

    for k in 0..len {
        let q_idx = q_pos + k;
        let t_idx = t_match_end - k;

        if k < len.saturating_sub(1) {
            let q_b1 = Base::from_byte(q_seq[q_idx]).idx();
            let q_b2 = Base::from_byte(q_seq[q_idx + 1]).idx();
            let t_b1 = Base::from_byte(t_seq[t_idx]).idx();
            let t_b2 = Base::from_byte(t_seq[t_match_end - (k + 1)]).idx();
            seed_energy += DSM_T04_POS[q_b1][q_b2][t_b1][t_b2] as f64;
        }

        // Build interaction string
        let qc = q_seq[q_idx];
        let tc = t_seq[t_idx];
        seed_int_str.push(Base::from_byte(qc).pairing_class(Base::from_byte(tc)));
    }

    // DEBUG: Trace internal_mismatch_20nt
    if opts.max_extension == 20 && seed_int_str.starts_with('U') {
        println!(
            "DEBUG_TRACE: seed_int_str={} q_pos={} t_pos={}",
            seed_int_str, q_pos, t_pos
        );
        for k in 0..len {
            let q_idx = q_pos + k;
            let t_idx = t_match_end - k;
            let qc = q_seq[q_idx];
            let tc = t_seq[t_idx];
            let p = Base::from_byte(qc).pairing_class(Base::from_byte(tc));
            println!(
                "  k={} q[{}]={} t[{}]={} pair={}",
                k, q_idx, qc as char, t_idx, tc as char, p
            );
        }
    }
    // Seed interaction string should match C output; omit any extra markers

    let max_ext = opts.max_extension as usize;
    let safe_ext = max_ext.min(MAX_DP_EXT);

    // DP Left: Extend Query Left (5'), Target Right (3')
    // Start at: q_pos (5' Q), t_match_end (3' T)
    let left_res = dp_left(ctx, q_seq, t_seq, q_pos, t_match_end, safe_ext);

    // DP Right: Extend Query Right (3'), Target Left (5')
    // Start at: q_pos + len - 1 (3' Q), t_pos (5' T)
    let right_res = dp_right(ctx, q_seq, t_seq, q_pos + len - 1, t_pos, safe_ext);

    let final_score =
        (seed_energy + left_res.score as f64 + right_res.score as f64 - 559.0) / -100.0;

    // Debug for first few calls in a clean way
    static DEBUG_EXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let c = DEBUG_EXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);

    // DEBUG: Trace internal_mismatch_20nt
    if opts.max_extension == 20 {
        let l_trace_str: String = left_res.trace.chars().rev().collect();
        let r_trace_str: String = right_res.trace.chars().rev().collect();
        let full_int = format!("{}{}{}", l_trace_str, seed_int_str, r_trace_str);

        println!(
            "DEBUG_TRACE: full={} seed={} l_trace={} r_trace={} score={:.2}",
            full_int, seed_int_str, left_res.trace, right_res.trace, final_score
        );

        // If we have a U at start, dump the sequences
        if full_int.starts_with('U') {
            println!("DEBUG_TRACE:   -> HIT INTERESTING CASE: U at start");
            println!(
                "DEBUG_TRACE:   -> q_pos={} t_pos={} seed_len={}",
                q_pos, t_pos, len
            );
            println!(
                "DEBUG_TRACE:   -> LENGTHS: full={} seed={} l_trace={} r_trace={}",
                full_int.len(),
                seed_int_str.len(),
                left_res.trace.len(),
                right_res.trace.len()
            );

            // Analyze l_trace bases (length 1 for now)
            let q_idx = q_pos.wrapping_sub(1);
            let t_idx = t_pos + len;
            let qc = if q_idx < q_seq.len() {
                q_seq[q_idx] as char
            } else {
                '?'
            };
            let tc = if t_idx < t_seq.len() {
                t_seq[t_idx] as char
            } else {
                '?'
            };
            let p_char = Base::from_byte(qc as u8).pairing_class(Base::from_byte(tc as u8));
            println!(
                "DEBUG_TRACE:   -> L_EXT check: q[{}]={} vs t[{}]={} -> PairClass={}",
                q_idx, qc, t_idx, tc, p_char
            );
        }
    }

    if final_score < -20.0 {
        println!(
            "C_DEBUG: extend_seed: seed_energy={:.0}, l_score={}, r_score={}, raw_total={:.0}, final={:.2}",
            seed_energy,
            left_res.score,
            right_res.score,
            seed_energy + left_res.score as f64 + right_res.score as f64,
            final_score
        );
    }

    let l_trace_str: String = left_res.trace.chars().collect();
    let r_trace_str: String = right_res.trace.chars().rev().collect();

    let full_interaction = format!("{}{}{}", l_trace_str, seed_int_str, r_trace_str);

    Some(SeedExtension {
        score: final_score,
        interaction: full_interaction.clone(),
        l_trace: l_trace_str,
        r_trace: r_trace_str,
        l_q: left_res.ext_q_len,
        l_t: left_res.ext_t_len,
        r_q: right_res.ext_q_len,
        r_t: right_res.ext_t_len,
    })
}

fn dp_left(
    ctx: &mut DpContext,
    q_seq: &[u8],
    t_seq: &[u8],
    q_start: usize,
    t_start: usize,
    max_ext: usize,
) -> DpResult {
    // DP Left: Extend Query toward 5' (decreasing index), Target toward 3' (increasing index).
    // q_start: index of first base in seed (query 5' end of seed)
    // t_start: index of last base in seed (target 3' end of seed, due to antiparallel binding)
    //
    // q_len: number of query positions available for extension (indices q_start down to 0)
    // t_len: number of target positions available for extension (indices t_start+1 to end)
    let q_len = (q_start + 1).min(max_ext);
    let t_len = (t_seq.len() - t_start - 1).min(max_ext);

    let s_mat = &DSM_T04_POS;

    // Index mapping functions for the DP:
    // q_char(i) -> nucleotide at q_start - i (query extends left, index decreases)
    // t_comp(j) -> nucleotide at t_start + j (target extends right, index increases)
    let q_char = |i: usize| {
        if i > q_start {
            0
        } else {
            Base::from_byte(q_seq[q_start - i]).idx()
        }
    };
    let t_comp = |j: usize| {
        if t_start + j >= t_seq.len() {
            0
        } else {
            Base::from_byte(t_seq[t_start + j]).idx()
        }
    };

    // Initial score: terminal penalty for the seed boundary base pair.
    // This matches C code: best_e = (*S)[GAP][Q(0)][GAP][T(0)]

    let mut best_e = s_mat[GAP_IDX][q_char(0)][GAP_IDX][t_comp(0)] as i32;
    let mut best_i = 0;
    let mut best_j = 0;

    // C: if (lq <= 1 || lt <= 1) return best_e;
    if q_len <= 1 || t_len <= 1 {
        return DpResult {
            score: best_e,
            ext_q_len: best_i,
            ext_t_len: best_j,
            trace: String::new(),
        };
    }

    // Resize context grids
    ctx.resize(t_len + 1, q_len + 1);

    // Create mutable aliases for easier access
    // Rust allows borrowing disjoint fields mutably
    let DpContext {
        m,
        bq,
        bt,
        tb_m,
        tb_bq,
        tb_bt,
    } = ctx;

    m.set(0, 0, Some(0));

    // Init (0,1), (1,0), (1,1)
    bt.set(
        0,
        1,
        Some(s_mat[GAP_IDX][q_char(0)][t_comp(1)][t_comp(0)] as i32),
    );
    bq.set(
        1,
        0,
        Some(s_mat[q_char(1)][q_char(0)][GAP_IDX][t_comp(0)] as i32),
    );
    m.set(
        1,
        1,
        Some(s_mat[q_char(1)][q_char(0)][t_comp(1)][t_comp(0)] as i32),
    );

    if let Some(m11) = m.get(1, 1) {
        let val = m11 + s_mat[GAP_IDX][q_char(1)][GAP_IDX][t_comp(1)] as i32;
        if val > best_e {
            best_e = val;
            best_i = 1;
            best_j = 1;
        }
    }

    // Row 0 (j from 2 to t_len-1)
    for j in 2..t_len {
        if let Some(prev_bt) = bt.get(0, j - 1) {
            bt.set(
                0,
                j,
                Some(prev_bt + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32),
            );
            tb_bt.set(0, j, TraceStep::GapT);

            if q_len >= 1 {
                let new_m = prev_bt + s_mat[q_char(1)][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32;
                m.set(1, j, Some(new_m));
                tb_m.set(1, j, TraceStep::GapT);
                let val = new_m + s_mat[GAP_IDX][q_char(1)][GAP_IDX][t_comp(j)] as i32;
                if val > best_e {
                    best_e = val;
                    best_i = 1;
                    best_j = j;
                }
            }
        }
    }

    // Col 0 (i from 2 to q_len-1)
    for i in 2..q_len {
        if let Some(prev_bq) = bq.get(i - 1, 0) {
            bq.set(
                i,
                0,
                Some(prev_bq + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32),
            );
            tb_bq.set(i, 0, TraceStep::GapQ);

            if t_len >= 1 {
                let new_m = prev_bq + s_mat[q_char(i)][q_char(i - 1)][t_comp(1)][GAP_IDX] as i32;
                m.set(i, 1, Some(new_m));
                tb_m.set(i, 1, TraceStep::GapQ);
                let val = new_m + s_mat[GAP_IDX][q_char(i)][GAP_IDX][t_comp(1)] as i32;
                if val > best_e {
                    best_e = val;
                    best_i = i;
                    best_j = 1;
                }
            }
        }
    }

    // 2,2 Init (needs at least length 3 to have index 2)
    if q_len >= 3 && t_len >= 3 {
        if let Some(m11) = m.get(1, 1) {
            bt.set(
                1,
                2,
                Some(m11 + s_mat[GAP_IDX][q_char(1)][t_comp(2)][t_comp(1)] as i32),
            );
            tb_bt.set(1, 2, TraceStep::Match);

            bq.set(
                2,
                1,
                Some(m11 + s_mat[q_char(2)][q_char(1)][GAP_IDX][t_comp(1)] as i32),
            );
            tb_bq.set(2, 1, TraceStep::Match);

            let m22 = m11 + s_mat[q_char(2)][q_char(1)][t_comp(2)][t_comp(1)] as i32;
            m.set(2, 2, Some(m22));
            tb_m.set(2, 2, TraceStep::Match);

            let val = m22 + s_mat[GAP_IDX][q_char(2)][GAP_IDX][t_comp(2)] as i32;
            if val > best_e {
                best_e = val;
                best_i = 2;
                best_j = 2;
            }
        }
        if let Some(m12) = m.get(1, 2) {
            bq.set(
                2,
                2,
                Some(m12 + s_mat[q_char(2)][q_char(1)][GAP_IDX][t_comp(2)] as i32),
            );
            tb_bq.set(2, 2, TraceStep::Match);
        }
        if let Some(m21) = m.get(2, 1) {
            bt.set(
                2,
                2,
                Some(m21 + s_mat[GAP_IDX][q_char(2)][t_comp(2)][t_comp(1)] as i32),
            );
            tb_bt.set(2, 2, TraceStep::Match);
        }
    }

    // Main DP
    for i in 2..q_len {
        for j in 2..t_len {
            if i == 2 && j == 2 {
                continue;
            }

            // Calc M[i,j] - pick best of three sources
            let s_mm = m
                .get(i - 1, j - 1)
                .map(|v| v + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][t_comp(j - 1)] as i32);
            let s_mq = bq
                .get(i - 1, j - 1)
                .map(|v| v + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][GAP_IDX] as i32);
            let s_mt = bt
                .get(i - 1, j - 1)
                .map(|v| v + s_mat[q_char(i)][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32);

            // Find best value and corresponding traceback step
            let (val_m, step_m) = [
                (s_mm, TraceStep::Match),
                (s_mq, TraceStep::GapQ),
                (s_mt, TraceStep::GapT),
            ]
            .into_iter()
            .filter_map(|(opt, step)| opt.map(|v| (v, step)))
            .max_by_key(|(v, _)| *v)
            .unwrap_or((i32::MIN, TraceStep::Stop));

            let val_m = if val_m == i32::MIN { None } else { Some(val_m) };
            m.set(i, j, val_m);
            tb_m.set(i, j, step_m);

            if let Some(v) = val_m {
                let curr_e = v + s_mat[GAP_IDX][q_char(i)][GAP_IDX][t_comp(j)] as i32;
                if curr_e > best_e {
                    best_e = curr_e;
                    best_i = i;
                    best_j = j;
                }
            }

            // Calc Bq[i,j]
            if i > 2 || (i == 2 && j > 2) {
                let s_qm = m
                    .get(i - 1, j)
                    .map(|v| v + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][t_comp(j)] as i32);
                let s_qq = bq
                    .get(i - 1, j)
                    .map(|v| v + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32);

                // Priority to Match (opening) if tie
                match (s_qq, s_qm) {
                    (Some(qq), Some(qm)) if qq > qm => {
                        bq.set(i, j, Some(qq));
                        tb_bq.set(i, j, TraceStep::GapQ);
                    }
                    (_, Some(qm)) => {
                        bq.set(i, j, Some(qm));
                        tb_bq.set(i, j, TraceStep::Match);
                    }
                    (Some(qq), None) => {
                        bq.set(i, j, Some(qq));
                        tb_bq.set(i, j, TraceStep::GapQ);
                    }
                    _ => {}
                }
            }

            // Calc Bt[i,j]
            if j > 2 || (j == 2 && i > 2) {
                let s_tm = m
                    .get(i, j - 1)
                    .map(|v| v + s_mat[GAP_IDX][q_char(i)][t_comp(j)][t_comp(j - 1)] as i32);
                let s_tt = bt
                    .get(i, j - 1)
                    .map(|v| v + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32);

                // Priority to Match (opening) if tie
                match (s_tt, s_tm) {
                    (Some(tt), Some(tm)) if tt > tm => {
                        bt.set(i, j, Some(tt));
                        tb_bt.set(i, j, TraceStep::GapT);
                    }
                    (_, Some(tm)) => {
                        bt.set(i, j, Some(tm));
                        tb_bt.set(i, j, TraceStep::Match);
                    }
                    (Some(tt), None) => {
                        bt.set(i, j, Some(tt));
                        tb_bt.set(i, j, TraceStep::GapT);
                    }
                    _ => {}
                }
            }
        }
    }

    // Backtracking: trace from best position back to seed boundary.
    // Start in Match state since best_e is always computed from M matrix.
    let mut i = best_i;
    let mut j = best_j;
    let mut fp = String::new();
    let mut state = DpState::Match;

    while i > 0 || j > 0 {
        let _qc = if i <= q_len { q_seq[q_start - i] } else { b'N' };

        let qc_byte = if i <= q_start {
            Base::from_byte(q_seq[q_start - i]).idx()
        } else {
            GAP_IDX
        };
        let t_char_idx = if t_start + j < t_seq.len() {
            Base::from_byte(t_seq[t_start + j]).idx()
        } else {
            GAP_IDX
        };

        match state {
            DpState::Match => {
                // Current state is Match (M[i,j])
                // We emit the character pair corresponding to this match/mismatch
                let p_char = Base::from_idx(qc_byte).pairing_class(Base::from_idx(t_char_idx));
                if p_char == 'U' && i == 1 {
                    // Trace the specific U mismatch at boundary
                    println!(
                        "DEBUG_TRACE: DEBUG_DP_LEFT: i={} j={} qc={} tc={} pair={}",
                        i,
                        j,
                        Base::from_idx(qc_byte).as_char(),
                        Base::from_idx(t_char_idx).as_char(),
                        p_char
                    );
                }
                fp.push(p_char);

                if i == 0 || j == 0 {
                    // Should not happen for Match state unless logic is wrong
                    // In rigorous check: if i=0, we can't be in Match.
                    // But if we are, we break.
                    break;
                }

                let step = tb_m.get(i, j);
                i -= 1;
                j -= 1;

                match step {
                    TraceStep::Stop => break,
                    TraceStep::Match => state = DpState::Match,
                    TraceStep::GapQ => state = DpState::GapQ,
                    TraceStep::GapT => state = DpState::GapT,
                }
            }
            DpState::GapQ => {
                // Current state is GapQ (Bq[i,j])
                fp.push('Q');

                let step = tb_bq.get(i, j);
                if i > 0 {
                    i -= 1;
                } else {
                    break;
                }

                match step {
                    TraceStep::Stop => break,
                    TraceStep::Match => state = DpState::Match,
                    TraceStep::GapQ => state = DpState::GapQ,
                    TraceStep::GapT => state = DpState::GapT, // Should not occur
                }
            }
            DpState::GapT => {
                // Current state is GapT (Bt[i,j])
                fp.push('T');

                let step = tb_bt.get(i, j);
                if j > 0 {
                    j -= 1;
                } else {
                    break;
                }

                match step {
                    TraceStep::Stop => break,
                    TraceStep::Match => state = DpState::Match,
                    TraceStep::GapT => state = DpState::GapT,
                    TraceStep::GapQ => state = DpState::GapQ, // Should not occur
                }
            }
        }
    }

    DpResult {
        score: best_e,
        ext_q_len: best_i,
        ext_t_len: best_j,
        trace: fp,
    }
}

fn dp_right(
    ctx: &mut DpContext,
    q_seq: &[u8],
    t_seq: &[u8],
    q_end: usize,
    t_end: usize,
    max_ext: usize,
) -> DpResult {
    // DP Right: Extend Query toward 3' (increasing index), Target toward 5' (decreasing index).
    // q_end: index of last base in seed (query 3' end of seed)
    // t_end: index of first base in seed (target 5' end of seed, due to antiparallel binding)
    //
    // q_len: number of query positions available (indices q_end to end of sequence)
    // t_len: number of target positions available (indices t_end down to 0)
    let q_len = (q_seq.len() - q_end).min(max_ext);
    let t_len = (t_end + 1).min(max_ext);

    let s_mat = &DSM_T04_POS;

    // Index mapping functions for the DP:
    // q_char(i) -> nucleotide at q_end + i (query extends right, index increases)
    // t_comp(j) -> nucleotide at t_end - j (target extends left, index decreases)

    let q_char = |i: usize| {
        if q_end + i >= q_seq.len() {
            0
        } else {
            Base::from_byte(q_seq[q_end + i]).idx()
        }
    };
    let t_comp = |j: usize| {
        if j > t_end {
            0
        } else {
            Base::from_byte(t_seq[t_end - j]).idx()
        }
    };

    // Initial score from seed boundary
    // C DP_right uses Q(0)->Gap, T(0)->Gap logic ([Q][Gap][T][Gap])
    let mut best_e = s_mat[q_char(0)][GAP_IDX][t_comp(0)][GAP_IDX] as i32;
    let mut best_i = 0;
    let mut best_j = 0;

    // C: if (lq <= 1 || lt <= 1) return best_e;
    if q_len <= 1 || t_len <= 1 {
        return DpResult {
            score: best_e,
            ext_q_len: best_i,
            ext_t_len: best_j,
            trace: String::new(),
        };
    }

    // Resize context grids
    ctx.resize(t_len + 1, q_len + 1);

    // Create mutable aliases
    let DpContext {
        m,
        bq,
        bt,
        tb_m,
        tb_bq,
        tb_bt,
    } = ctx;

    m.set(0, 0, Some(0));

    // Init (0,1) Bt, (1,0) Bq, (1,1) M
    bt.set(
        0,
        1,
        Some(s_mat[q_char(0)][GAP_IDX][t_comp(0)][t_comp(1)] as i32),
    );
    bq.set(
        1,
        0,
        Some(s_mat[q_char(0)][q_char(1)][t_comp(0)][GAP_IDX] as i32),
    );
    m.set(
        1,
        1,
        Some(s_mat[q_char(0)][q_char(1)][t_comp(0)][t_comp(1)] as i32),
    );

    if let Some(m11) = m.get(1, 1) {
        let val = m11 + s_mat[q_char(1)][GAP_IDX][t_comp(1)][GAP_IDX] as i32;
        if val > best_e {
            best_e = val;
            best_i = 1;
            best_j = 1;
        }
    }

    // Row 0 (j from 2 to t_len-1)
    for j in 2..t_len {
        if let Some(prev_bt) = bt.get(0, j - 1) {
            bt.set(
                0,
                j,
                Some(prev_bt + s_mat[GAP_IDX][GAP_IDX][t_comp(j - 1)][t_comp(j)] as i32),
            );
            tb_bt.set(0, j, TraceStep::GapT);

            if q_len >= 1 {
                let new_m = prev_bt + s_mat[GAP_IDX][q_char(1)][t_comp(j - 1)][t_comp(j)] as i32;
                m.set(1, j, Some(new_m));
                tb_m.set(1, j, TraceStep::GapT);
                let val = new_m + s_mat[q_char(1)][GAP_IDX][t_comp(j)][GAP_IDX] as i32;
                if val > best_e {
                    best_e = val;
                    best_i = 1;
                    best_j = j;
                }
            }
        }
    }

    // Col 0 (i from 2 to q_len-1)
    for i in 2..q_len {
        if let Some(prev_bq) = bq.get(i - 1, 0) {
            bq.set(
                i,
                0,
                Some(prev_bq + s_mat[q_char(i - 1)][q_char(i)][GAP_IDX][GAP_IDX] as i32),
            );
            tb_bq.set(i, 0, TraceStep::GapQ);

            if t_len >= 1 {
                let new_m = prev_bq + s_mat[q_char(i - 1)][q_char(i)][GAP_IDX][t_comp(1)] as i32;
                m.set(i, 1, Some(new_m));
                tb_m.set(i, 1, TraceStep::GapQ);
                let val = new_m + s_mat[q_char(i)][GAP_IDX][t_comp(1)][GAP_IDX] as i32;
                if val > best_e {
                    best_e = val;
                    best_i = i;
                    best_j = 1;
                }
            }
        }
    }

    // 2,2 Init (needs at least length 3 to have index 2)
    if q_len >= 3 && t_len >= 3 {
        if let Some(m11) = m.get(1, 1) {
            bt.set(
                1,
                2,
                Some(m11 + s_mat[q_char(1)][GAP_IDX][t_comp(1)][t_comp(2)] as i32),
            );
            tb_bt.set(1, 2, TraceStep::Match);

            bq.set(
                2,
                1,
                Some(m11 + s_mat[q_char(1)][q_char(2)][t_comp(1)][GAP_IDX] as i32),
            );
            tb_bq.set(2, 1, TraceStep::Match);

            let m22 = m11 + s_mat[q_char(1)][q_char(2)][t_comp(1)][t_comp(2)] as i32;
            m.set(2, 2, Some(m22));
            tb_m.set(2, 2, TraceStep::Match);

            let val = m22 + s_mat[q_char(2)][GAP_IDX][t_comp(2)][GAP_IDX] as i32;
            if val > best_e {
                best_e = val;
                best_i = 2;
                best_j = 2;
            }
        }
        if let Some(m12) = m.get(1, 2) {
            bq.set(
                2,
                2,
                Some(m12 + s_mat[q_char(1)][q_char(2)][t_comp(2)][GAP_IDX] as i32),
            );
            tb_bq.set(2, 2, TraceStep::Match);
        }
        if let Some(m21) = m.get(2, 1) {
            bt.set(
                2,
                2,
                Some(m21 + s_mat[q_char(2)][GAP_IDX][t_comp(1)][t_comp(2)] as i32),
            );
            tb_bt.set(2, 2, TraceStep::Match);
        }
    }

    // Main DP
    for i in 2..q_len {
        for j in 2..t_len {
            if i == 2 && j == 2 {
                continue;
            }

            // Calc M[i,j] - pick best of three sources
            let s_mm = m
                .get(i - 1, j - 1)
                .map(|v| v + s_mat[q_char(i - 1)][q_char(i)][t_comp(j - 1)][t_comp(j)] as i32);
            let s_mq = bq
                .get(i - 1, j - 1)
                .map(|v| v + s_mat[q_char(i - 1)][q_char(i)][GAP_IDX][t_comp(j)] as i32);
            let s_mt = bt
                .get(i - 1, j - 1)
                .map(|v| v + s_mat[GAP_IDX][q_char(i)][t_comp(j - 1)][t_comp(j)] as i32);

            // Find best value and corresponding traceback step
            let (val_m, step_m) = [
                (s_mm, TraceStep::Match),
                (s_mq, TraceStep::GapQ),
                (s_mt, TraceStep::GapT),
            ]
            .into_iter()
            .filter_map(|(opt, step)| opt.map(|v| (v, step)))
            .max_by_key(|(v, _)| *v)
            .unwrap_or((i32::MIN, TraceStep::Stop));

            let val_m = if val_m == i32::MIN { None } else { Some(val_m) };
            m.set(i, j, val_m);
            tb_m.set(i, j, step_m);

            if let Some(v) = val_m {
                let term = s_mat[q_char(i)][GAP_IDX][t_comp(j)][GAP_IDX] as i32;
                if v + term > best_e {
                    best_e = v + term;
                    best_i = i;
                    best_j = j;
                }
            }

            // Calc Bq[i,j]
            if i > 2 || (i == 2 && j > 2) {
                let s_qm = m
                    .get(i - 1, j)
                    .map(|v| v + s_mat[q_char(i - 1)][q_char(i)][t_comp(j)][GAP_IDX] as i32);
                let s_qq = bq
                    .get(i - 1, j)
                    .map(|v| v + s_mat[q_char(i - 1)][q_char(i)][GAP_IDX][GAP_IDX] as i32);

                match (s_qq, s_qm) {
                    (Some(qq), Some(qm)) if qq > qm => {
                        bq.set(i, j, Some(qq));
                        tb_bq.set(i, j, TraceStep::GapQ);
                    }
                    (_, Some(qm)) => {
                        bq.set(i, j, Some(qm));
                        tb_bq.set(i, j, TraceStep::Match);
                    }
                    (Some(qq), None) => {
                        bq.set(i, j, Some(qq));
                        tb_bq.set(i, j, TraceStep::GapQ);
                    }
                    _ => {}
                }
            }

            // Calc Bt[i,j]
            if j > 2 || (j == 2 && i > 2) {
                let s_tm = m
                    .get(i, j - 1)
                    .map(|v| v + s_mat[q_char(i)][GAP_IDX][t_comp(j - 1)][t_comp(j)] as i32);
                let s_tt = bt
                    .get(i, j - 1)
                    .map(|v| v + s_mat[GAP_IDX][GAP_IDX][t_comp(j - 1)][t_comp(j)] as i32);

                match (s_tt, s_tm) {
                    (Some(tt), Some(tm)) if tt > tm => {
                        bt.set(i, j, Some(tt));
                        tb_bt.set(i, j, TraceStep::GapT);
                    }
                    (_, Some(tm)) => {
                        bt.set(i, j, Some(tm));
                        tb_bt.set(i, j, TraceStep::Match);
                    }
                    (Some(tt), None) => {
                        bt.set(i, j, Some(tt));
                        tb_bt.set(i, j, TraceStep::GapT);
                    }
                    _ => {}
                }
            }
        }
    }

    // Backtracking: trace from best position back to seed boundary.
    // Start in Match state since best_e is always computed from M matrix.
    let mut i = best_i;
    let mut j = best_j;
    let mut fp = String::new();
    let mut state = DpState::Match;

    while i > 0 || j > 0 {
        let qc_byte = if q_end + i < q_seq.len() {
            Base::from_byte(q_seq[q_end + i]).idx()
        } else {
            GAP_IDX
        };
        let t_char_idx = if t_end >= j {
            Base::from_byte(t_seq[t_end - j]).idx()
        } else {
            GAP_IDX
        };

        match state {
            DpState::Match => {
                // Current state is Match (M[i,j])
                fp.push(Base::from_idx(qc_byte).pairing_class(Base::from_idx(t_char_idx)));

                if i == 0 || j == 0 {
                    break;
                }

                let step = tb_m.get(i, j);
                i -= 1;
                j -= 1;

                match step {
                    TraceStep::Stop => break,
                    TraceStep::Match => state = DpState::Match,
                    TraceStep::GapQ => state = DpState::GapQ,
                    TraceStep::GapT => state = DpState::GapT,
                }
            }
            DpState::GapQ => {
                // Current state is GapQ (Bq[i,j])
                fp.push('Q');

                let step = tb_bq.get(i, j);
                if i > 0 {
                    i -= 1;
                } else {
                    break;
                }
                // j stays same

                match step {
                    TraceStep::Stop => break,
                    TraceStep::Match => state = DpState::Match,
                    TraceStep::GapQ => state = DpState::GapQ,
                    TraceStep::GapT => state = DpState::GapT, // Should not occur
                }
            }
            DpState::GapT => {
                // Current state is GapT (Bt[i,j])
                fp.push('T');

                let step = tb_bt.get(i, j);
                if j > 0 {
                    j -= 1;
                } else {
                    break;
                }
                // i stays same

                match step {
                    TraceStep::Stop => break,
                    TraceStep::Match => state = DpState::Match,
                    TraceStep::GapT => state = DpState::GapT,
                    TraceStep::GapQ => state = DpState::GapQ, // Should not occur
                }
            }
        }
    }

    DpResult {
        score: best_e,
        ext_q_len: best_i,
        ext_t_len: best_j,
        trace: fp,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]

    fn test_dp_left_returns_score() {
        // Verifies dp_left runs without panicking and returns a score tuple
        let mut ctx = DpContext::new(10, 10);
        let q = b"AAAA";
        let t = b"UUUU";
        let res = dp_left(&mut ctx, q, t, 3, 3, 10);
        // The actual score depends on the scoring matrix and sequence alignment
        // For now, just verify the function returns without panicking
        println!(
            "dp_left score: {}, i: {}, j: {}",
            res.score, res.ext_q_len, res.ext_t_len
        );
    }

    #[test]
    fn test_dp_right_returns_score() {
        // Verifies dp_right runs without panicking and returns a score tuple
        let mut ctx = DpContext::new(10, 10);
        let q = b"AAAA";
        let t = b"UUUU";
        let res = dp_right(&mut ctx, q, t, 0, 0, 10);
        println!(
            "dp_right score: {}, i: {}, j: {}",
            res.score, res.ext_q_len, res.ext_t_len
        );
    }
}
