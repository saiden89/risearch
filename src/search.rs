use anyhow::{Context, Result, anyhow};
use log::{debug, info, trace, warn};
use std::io::Write;
use std::path::Path;
use std::str::FromStr;

use crate::dsm::{Base, DSM_T04_POS, PAIR_MAT};
use crate::sa::IndexFile;
use crate::seed::SeedSpec;
use crate::{SearchArgs, SeedPairing, Strand};

use std::collections::HashMap;

const MAX_DP_EXT: usize = 30;
const GAP_IDX: usize = Base::Gap as usize;

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

    // Deduplication (deduplicate_hits)
    DedupExactMatch,
    DedupContainedByShadow,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pairing {
    Match(u8, u8),    // e.g. (b'G', b'C')
    Wobble(u8, u8),   // e.g. (b'G', b'U')
    Mismatch(u8, u8), // e.g. (b'A', b'A')
    GapQuery(u8),     // Gap in Query, Base in Target (stored as u8)
    GapTarget(u8),    // Gap in Target, Base in Query (stored as u8)
}

impl Pairing {
    pub fn from_bases(q_byte: u8, t_byte: u8) -> Self {
        let q = Base::from_byte(q_byte);
        let t = Base::from_byte(t_byte);
        // Base::pairing_class logic check:
        // G-C -> |, G-U -> :, G-A -> ' ', A-U -> |
        // We match logic of Base::pairing_class but return Enum
        match (q, t) {
            (Base::G, Base::C) | (Base::C, Base::G) | (Base::A, Base::U) | (Base::U, Base::A) => {
                Pairing::Match(q_byte, t_byte)
            }

            (Base::G, Base::U) | (Base::U, Base::G) => Pairing::Wobble(q_byte, t_byte),

            _ => Pairing::Mismatch(q_byte, t_byte),
        }
    }

    pub fn to_char(&self) -> char {
        match self {
            Pairing::Match(_, _) => 'P',
            Pairing::Wobble(_, _) => 'W',
            Pairing::Mismatch(_, _) => 'U',
            Pairing::GapQuery(_) => 'T', // Gap in Query = Target Bulge ('T')
            Pairing::GapTarget(_) => 'Q', // Gap in Target = Query Bulge ('Q')
        }
    }

    pub fn target_char(&self) -> char {
        match self {
            Pairing::Match(_, t)
            | Pairing::Wobble(_, t)
            | Pairing::Mismatch(_, t)
            | Pairing::GapQuery(t) => *t as char,
            Pairing::GapTarget(_) => '-',
        }
    }
}

use std::ops::Range;

/// Represents a full biological alignment between Query and Target.
/// Stores the sequence of interactions and metadata about the seed location.
#[derive(Debug, Clone)]
pub struct Alignment {
    /// The complete sequence of pairing steps (5' -> 3' of Query).
    steps: Vec<Pairing>,

    /// The range of indices in `steps` that corresponds to the initial Seed match.
    /// This allows easy extraction of the "core" interaction vs extensions.
    seed_range: Range<usize>,
}

impl Alignment {
    /// Constructor from the three phases of extension.
    /// This fits naturally into `extend_seed` which generates these 3 parts.
    pub fn new(left: Vec<Pairing>, seed: Vec<Pairing>, right: Vec<Pairing>) -> Self {
        let left_len = left.len();
        let seed_len = seed.len();

        let mut steps = Vec::with_capacity(left_len + seed_len + right.len());
        steps.extend(left);
        steps.extend(seed);
        steps.extend(right);

        Self {
            steps,
            seed_range: left_len..(left_len + seed_len),
        }
    }

    /// Returns the full alignment steps
    pub fn steps(&self) -> &[Pairing] {
        &self.steps
    }

    /// Returns only the seed region steps
    pub fn seed(&self) -> &[Pairing] {
        &self.steps[self.seed_range.clone()]
    }

    /// Returns the 5' extension (Left of seed)
    pub fn left_extension(&self) -> &[Pairing] {
        &self.steps[..self.seed_range.start]
    }

    /// Returns the 3' extension (Right of seed)
    pub fn right_extension(&self) -> &[Pairing] {
        &self.steps[self.seed_range.end..]
    }

    /// Generates the interaction string (e.g. "||| :::")
    pub fn fingerprint(&self) -> String {
        self.steps.iter().map(|p| p.to_char()).collect()
    }

    /// Generates the target sequence string (e.g. "accu--cg")
    pub fn target_sequence(&self) -> String {
        self.steps.iter().map(|p| p.target_char()).collect()
    }
}
/// Statistics for search filtering
#[derive(Debug, Default)]
pub struct SearchStats {
    pub seeds_tried: usize,
    pub candidates_processed: usize,
    pub hits_before_dedup: usize,
    pub hits_final: usize,
    pub filtered: HashMap<FilterReason, usize>,
}

impl SearchStats {
    pub fn record_filter(&mut self, reason: FilterReason) {
        *self.filtered.entry(reason).or_insert(0) += 1;
    }
}

pub trait Sequence {
    fn reverse_complement_dna(&self) -> Vec<u8>;
    fn reverse_complement_rna(&self) -> Vec<u8>;
}

impl Sequence for [u8] {
    fn reverse_complement_dna(&self) -> Vec<u8> {
        self.iter()
            .rev()
            .map(|&b| {
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
    pub strand: Strand,
}

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SearchHit {
    pub query_id: String,
    pub target_id: String,

    pub q_start: usize,
    pub q_end: usize,
    pub t_start: usize,
    pub t_end: usize,
    pub output_t_start: usize, // 1-based, strand-aware
    pub output_t_end: usize,   // 1-based, strand-aware
    pub strand: char,
    pub energy: f64,
    pub alignment: Alignment,
    pub flank_5: String,
    pub flank_3: String,
}

impl SearchHit {
    pub fn write(&self, w: &mut dyn Write) -> std::io::Result<()> {
        // Normalize strings for output (T->U)
        // Derive from alignment
        let norm_fp = self.alignment.fingerprint();
        let norm_ts = self
            .alignment
            .target_sequence()
            .replace('T', "U")
            .replace('t', "u");

        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}",
            self.query_id
                .split_whitespace()
                .next()
                .unwrap_or(&self.query_id),
            self.q_start + 1,
            self.q_end + 1,
            self.target_id
                .split_whitespace()
                .next()
                .unwrap_or(&self.target_id),
            self.output_t_start,
            self.output_t_end,
            self.strand,
            self.energy,
            norm_fp,
            norm_ts,
            self.flank_5,
            self.flank_3
        )
    }
}

// Reimplementing mapping locally for safety and speed

pub struct SaIndex<'a> {
    pub index: &'a IndexFile,
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

        debug!(
            "FIND_CAND: seed={} pairing={:?}",
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

            // 2. Search REVERSE COMPLEMENT
            let rc_seq = seq_idx.sequence.reverse_complement_dna();
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
                "FIND_CAND: seq_idx={} name={} fwd_hits={} rc_hits={}",
                i, &seq_idx.name, fwd_count, rc_count
            );
        }

        debug!("FIND_CAND: total_candidates={}", candidates.len());
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
        pairing: SeedPairing,
        seq_idx: usize,
        strand: Strand,
        candidates: &mut Vec<SeedCandidate>,
    ) {
        trace!(
            "SA_SEARCH: START seed={} pairing={:?} idx={} strand={:?}",
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
                    "SA_SEARCH: Wobble Check offset={} char={} pairing={:?} wobble_char={:?}",
                    offset, target_char as char, pairing, wobble_char
                );

                if let Some(wc) = wobble_char {
                    let (ws, we) = self.get_sa_interval(sa, text, start, end, offset, wc);
                    if ws < we {
                        trace!(
                            "SA_SEARCH: Wobble Found! offset={} range={}-{}",
                            offset, ws, we
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
}

// --- DP Structure ---

/// DP alignment state - replaces magic integers (0=M, 1=Bq, 2=Bt)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DpState {
    Match, // M - Match/Mismatch state (paired bases)
    GapQ,  // Bq - Query Bulge state (gap in target, query base unpaired)
    GapT,  // Bt - Target Bulge state (gap in query, target base unpaired)
}

#[derive(Clone, Debug)]
pub struct Grid<T> {
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

    /// Access the diagonal neighbor (i-1, j-1)
    #[inline(always)]
    pub fn diag(&self, i: usize, j: usize) -> T {
        self.get(i - 1, j - 1)
    }

    /// Access the upper neighbor (i-1, j) - Gap in Target
    #[inline(always)]
    pub fn up(&self, i: usize, j: usize) -> T {
        self.get(i - 1, j)
    }

    /// Access the left neighbor (i, j-1) - Gap in Query
    #[inline(always)]
    pub fn left(&self, i: usize, j: usize) -> T {
        self.get(i, j - 1)
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
pub enum DpMove {
    Stop,
    Match, // From M(i-1, j-1)
    GapQ,  // From Bq(i, j)
    GapT,  // From Bt(i, j)
}

impl Default for DpMove {
    fn default() -> Self {
        Self::Stop
    }
}

pub struct DpContext {
    // Score matrices (None = not reachable, Some(score) = reachable with score)
    m: ScoreGrid,
    bq: ScoreGrid,
    bt: ScoreGrid,

    // Traceback matrices
    tb_m: Grid<DpMove>,
    tb_bq: Grid<DpMove>,
    tb_bt: Grid<DpMove>,
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

pub struct SearchContext<'a> {
    pub index: &'a SaIndex<'a>,
    pub args: &'a SearchArgs,
    pub dp_ctx: DpContext,
    pub stats: SearchStats,
}

impl<'a> SearchContext<'a> {
    pub fn new(index: &'a SaIndex<'a>, args: &'a SearchArgs) -> Self {
        Self {
            index,
            args,
            dp_ctx: DpContext::new(200, 200),
            stats: SearchStats::default(),
        }
    }
}

pub fn run_search(
    queries: &[(String, Vec<u8>)],
    index: &SaIndex<'_>,
    output: impl AsRef<Path>,
    opts: &SearchArgs,
) -> Result<()> {
    info!(
        "Starting search: {} queries, seed={:?}, max_ext={}, delta_g={}",
        queries.len(),
        opts.seed.seed,
        opts.extend.max_extension,
        opts.extend.delta_g
    );

    // Create output writer
    let mut writer: Box<dyn Write> = if output.as_ref() == Path::new("-") {
        Box::new(std::io::stdout())
    } else {
        Box::new(std::fs::File::create(output.as_ref()).context("Failed to create output file")?)
    };

    debug!("SEARCH: output={:?}", output.as_ref());

    // Create Search Context
    let mut ctx = SearchContext::new(index, opts);

    let mut all_hits = Vec::new();

    for (q_id, q_seq) in queries {
        debug!("QUERY: id={} len={}", q_id, q_seq.len());
        trace!("QUERY_SEQ: {}", String::from_utf8_lossy(q_seq));

        // Find seeds
        let seeds = find_seeds_for_query(q_seq, &mut ctx)?;
        debug!("QUERY: {} candidates found", seeds.len());
        ctx.stats.candidates_processed += seeds.len();

        for candidate in &seeds {
            trace!(
                "CANDIDATE: q_pos={} t_idx={} t_start={} len={} strand={:?}",
                candidate.query_pos,
                candidate.target_idx,
                candidate.target_start,
                candidate.len,
                candidate.strand
            );
            if let Some(hit) = process_candidate(q_id, q_seq, candidate, &mut ctx) {
                trace!(
                    "HIT_ACCEPTED: q={}-{} t={}-{} E={:.2}",
                    hit.q_start, hit.q_end, hit.t_start, hit.t_end, hit.energy
                );
                all_hits.push(hit);
            }
        }
    }

    // Deduplicate logic
    ctx.stats.hits_before_dedup = all_hits.len();
    let deduped = deduplicate_hits(all_hits, &mut ctx.stats);
    ctx.stats.hits_final = deduped.len();

    // Log filter stats
    info!(
        "Search complete: {} hits ({} before dedup)",
        ctx.stats.hits_final, ctx.stats.hits_before_dedup
    );
    if !ctx.stats.filtered.is_empty() {
        let filter_summary: Vec<String> = ctx
            .stats
            .filtered
            .iter()
            .map(|(r, c)| format!("{:?}={}", r, c))
            .collect();
        info!("Filtered: {}", filter_summary.join(", "));
    }

    for hit in deduped {
        hit.write(&mut writer)?;
    }

    Ok(())
}

/// Check if hit `k` shadows hit `h` (k is better and contains h)
fn shadows(k: &SearchHit, h: &SearchHit) -> Option<FilterReason> {
    // Exact match - identical coordinates
    if k.q_start == h.q_start && k.q_end == h.q_end && k.t_start == h.t_start && k.t_end == h.t_end
    {
        return Some(FilterReason::DedupExactMatch);
    }

    // Check containment
    let q_contained = k.q_start <= h.q_start && k.q_end >= h.q_end;
    let t_contained = k.t_start <= h.t_start && k.t_end >= h.t_end;

    if !q_contained || !t_contained {
        return None;
    }

    // h must start strictly after k (not share same start)
    if h.q_start == k.q_start {
        return None;
    }

    Some(FilterReason::DedupContainedByShadow)
}

fn deduplicate_hits(mut hits: Vec<SearchHit>, stats: &mut SearchStats) -> Vec<SearchHit> {
    debug!("DEDUP: input_count={}", hits.len());

    if hits.is_empty() {
        return hits;
    }

    // Sort by Energy ascending (best first)
    hits.sort_by(|a, b| {
        a.energy
            .partial_cmp(&b.energy)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.q_start.cmp(&b.q_start))
    });

    let mut kept: Vec<SearchHit> = Vec::new();

    for h in hits {
        if let Some(reason) = kept.iter().find_map(|k| shadows(k, &h)) {
            stats.record_filter(reason);
            trace!(
                "DEDUP: FILTERED q{}-{}:t{}-{} reason={:?}",
                h.q_start, h.q_end, h.t_start, h.t_end, reason
            );
        } else {
            kept.push(h);
        }
    }

    debug!("DEDUP: output_count={}", kept.len());
    kept
}

//TODO: improve performance via better search strategies!

fn find_seeds_for_query(q_seq: &[u8], ctx: &mut SearchContext<'_>) -> Result<Vec<SeedCandidate>> {
    let seed_spec_str = ctx.args.seed.seed.as_deref().unwrap_or("17");
    let seed_len_specs = SeedSpec::from_str(seed_spec_str)
        .map_err(|e| anyhow!("Failed to parse seed specification: {}", e))?
        .normalize(q_seq.len())
        .map_err(|e| anyhow!("Invalid seed spec for query length: {}", e))?;

    let mut candidates = Vec::new();
    let q_len = q_seq.len();

    let (start, end, mi_len) = seed_len_specs;
    let start0 = start - 1;
    let end0 = end - 1;

    debug!(
        "SEEDS: spec={} q_len={} range=({},{}) mi_len={} pairing={:?}",
        seed_spec_str, q_len, start, end, mi_len, ctx.args.seed.pairing
    );

    if start0 + mi_len > q_len {
        warn!(
            "SEEDS: seed range too long for query: start0={} mi_len={} q_len={}",
            start0, mi_len, q_len
        );
        return Ok(candidates);
    }

    let last_start = end0.saturating_sub(mi_len - 1);
    let n_skipped = 0usize;
    let pairing = ctx.args.seed.pairing;

    for q_pos in start0..=last_start {
        // Max seed length from this position
        let max_seed_len = (end0 + 1).saturating_sub(q_pos).min(q_len - q_pos);

        for seed_len in mi_len..=max_seed_len {
            let seed_seq = &q_seq[q_pos..q_pos + seed_len];
            if seed_seq.contains(&b'N') || seed_seq.contains(&b'n') {
                ctx.stats.record_filter(FilterReason::SeedContainsN);
                trace!(
                    "SEEDS: FILTERED q_pos={} len={} reason=SeedContainsN",
                    q_pos, seed_len
                );
                continue;
            }

            // Index search (RC of seed)
            let seed_rc = seed_seq.reverse_complement_rna();
            // find_candidates returns Vec<SeedCandidate> now

            let mut hits = ctx.index.find_candidates(&seed_rc, pairing);

            trace!(
                "SEEDS: q_pos={} len={} seed={} rc={} hits={}",
                q_pos,
                seed_len,
                String::from_utf8_lossy(seed_seq),
                String::from_utf8_lossy(&seed_rc),
                hits.len()
            );

            for h in &mut hits {
                h.query_pos = q_pos;
                h.len = seed_len;
                candidates.push(SeedCandidate {
                    query_pos: q_pos,
                    target_idx: h.target_idx,
                    target_start: h.target_start,
                    len: seed_len,
                    strand: h.strand,
                });
            }
        }
    }

    debug!(
        "SEEDS: generated {} candidates ({} skipped for N)",
        candidates.len(),
        n_skipped
    );

    Ok(candidates)
}

fn process_candidate(
    q_id: &str,
    q_seq: &[u8],
    candidate: &SeedCandidate,
    ctx: &mut SearchContext<'_>,
) -> Option<SearchHit> {
    let t_idx = candidate.target_idx;
    let t_seq_cow = match candidate.strand {
        Strand::Reverse => std::borrow::Cow::Owned(ctx.index.get_sequence_rc(t_idx)),
        Strand::Forward => std::borrow::Cow::Borrowed(ctx.index.get_sequence(t_idx)),
    };
    let t_seq = &t_seq_cow;
    let t_start_idx = candidate.target_start;
    let seed_len = candidate.len;
    let q_pos = candidate.query_pos;

    trace!(
        "PROC_CAND: q_id={} t_idx={} q_pos={} t_start={} seed_len={} strand={:?}",
        q_id, t_idx, q_pos, t_start_idx, seed_len, candidate.strand
    );

    if t_start_idx + seed_len > t_seq.len() {
        ctx.stats.record_filter(FilterReason::SeedOutOfBounds);
        warn!(
            "PROC_CAND: FILTERED reason=SeedOutOfBounds t_start={} seed_len={} t_len={}",
            t_start_idx,
            seed_len,
            t_seq.len()
        );
        return None;
    }

    // Call extend_seed (or essentially reproduce its valuable logic).
    let extension_result = extend_seed(ctx, q_seq, t_seq, candidate)?;

    let ext = extension_result;
    let score = ext.score;

    if score > ctx.args.extend.delta_g {
        ctx.stats.record_filter(FilterReason::EnergyAboveThreshold);
        trace!(
            "PROC_CAND: FILTERED reason=EnergyAboveThreshold score={:.2} > delta_g={}",
            score, ctx.args.extend.delta_g
        );
        return None;
    }

    let l_q = ext.l_q;
    let l_t = ext.l_t;
    let r_q = ext.r_q;
    let r_t = ext.r_t;
    let seed_q = q_pos;
    let seed_t = t_start_idx;

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
        target_id: ctx.index.get_id(t_idx).to_string(),

        q_start: final_q_start,
        q_end: final_q_end,
        t_start: final_t_start,
        t_end: final_t_end,
        output_t_start: out_t_start,
        output_t_end: out_t_end,
        strand: strand_char,
        energy: score,
        alignment: ext.alignment,
        flank_5,
        flank_3: flank_3,
    })
}

pub struct DpResult {
    pub score: i32,
    pub ext_q_len: usize,
    pub ext_t_len: usize,
    pub trace: Vec<DpMove>,
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
    ctx: &mut SearchContext<'_>,
    q_seq: &[u8],
    t_seq: &[u8],
    candidate: &SeedCandidate,
) -> Option<ExtensionResult> {
    let q_pos = candidate.query_pos;
    let t_pos = candidate.target_start;
    let len = candidate.len;
    let opts = &ctx.args.extend;

    // MAXIMALITY CHECK
    // Skip non-maximal seeds: if the seed can be extended by a valid base pair
    // on either end, it's a sub-seed of a longer match and will

    trace!(
        "MAXIMALITY: ENTERING extend_seed q_pos={} t_pos={} len={} delta_g={}",
        q_pos, t_pos, len, opts.delta_g
    );

    // 1. Left extendable?
    if q_pos > 0 && t_pos + len < t_seq.len() {
        let q_prev = Base::from_byte(q_seq[q_pos - 1]).idx();
        let t_next = Base::from_byte(t_seq[t_pos + len]).idx();
        let p_class = PAIR_MAT[q_prev][t_next];

        trace!(
            "MAXIMALITY: Left Check q_pos={} t_pos={} len={} q_prev={} t_next={} pair={}",
            q_pos, t_pos, len, q_prev, t_next, p_class
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
                    "MAXIMALITY: FILTERED reason=MaximalityLeft q_pos={} t_pos={} len={} pair={}",
                    q_pos, t_pos, len, p_class
                );
                return None;
            }
        }
    }

    // 2. Right extendable? Check if q[q_pos+len] pairs with t[t_pos-1]
    if q_pos + len < q_seq.len() && t_pos > 0 {
        let q_next = Base::from_byte(q_seq[q_pos + len]).idx();
        let t_prev = Base::from_byte(t_seq[t_pos - 1]).idx();
        let p_class = PAIR_MAT[q_next][t_prev];

        trace!(
            "MAXIMALITY: Right Check q_pos={} t_pos={} len={} q_next={} t_prev={} pair={}",
            q_pos, t_pos, len, q_next, t_prev, p_class
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
                    "MAXIMALITY: FILTERED reason=MaximalityRight q_pos={} t_pos={} len={} pair={}",
                    q_pos, t_pos, len, p_class
                );
                return None;
            }
        }
    }

    trace!(
        "MAXIMALITY: Accepted Seed: q_pos={} t_pos={} len={} delta_g={}",
        q_pos, t_pos, len, opts.delta_g
    );

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

    let max_ext = opts.max_extension as usize;
    let safe_ext = max_ext.min(MAX_DP_EXT);

    // DP Left: Extend Query Left (5'), Target Right (3')
    // Start at: q_pos (5' Q), t_match_end (3' T)
    let left_res = dp_left(&mut ctx.dp_ctx, q_seq, t_seq, q_pos, t_match_end, safe_ext);

    // DP Right: Extend Query Right (3'), Target Left (5')
    // Start at: q_pos + len - 1 (3' Q), t_pos (5' T)
    let right_res = dp_right(
        &mut ctx.dp_ctx,
        q_seq,
        t_seq,
        q_pos + len - 1,
        t_pos,
        safe_ext,
    );

    let final_score =
        (seed_energy + left_res.score as f64 + right_res.score as f64 - 559.0) / -100.0;

    debug!(
        "EXTEND: score={:.2} (seed={:.2} L={} R={}) L_len={}/{} R_len={}/{}",
        final_score,
        seed_energy / -100.0,
        left_res.score,
        right_res.score,
        left_res.ext_q_len,
        left_res.ext_t_len,
        right_res.ext_q_len,
        right_res.ext_t_len
    );

    let mut left_alignment = Vec::new();

    // 1. Left Trace (Query 5' -> Seed)
    // l_trace is 5'->3' (from far left to seed start)
    // Coordinates: q_start - i, t_start + j
    // trace_vec[k] corresponds to step from (i,j) to (i-1, j) etc.
    // We need to replay properly.
    // Ideally we reconstruct by iterating trace and tracking (i, j).
    // Start at best_i, best_j.
    {
        let mut i = left_res.ext_q_len;
        let mut j = left_res.ext_t_len;
        // left_res.trace is ordered from [step at best_i] ... [step at 1].
        // So iterating it naturally goes from 5' end toward seed.

        for step in &left_res.trace {
            match step {
                DpMove::Match | DpMove::Stop => {
                    // Consumes both
                    let q_b = if i > 0 {
                        q_seq[candidate.query_pos - i]
                    } else {
                        b'N'
                    };
                    let t_b = if t_match_end + j < t_seq.len() {
                        t_seq[t_match_end + j]
                    } else {
                        b'N'
                    };
                    left_alignment.push(Pairing::from_bases(q_b, t_b));
                    if i > 0 {
                        i -= 1;
                    }
                    if j > 0 {
                        j -= 1;
                    }
                }
                DpMove::GapQ => {
                    // Gap in Query? No, GapQ means Bq matrix (Query has base, Target has Gap).
                    // In dp_left: trace "GapQ" means we came from Bq.
                    // Bq state: "Insertion in Query" (vs Target).
                    // So Query has Base. Target has Gap.
                    // So Pairing::GapTarget.

                    let q_b = if i > 0 {
                        q_seq[candidate.query_pos - i]
                    } else {
                        b'N'
                    };
                    left_alignment.push(Pairing::GapTarget(q_b));
                    if i > 0 {
                        i -= 1;
                    }
                }
                DpMove::GapT => {
                    // GapT means Bt matrix. Target has base. Query has Gap.
                    // So Pairing::GapQuery.

                    let t_b = if t_match_end + j < t_seq.len() {
                        t_seq[t_match_end + j]
                    } else {
                        b'N'
                    };
                    left_alignment.push(Pairing::GapQuery(t_b));
                    if j > 0 {
                        j -= 1;
                    }
                }
            }
        }
    }

    // 2. Seed itself
    let mut seed_alignment = Vec::new();
    for n in 0..len {
        let q_idx = candidate.query_pos + n;
        // candidate.target_start is start (lowest index).
        // Since seed is antiparallel: Q binds T.
        // Q: 5'->3' (idx +n).
        // T: 3'->5' (idx -n).
        // t_match_end is the 3' end of target site (matches 5' of query).
        // So t_match_end corresponds to q_pos.
        // t_match_end - n corresponds to q_pos + n.
        let t_idx = if t_match_end >= n { t_match_end - n } else { 0 };

        let q_b = q_seq[q_idx];
        let t_b = t_seq[t_idx];
        seed_alignment.push(Pairing::from_bases(q_b, t_b));
    }

    // 3. Right Trace (Query 3' -> end)
    let mut right_alignment = Vec::new();
    {
        let mut curr_i = 0;
        let mut curr_j = 0;

        // We iterate reversed trace (Start -> End)
        for step in right_res.trace.iter().rev() {
            match step {
                DpMove::Match | DpMove::Stop => {
                    curr_i += 1;
                    curr_j += 1;
                    let q_b = if candidate.query_pos + len - 1 + curr_i < q_seq.len() {
                        q_seq[candidate.query_pos + len - 1 + curr_i]
                    } else {
                        b'N'
                    };

                    // t_pos is 5' end of seed (matches 3' of query).
                    // Matches q_pos + len - 1.
                    // As we extend right (3' of query), we extend left (5' of target).
                    // So Target Index decreases.
                    // t_b = t_seq[t_pos - curr_j].
                    let t_b = if t_pos >= curr_j {
                        t_seq[t_pos - curr_j]
                    } else {
                        b'N'
                    };
                    right_alignment.push(Pairing::from_bases(q_b, t_b));
                }
                DpMove::GapQ => {
                    // DpMove::GapQ in dp_right.
                    // i increases (Query). j same.
                    // Query has Base. Target Gap.
                    // Pairing::GapTarget.
                    curr_i += 1;
                    let q_b = if candidate.query_pos + len - 1 + curr_i < q_seq.len() {
                        q_seq[candidate.query_pos + len - 1 + curr_i]
                    } else {
                        b'N'
                    };
                    right_alignment.push(Pairing::GapTarget(q_b));
                }
                DpMove::GapT => {
                    // DpMove::GapT in dp_right.
                    // j increases (Target 5'). i same.
                    // Target Base. Query Gap.
                    curr_j += 1;
                    let t_b = if t_pos >= curr_j {
                        t_seq[t_pos - curr_j]
                    } else {
                        b'N'
                    };
                    right_alignment.push(Pairing::GapQuery(t_b));
                }
            }
        }
    }

    let alignment = Alignment::new(left_alignment, seed_alignment, right_alignment);

    Some(ExtensionResult {
        score: final_score,
        alignment,

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

    trace!(
        "DP_LEFT: START q_len={} t_len={} initial_best_e={}",
        q_len, t_len, best_e
    );

    // C: if (lq <= 1 || lt <= 1) return best_e;
    if q_len <= 1 || t_len <= 1 {
        return DpResult {
            score: best_e,
            ext_q_len: best_i,
            ext_t_len: best_j,
            trace: Vec::new(),
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
        if let Some(prev_bt) = bt.left(0, j) {
            bt.set(
                0,
                j,
                Some(prev_bt + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32),
            );
            tb_bt.set(0, j, DpMove::GapT);

            if q_len >= 1 {
                let new_m = prev_bt + s_mat[q_char(1)][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32;
                m.set(1, j, Some(new_m));
                tb_m.set(1, j, DpMove::GapT);
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
        if let Some(prev_bq) = bq.up(i, 0) {
            bq.set(
                i,
                0,
                Some(prev_bq + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32),
            );
            tb_bq.set(i, 0, DpMove::GapQ);

            if t_len >= 1 {
                let new_m = prev_bq + s_mat[q_char(i)][q_char(i - 1)][t_comp(1)][GAP_IDX] as i32;
                m.set(i, 1, Some(new_m));
                tb_m.set(i, 1, DpMove::GapQ);
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
        if let Some(m11) = m.diag(2, 2) {
            bt.set(
                1,
                2,
                Some(m11 + s_mat[GAP_IDX][q_char(1)][t_comp(2)][t_comp(1)] as i32),
            );
            tb_bt.set(1, 2, DpMove::Match);

            bq.set(
                2,
                1,
                Some(m11 + s_mat[q_char(2)][q_char(1)][GAP_IDX][t_comp(1)] as i32),
            );
            tb_bq.set(2, 1, DpMove::Match);

            let m22 = m11 + s_mat[q_char(2)][q_char(1)][t_comp(2)][t_comp(1)] as i32;
            m.set(2, 2, Some(m22));
            tb_m.set(2, 2, DpMove::Match);

            let val = m22 + s_mat[GAP_IDX][q_char(2)][GAP_IDX][t_comp(2)] as i32;
            if val > best_e {
                best_e = val;
                best_i = 2;
                best_j = 2;
            }
        }
        if let Some(m12) = m.up(2, 2) {
            bq.set(
                2,
                2,
                Some(m12 + s_mat[q_char(2)][q_char(1)][GAP_IDX][t_comp(2)] as i32),
            );
            tb_bq.set(2, 2, DpMove::Match);
        }
        if let Some(m21) = m.left(2, 2) {
            bt.set(
                2,
                2,
                Some(m21 + s_mat[GAP_IDX][q_char(2)][t_comp(2)][t_comp(1)] as i32),
            );
            tb_bt.set(2, 2, DpMove::Match);
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
                .diag(i, j)
                .map(|v| v + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][t_comp(j - 1)] as i32);
            let s_mq = bq
                .diag(i, j)
                .map(|v| v + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][GAP_IDX] as i32);
            let s_mt = bt
                .diag(i, j)
                .map(|v| v + s_mat[q_char(i)][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32);

            // Find best value and corresponding traceback step
            let (val_m, step_m) = [
                (s_mm, DpMove::Match),
                (s_mq, DpMove::GapQ),
                (s_mt, DpMove::GapT),
            ]
            .into_iter()
            .filter_map(|(opt, step)| opt.map(|v| (v, step)))
            .max_by_key(|(v, _)| *v)
            .unwrap_or((i32::MIN, DpMove::Stop));

            let val_m = if val_m == i32::MIN { None } else { Some(val_m) };
            m.set(i, j, val_m);
            tb_m.set(i, j, step_m);

            trace!(
                "DP_LEFT_CELL: i={} j={} M={} (mm={} mq={} mt={}) step={:?}",
                i,
                j,
                val_m.map_or("-".to_string(), |v| v.to_string()),
                s_mm.map_or("-".to_string(), |v| v.to_string()),
                s_mq.map_or("-".to_string(), |v| v.to_string()),
                s_mt.map_or("-".to_string(), |v| v.to_string()),
                step_m
            );

            if let Some(v) = val_m {
                let curr_e = v + s_mat[GAP_IDX][q_char(i)][GAP_IDX][t_comp(j)] as i32;
                if curr_e > best_e {
                    trace!(
                        "DP_LEFT: UPDATE best: i={} j={} curr_e={} (was {})",
                        i, j, curr_e, best_e
                    );
                    best_e = curr_e;
                    best_i = i;
                    best_j = j;
                }
            }

            // Calc Bq[i,j]
            if i > 2 || (i == 2 && j > 2) {
                let s_qm = m
                    .up(i, j)
                    .map(|v| v + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][t_comp(j)] as i32);
                let s_qq = bq
                    .up(i, j)
                    .map(|v| v + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32);

                // Priority to Match (opening) if tie
                match (s_qq, s_qm) {
                    (Some(qq), Some(qm)) if qq > qm => {
                        bq.set(i, j, Some(qq));
                        tb_bq.set(i, j, DpMove::GapQ);
                        trace!("DP_LEFT_CELL: i={} j={} Bq={} (from GapQ)", i, j, qq);
                    }
                    (_, Some(qm)) => {
                        bq.set(i, j, Some(qm));
                        tb_bq.set(i, j, DpMove::Match);
                        trace!("DP_LEFT_CELL: i={} j={} Bq={} (from Match)", i, j, qm);
                    }
                    (Some(qq), None) => {
                        bq.set(i, j, Some(qq));
                        tb_bq.set(i, j, DpMove::GapQ);
                        trace!("DP_LEFT_CELL: i={} j={} Bq={} (from GapQ, no M)", i, j, qq);
                    }
                    _ => {}
                }
            }

            // Calc Bt[i,j]
            if j > 2 || (j == 2 && i > 2) {
                let s_tm = m
                    .left(i, j)
                    .map(|v| v + s_mat[GAP_IDX][q_char(i)][t_comp(j)][t_comp(j - 1)] as i32);
                let s_tt = bt
                    .left(i, j)
                    .map(|v| v + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32);

                // Priority to Match (opening) if tie
                match (s_tt, s_tm) {
                    (Some(tt), Some(tm)) if tt > tm => {
                        bt.set(i, j, Some(tt));
                        tb_bt.set(i, j, DpMove::GapT);
                        trace!("DP_LEFT_CELL: i={} j={} Bt={} (from GapT)", i, j, tt);
                    }
                    (_, Some(tm)) => {
                        bt.set(i, j, Some(tm));
                        tb_bt.set(i, j, DpMove::Match);
                        trace!("DP_LEFT_CELL: i={} j={} Bt={} (from Match)", i, j, tm);
                    }
                    (Some(tt), None) => {
                        bt.set(i, j, Some(tt));
                        tb_bt.set(i, j, DpMove::GapT);
                        trace!("DP_LEFT_CELL: i={} j={} Bt={} (from GapT, no M)", i, j, tt);
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
    let mut trace_vec = Vec::new();
    let mut state = DpState::Match;

    while i > 0 || j > 0 {
        match state {
            DpState::Match => {
                if i == 0 || j == 0 {
                    // End of traceback
                    break;
                }

                let step = tb_m.get(i, j);
                i -= 1;
                j -= 1;

                trace_vec.push(step); // Push DpMove

                match step {
                    DpMove::Stop => break,
                    DpMove::Match => state = DpState::Match,
                    DpMove::GapQ => state = DpState::GapQ,
                    DpMove::GapT => state = DpState::GapT,
                }
            }
            DpState::GapQ => {
                let step = tb_bq.get(i, j);
                trace_vec.push(step);

                if i > 0 {
                    i -= 1;
                } else {
                    break;
                }

                match step {
                    DpMove::Stop => break,
                    DpMove::Match => state = DpState::Match,
                    DpMove::GapQ => state = DpState::GapQ,
                    DpMove::GapT => state = DpState::GapT, // Should not occur
                }
            }
            DpState::GapT => {
                let step = tb_bt.get(i, j);
                trace_vec.push(step);

                if j > 0 {
                    j -= 1;
                } else {
                    break;
                }

                match step {
                    DpMove::Stop => break,
                    DpMove::Match => state = DpState::Match,
                    DpMove::GapT => state = DpState::GapT,
                    DpMove::GapQ => state = DpState::GapQ, // Should not occur
                }
            }
        }
    }

    trace!(
        "DP_LEFT: END score={:.2} q_ext={} t_ext={} trace_len={}",
        best_e,
        best_i,
        best_j,
        trace_vec.len()
    );

    DpResult {
        score: best_e,
        ext_q_len: best_i,
        ext_t_len: best_j,
        trace: trace_vec,
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
            trace: Vec::new(),
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
        if let Some(prev_bt) = bt.left(0, j) {
            bt.set(
                0,
                j,
                Some(prev_bt + s_mat[GAP_IDX][GAP_IDX][t_comp(j - 1)][t_comp(j)] as i32),
            );
            tb_bt.set(0, j, DpMove::GapT);

            if q_len >= 1 {
                let new_m = prev_bt + s_mat[GAP_IDX][q_char(1)][t_comp(j - 1)][t_comp(j)] as i32;
                m.set(1, j, Some(new_m));
                tb_m.set(1, j, DpMove::GapT);
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
        if let Some(prev_bq) = bq.up(i, 0) {
            bq.set(
                i,
                0,
                Some(prev_bq + s_mat[q_char(i - 1)][q_char(i)][GAP_IDX][GAP_IDX] as i32),
            );
            tb_bq.set(i, 0, DpMove::GapQ);

            if t_len >= 1 {
                let new_m = prev_bq + s_mat[q_char(i - 1)][q_char(i)][GAP_IDX][t_comp(1)] as i32;
                m.set(i, 1, Some(new_m));
                tb_m.set(i, 1, DpMove::GapQ);
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
        if let Some(m11) = m.diag(2, 2) {
            bt.set(
                1,
                2,
                Some(m11 + s_mat[q_char(1)][GAP_IDX][t_comp(1)][t_comp(2)] as i32),
            );
            tb_bt.set(1, 2, DpMove::Match);

            bq.set(
                2,
                1,
                Some(m11 + s_mat[q_char(1)][q_char(2)][t_comp(1)][GAP_IDX] as i32),
            );
            tb_bq.set(2, 1, DpMove::Match);

            let m22 = m11 + s_mat[q_char(1)][q_char(2)][t_comp(1)][t_comp(2)] as i32;
            m.set(2, 2, Some(m22));
            tb_m.set(2, 2, DpMove::Match);

            let val = m22 + s_mat[q_char(2)][GAP_IDX][t_comp(2)][GAP_IDX] as i32;
            if val > best_e {
                best_e = val;
                best_i = 2;
                best_j = 2;
            }
        }
        if let Some(m12) = m.up(2, 2) {
            bq.set(
                2,
                2,
                Some(m12 + s_mat[q_char(1)][q_char(2)][t_comp(2)][GAP_IDX] as i32),
            );
            tb_bq.set(2, 2, DpMove::Match);
        }
        if let Some(m21) = m.left(2, 2) {
            bt.set(
                2,
                2,
                Some(m21 + s_mat[q_char(2)][GAP_IDX][t_comp(1)][t_comp(2)] as i32),
            );
            tb_bt.set(2, 2, DpMove::Match);
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
                .diag(i, j)
                .map(|v| v + s_mat[q_char(i - 1)][q_char(i)][t_comp(j - 1)][t_comp(j)] as i32);
            let s_mq = bq
                .diag(i, j)
                .map(|v| v + s_mat[q_char(i - 1)][q_char(i)][GAP_IDX][t_comp(j)] as i32);
            let s_mt = bt
                .diag(i, j)
                .map(|v| v + s_mat[GAP_IDX][q_char(i)][t_comp(j - 1)][t_comp(j)] as i32);

            // Find best value and corresponding traceback step
            let (val_m, step_m) = [
                (s_mm, DpMove::Match),
                (s_mq, DpMove::GapQ),
                (s_mt, DpMove::GapT),
            ]
            .into_iter()
            .filter_map(|(opt, step)| opt.map(|v| (v, step)))
            .max_by_key(|(v, _)| *v)
            .unwrap_or((i32::MIN, DpMove::Stop));

            let val_m = if val_m == i32::MIN { None } else { Some(val_m) };
            m.set(i, j, val_m);
            tb_m.set(i, j, step_m);

            trace!(
                "DP_RIGHT_CELL: i={} j={} M={} (mm={} mq={} mt={}) step={:?}",
                i,
                j,
                val_m.map_or("-".to_string(), |v| v.to_string()),
                s_mm.map_or("-".to_string(), |v| v.to_string()),
                s_mq.map_or("-".to_string(), |v| v.to_string()),
                s_mt.map_or("-".to_string(), |v| v.to_string()),
                step_m
            );

            if let Some(v) = val_m {
                let term = s_mat[q_char(i)][GAP_IDX][t_comp(j)][GAP_IDX] as i32;
                if v + term > best_e {
                    trace!(
                        "DP_RIGHT: UPDATE best: i={} j={} curr_e={} (was {})",
                        i,
                        j,
                        v + term,
                        best_e
                    );
                    best_e = v + term;
                    best_i = i;
                    best_j = j;
                }
            }

            // Calc Bq[i,j]
            if i > 2 || (i == 2 && j > 2) {
                let s_qm = m
                    .up(i, j)
                    .map(|v| v + s_mat[q_char(i - 1)][q_char(i)][t_comp(j)][GAP_IDX] as i32);
                let s_qq = bq
                    .up(i, j)
                    .map(|v| v + s_mat[q_char(i - 1)][q_char(i)][GAP_IDX][GAP_IDX] as i32);

                match (s_qq, s_qm) {
                    (Some(qq), Some(qm)) if qq > qm => {
                        bq.set(i, j, Some(qq));
                        tb_bq.set(i, j, DpMove::GapQ);
                    }
                    (_, Some(qm)) => {
                        bq.set(i, j, Some(qm));
                        tb_bq.set(i, j, DpMove::Match);
                    }
                    (Some(qq), None) => {
                        bq.set(i, j, Some(qq));
                        tb_bq.set(i, j, DpMove::GapQ);
                    }
                    _ => {}
                }
            }

            // Calc Bt[i,j]
            if j > 2 || (j == 2 && i > 2) {
                let s_tm = m
                    .left(i, j)
                    .map(|v| v + s_mat[q_char(i)][GAP_IDX][t_comp(j - 1)][t_comp(j)] as i32);
                let s_tt = bt
                    .left(i, j)
                    .map(|v| v + s_mat[GAP_IDX][GAP_IDX][t_comp(j - 1)][t_comp(j)] as i32);

                match (s_tt, s_tm) {
                    (Some(tt), Some(tm)) if tt > tm => {
                        bt.set(i, j, Some(tt));
                        tb_bt.set(i, j, DpMove::GapT);
                    }
                    (_, Some(tm)) => {
                        bt.set(i, j, Some(tm));
                        tb_bt.set(i, j, DpMove::Match);
                    }
                    (Some(tt), None) => {
                        bt.set(i, j, Some(tt));
                        tb_bt.set(i, j, DpMove::GapT);
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
    let mut trace_vec = Vec::new(); // Changed to Vec
    let mut state = DpState::Match;

    while i > 0 || j > 0 {
        match state {
            DpState::Match => {
                // Current state is Match (M[i,j])
                if i == 0 || j == 0 {
                    break;
                }

                let step = tb_m.get(i, j);
                i -= 1;
                j -= 1;

                trace_vec.push(step);

                match step {
                    DpMove::Stop => break,
                    DpMove::Match => state = DpState::Match,
                    DpMove::GapQ => state = DpState::GapQ,
                    DpMove::GapT => state = DpState::GapT,
                }
            }
            DpState::GapQ => {
                // Current state is GapQ (Bq[i,j])
                let step = tb_bq.get(i, j);
                trace_vec.push(step);

                if i > 0 {
                    i -= 1;
                } else {
                    break;
                }
                // j stays same

                match step {
                    DpMove::Stop => break,
                    DpMove::Match => state = DpState::Match,
                    DpMove::GapQ => state = DpState::GapQ,
                    DpMove::GapT => state = DpState::GapT, // Should not occur
                }
            }
            DpState::GapT => {
                // Current state is GapT (Bt[i,j])
                let step = tb_bt.get(i, j);
                trace_vec.push(step);

                if j > 0 {
                    j -= 1;
                } else {
                    break;
                }
                // i stays same

                match step {
                    DpMove::Stop => break,
                    DpMove::Match => state = DpState::Match,
                    DpMove::GapT => state = DpState::GapT,
                    DpMove::GapQ => state = DpState::GapQ, // Should not occur
                }
            }
        }
    }

    DpResult {
        score: best_e,
        ext_q_len: best_i,
        ext_t_len: best_j,
        trace: trace_vec,
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
