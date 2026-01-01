use anyhow::{Context, Result, anyhow};
use log::{debug, trace};
use std::io::Write;
use std::path::Path;
use std::str::FromStr;

use crate::SearchArgs;
use crate::dsm::{DSM_T04_POS, PAIR_MAT};
use crate::sa::IndexFile;
use crate::seed::SeedSpec;

const NA_VAL: i32 = i32::MIN / 2;
/// Large offset to ensure NA_VAL comparisons are deterministically false during backtracking
const NA_FALLBACK: i32 = NA_VAL - 999999;
const MAX_DP_EXT: usize = 30;
const GAP_IDX: usize = 0;

// Reimplementing mapping locally for safety and speed
const NUCL_MAP: [usize; 256] = {
    let mut table = [5; 256];
    table[b'A' as usize] = 1;
    table[b'a' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'C' as usize] = 3;
    table[b'c' as usize] = 3;
    table[b'U' as usize] = 4;
    table[b'u' as usize] = 4;
    table[b'T' as usize] = 4;
    table[b't' as usize] = 4; // T handled as U
    table[b'-' as usize] = 0; // gap
    table
};

/// Maps a nucleotide byte to its DSM index
/// - 0 = gap, 1 = A, 2 = G, 3 = C, 4 = U/T, 5 = other
#[inline(always)]
fn dsm_idx(b: u8) -> usize {
    NUCL_MAP[b as usize]
}

/// Trait to abstract over different index types (Suffix Array)
pub trait RisearchIndexTrait {
    /// Returns (seq_idx, position, is_reverse_complement)
    fn find_candidates(&self, seed: &[u8], wobble: bool) -> Vec<(usize, usize, bool)>;
    fn get_sequence(&self, seq_idx: usize) -> &[u8];
    fn get_sequence_rc(&self, seq_idx: usize) -> Vec<u8>; // RC of sequence
    fn get_id(&self, seq_idx: usize) -> &str;
    fn get_sequence_len(&self, seq_idx: usize) -> usize;
}

pub struct SaIndex<'a> {
    pub index: &'a IndexFile,
}

impl<'a> RisearchIndexTrait for SaIndex<'a> {
    fn find_candidates(&self, seed: &[u8], wobble: bool) -> Vec<(usize, usize, bool)> {
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
            let rc_seq = reverse_complement_dna(&seq_idx.sequence);
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

    fn get_sequence(&self, seq_idx: usize) -> &[u8] {
        &self.index.sequences[seq_idx].sequence
    }

    fn get_sequence_rc(&self, seq_idx: usize) -> Vec<u8> {
        reverse_complement_dna(&self.index.sequences[seq_idx].sequence)
    }

    fn get_id(&self, seq_idx: usize) -> &str {
        &self.index.sequences[seq_idx].name
    }

    fn get_sequence_len(&self, seq_idx: usize) -> usize {
        self.index.sequences[seq_idx].sequence.len()
    }
}

/// DNA reverse complement (for target sequences stored as DNA with T not U)
#[inline]
fn complement_dna(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b't',
        b'T' | b't' => b'a',
        b'C' | b'c' => b'g',
        b'G' | b'g' => b'c',
        _ => b'n',
    }
}

fn reverse_complement_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_dna(b)).collect()
}

impl<'a> SaIndex<'a> {
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
    Match, // M - Match/Mismatch
    GapQ,  // Bq - Bulge in Query
    GapT,  // Bt - Bulge in Target
}

struct Grid {
    data: Vec<i32>,
    width: usize,
    #[allow(dead_code)]
    height: usize,
}

impl Grid {
    fn new(width: usize, height: usize, default: i32) -> Self {
        Self {
            data: vec![default; width * height],
            width,
            height,
        }
    }

    #[inline(always)]
    fn get(&self, r: usize, c: usize) -> i32 {
        self.data[r * self.width + c]
    }

    #[inline(always)]
    fn set(&mut self, r: usize, c: usize, val: i32) {
        self.data[r * self.width + c] = val;
    }
}

pub fn run_search(
    queries: &[(String, Vec<u8>)],
    index: &impl RisearchIndexTrait,
    output: impl AsRef<Path>,
    opts: &SearchArgs,
) -> Result<()> {
    // Create output writer
    let mut writer: Box<dyn Write> = if output.as_ref() == Path::new("-") {
        Box::new(std::io::stdout())
    } else {
        Box::new(std::fs::File::create(output.as_ref()).context("Failed to create output file")?)
    };
    let seed_spec = opts.seed.as_deref().unwrap_or("6");
    let parsed_spec = SeedSpec::from_str(seed_spec)
        .map_err(|e| anyhow!("Failed to parse seed specification: {}", e))?;

    debug!("run_search called with {} queries", queries.len());

    for (q_id, q_seq) in queries {
        debug!("Processing query: {} (len={})", q_id, q_seq.len());
        if let Ok((start, end, len)) = parsed_spec.normalize(q_seq.len()) {
            trace!("Seed spec: start={}, end={}, len={}", start, end, len);
            let s0 = start - 1;
            let e0 = end - 1;
            if s0 + len > q_seq.len() {
                trace!(
                    "Skipping: s0+len ({}) > q_seq.len() ({})",
                    s0 + len,
                    q_seq.len()
                );
                continue;
            }
            let last_start = e0.saturating_sub(len - 1);
            trace!(
                "Seed positions: s0={}, e0={}, last_start={}",
                s0, e0, last_start
            );

            // Search for variable-length seeds at each position.
            // For each q_pos, we search for all seed lengths from min_len up to the
            // maximum possible length that fits within the query. Only maximal seeds
            // (those that cannot be extended further) will pass the maximality check.
            // This matches the C implementation's behavior where parallel SA traversal
            // naturally finds maximal seeds of varying lengths.
            let min_len = len;

            for q_pos in s0..=last_start {
                // Max seed length from this position that fits in query and seed region
                let max_seed_len = (e0 + 1).saturating_sub(q_pos).min(q_seq.len() - q_pos);

                // Search for seeds of all lengths from min_len to max_seed_len
                for seed_len in min_len..=max_seed_len {
                    let seed_seq = &q_seq[q_pos..q_pos + seed_len];
                    if seed_seq.contains(&b'N') || seed_seq.contains(&b'n') {
                        continue;
                    }

                    // Index search for Reverse Complement of Seed
                    let seed_rc = reverse_complement_rna(seed_seq);
                    // Pass wobble option
                    let hits = index.find_candidates(&seed_rc, opts.wobble);

                    if q_pos == s0 && seed_len == min_len {
                        trace!(
                            "First seed at q_pos={}: seed={:?}, rc={:?}, hits={}",
                            q_pos,
                            String::from_utf8_lossy(seed_seq),
                            String::from_utf8_lossy(&seed_rc),
                            hits.len()
                        );
                    }

                    for (t_idx, t_pos, is_antisense) in hits {
                        // is_antisense = true means hit found on RC of target (antisense strand)
                        // is_antisense = false means hit found on forward target (sense strand)

                        // Get the sequence we're aligning against
                        // Use Cow to avoid memory leak - owned Vec for RC, borrowed slice for forward
                        let t_seq_cow: std::borrow::Cow<'_, [u8]> = if is_antisense {
                            std::borrow::Cow::Owned(index.get_sequence_rc(t_idx))
                        } else {
                            std::borrow::Cow::Borrowed(index.get_sequence(t_idx))
                        };
                        let t_seq: &[u8] = &t_seq_cow;

                        if t_pos + seed_len > t_seq.len() {
                            continue;
                        }

                        // Maximality Check and Extension
                        // Returns None if the seed is skipped due to maximality (it can be extended left/right)
                        let Some((score, _)) =
                            extend_seed(q_seq, t_seq, q_pos, t_pos, seed_len, opts)
                        else {
                            continue;
                        };

                        // Debug: print hits with good scores
                        static BEST_SCORE: std::sync::atomic::AtomicI64 =
                            std::sync::atomic::AtomicI64::new(0);
                        let score_i64 = (score * 100.0) as i64;
                        if score_i64 < BEST_SCORE.load(std::sync::atomic::Ordering::Relaxed) {
                            BEST_SCORE.store(score_i64, std::sync::atomic::Ordering::Relaxed);
                            trace!(
                                "New best score: {:.2} at q_pos={}, t_idx={}, t_pos={}, is_antisense={}",
                                score, q_pos, t_idx, t_pos, is_antisense
                            );
                        }

                        // Debug: print first few scores
                        static DEBUG_COUNT: std::sync::atomic::AtomicUsize =
                            std::sync::atomic::AtomicUsize::new(0);
                        let count = DEBUG_COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if count < 10 {
                            trace!(
                                "Hit {}: q_pos={}, t_idx={}, t_pos={}, is_antisense={}, score={:.2}, threshold={:.2}",
                                count, q_pos, t_idx, t_pos, is_antisense, score, opts.delta_g
                            );
                        }

                        if score <= opts.delta_g {
                            let (q_start, q_end, t_start, t_end, l_q, l_t, r_q, r_t) =
                                find_best_extent(q_seq, t_seq, q_pos, t_pos, seed_len, opts);

                            // C convention for strand:
                            // '-' = query binds to SENSE (forward) strand of target mRNA
                            // '+' = query binds to ANTISENSE (reverse complement) strand
                            //
                            // For antisense hits (is_antisense=true):
                            //   We searched the RC sequence, positions are in RC coordinates
                            //   Convert back to forward strand: pos_fwd = len - pos_rc - 1
                            //   Since t_start/t_end are 0-based exclusive end: [t_start, t_end)
                            //   fwd_start = len - t_end, fwd_end = len - t_start (both 0-based)
                            //   Then +1 for 1-based output

                            let (output_t_start, output_t_end, strand) = if is_antisense {
                                let seq_len = index.get_sequence_len(t_idx);
                                let fwd_start = seq_len - t_end; // 1-based (Start = L - EndIdx)
                                let fwd_end = seq_len - t_start; // 1-based (End = L - StartIdx)
                                (fwd_start, fwd_end, '-')
                            } else {
                                (t_start + 1, t_end + 1, '+') // 0-based inclusive -> 1-based inclusive
                            };

                            print_detailed_output(
                                &mut writer,
                                q_id,
                                q_seq,
                                q_start,
                                q_end,
                                index.get_id(t_idx),
                                t_seq,
                                t_start,
                                t_end, // Use 0-based for slicing
                                output_t_start,
                                output_t_end, // Use 1-based for output
                                score,
                                strand,
                                l_q,
                                l_t,
                                r_q,
                                r_t,
                            )?;
                        }
                    }
                } // end for seed_len
            }
        }
    }
    Ok(())
}

fn extend_seed(
    q_seq: &[u8],
    t_seq: &[u8],
    q_pos: usize,
    t_pos: usize,
    len: usize,
    opts: &SearchArgs,
) -> Option<(f64, String)> {
    // MAXIMALITY CHECK
    // Check if the seed can be extended to the LEFT or RIGHT.
    // If it can, it is considered a sub-seed of a longer maximal seed, so we skip it.
    // This matches the logic in RIsearch2 (legacy C) `search.c`.

    // 1. Is Left Extendable?
    // Query extends Left (index -1). Target extends Right (index +len in T).
    // T is RC(Target), so increasing index moves 5'->3' on RC (which corresponds to 3'<-5' on Target).
    // Antiparallel binding: Q(Left) pairs with T(Right).
    // Check pair q[q_pos-1] vs t[t_pos+len].
    if q_pos > 0 && t_pos + len < t_seq.len() {
        let q_prev = dsm_idx(q_seq[q_pos - 1]);
        let t_next = dsm_idx(t_seq[t_pos + len]);
        if PAIR_MAT[q_prev][t_next] != 0 {
            // Extendable Left -> Skip
            return None;
        }
    }

    // 2. Is Right Extendable?
    // Query extends Right (index +len). Target extends Left (index -1).
    // Check pair q[q_pos+len] vs t[t_pos-1].
    // NOTE: C code logic skips if EITHER is extendable.
    if q_pos + len < q_seq.len() && t_pos > 0 {
        let q_next = dsm_idx(q_seq[q_pos + len]);
        let t_prev = dsm_idx(t_seq[t_pos - 1]);
        if PAIR_MAT[q_next][t_prev] != 0 {
            // Extendable Right -> Skip
            return None;
        }
    }

    let mut seed_energy = 0.0;

    // Seed energy - Antiparallel RNA Binding
    // Query goes 5'->3'. Target (RC) goes 5'->3'.
    // Binding is antiparallel: 5' Q pairs with 3' T.
    //
    // Q: q_pos (5') ... q_pos + len - 1 (3')
    // T: t_pos (5') ... t_pos + len - 1 (3')
    //
    // Pairing:
    // q[q_pos] (5') pairs with t[t_pos + len - 1] (3')
    // q[q_pos + k] pairs with t[t_pos + len - 1 - k]

    let t_match_end = t_pos + len - 1;
    let mut seed_int_str = String::with_capacity(len);

    for k in 0..len {
        let q_idx = q_pos + k;
        let t_idx = t_match_end - k;

        if k < len.saturating_sub(1) {
            let q_b1 = dsm_idx(q_seq[q_idx]);
            let q_b2 = dsm_idx(q_seq[q_idx + 1]);
            let t_b1 = dsm_idx(t_seq[t_idx]);
            let t_b2 = dsm_idx(t_seq[t_match_end - (k + 1)]);
            seed_energy += DSM_T04_POS[q_b1][q_b2][t_b1][t_b2] as f64;
        }

        // Build interaction string
        let qc = q_seq[q_idx];
        let tc = t_seq[t_idx];
        seed_int_str.push(get_fingerprint_char(qc, tc));
    }
    // Seed interaction string should match C output; omit any extra markers

    let max_ext = opts.max_extension as usize;
    let safe_ext = max_ext.min(MAX_DP_EXT);

    // DP Left: Extend Query Left (5'), Target Right (3')
    // Start at: q_pos (5' Q), t_match_end (3' T)
    let (l_score, _, _, l_trace) = dp_left(q_seq, t_seq, q_pos, t_match_end, safe_ext);

    // DP Right: Extend Query Right (3'), Target Left (5')
    // Start at: q_pos + len - 1 (3' Q), t_pos (5' T)
    let (r_score, _, _, r_trace) = dp_right(q_seq, t_seq, q_pos + len - 1, t_pos, safe_ext);

    let final_score = (seed_energy + l_score as f64 + r_score as f64 - 559.0) / -100.0;

    // Debug for first few calls in a clean way
    static DEBUG_EXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let c = DEBUG_EXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    if c < 5 {
        trace!(
            "extend_seed: seed_energy={:.0}, l_score={}, r_score={}, raw_total={:.0}, final={:.2}",
            seed_energy,
            l_score,
            r_score,
            seed_energy + l_score as f64 + r_score as f64,
            final_score
        );
    }

    // Combine output string here properly
    let l_str: String = l_trace.chars().rev().collect();
    let full_interaction = format!(
        "{}{}{}",
        l_str,
        seed_int_str,
        r_trace.chars().rev().collect::<String>()
    );

    Some((final_score, full_interaction))
}

fn dp_left(
    q_seq: &[u8],
    t_seq: &[u8],
    q_start: usize,
    t_start: usize,
    max_ext: usize,
) -> (i32, usize, usize, String) {
    // C code: rem_query_len_left = min(max_ext_len, q_start + 1)
    // The +1 includes the seed boundary position in the DP
    // q_len/t_len represent how many positions we can access (indices 0..q_len-1)
    // DP Left: Extend Query Left (towards 5'), Target Right (towards 3')
    // q_start: 5' end of seed (index decr)
    // t_start: 3' end of seed (index incr)

    // C: rem_query_len_left = min(max_ext_len, q_start + 1)
    let q_len = (q_start + 1).min(max_ext);
    // C: rem_target_len_left targets: t_start - sum_l[idx] + 1 ?
    // Target moves towards end of sequence.
    let t_len = (t_seq.len() - t_start - 1).min(max_ext);

    // ...

    let s_mat = &DSM_T04_POS;

    // Q uses q_start - i (Moving 3' -> 5' relative to q_seq 5'-3' storage?)
    // q_start is index. q_start - i.
    // If q_start is 5' end (lowest index).
    // Wait. If q_seq is 5'->3'.
    // q_start (lowest index).
    // q_start - i will be < 0 immediately if i >= 1.
    // So q_start must be HIGH index for q_start - i to work?
    //
    // Re-evaluating q_start.
    // `extend_seed` used `q_pos` (low index) to `q_pos + len - 1` (high index).
    // `dp_left`: extend query 5'.
    // We want to access `q_pos - 1`, `q_pos - 2`.
    // So `q_start` should be `q_pos`.
    // And `i` goes 1..len.
    // `q_start - i`.
    // q_len = (q_start + 1).min(max).
    // If q_start = 5. q_len = 6.
    // i=0 -> q[5]. i=1 -> q[4].
    //
    // YES. q_start - i.
    // So q_start MUST be passed as `q_pos` (index matching first char of seed).
    //
    // Target:
    // We want to extend Target Right (towards 3').
    // Access `t_3_end + 1`, `t_3_end + 2`.
    // So `t_start + j`.
    // So `t_start` MUST be passed as `t_pos + len - 1` (index matching last char of seed).

    let q_char = |i: usize| {
        if i > q_start {
            0
        } else {
            dsm_idx(q_seq[q_start - i])
        }
    };
    let t_comp = |j: usize| {
        if t_start + j >= t_seq.len() {
            0
        } else {
            dsm_idx(t_seq[t_start + j])
        }
    };

    // Initial score
    // C: best_e = (*S)[GAP][Q (0)][GAP][T (0)]
    // Q(0) = q[q_start].
    // T(0) = t[t_start].
    //
    // Note: T(0) is included.

    let mut best_e = s_mat[GAP_IDX][q_char(0)][GAP_IDX][t_comp(0)] as i32;
    let mut best_i = 0;
    let mut best_j = 0;

    // C: if (lq <= 1 || lt <= 1) return best_e;
    if q_len <= 1 || t_len <= 1 {
        return (best_e, best_i, best_j, String::new());
    }

    let mut m = Grid::new(t_len + 1, q_len + 1, NA_VAL);
    let mut bq = Grid::new(t_len + 1, q_len + 1, NA_VAL);
    let mut bt = Grid::new(t_len + 1, q_len + 1, NA_VAL);

    m.set(0, 0, 0);

    // Init (0,1), (1,0), (1,1)
    bt.set(0, 1, s_mat[GAP_IDX][q_char(0)][t_comp(1)][t_comp(0)] as i32);
    bq.set(1, 0, s_mat[q_char(1)][q_char(0)][GAP_IDX][t_comp(0)] as i32);
    m.set(
        1,
        1,
        s_mat[q_char(1)][q_char(0)][t_comp(1)][t_comp(0)] as i32,
    );

    let val = m.get(1, 1) + s_mat[GAP_IDX][q_char(1)][GAP_IDX][t_comp(1)] as i32;
    if val > best_e {
        best_e = val;
        best_i = 1;
        best_j = 1;
    }

    // Row 0 (j from 2 to t_len-1)
    for j in 2..t_len {
        let prev_bt = bt.get(0, j - 1);
        if prev_bt != NA_VAL {
            bt.set(
                0,
                j,
                prev_bt + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32,
            );
        }
        if q_len >= 1 && prev_bt != NA_VAL {
            m.set(
                1,
                j,
                prev_bt + s_mat[q_char(1)][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32,
            );
            let val = m.get(1, j) + s_mat[GAP_IDX][q_char(1)][GAP_IDX][t_comp(j)] as i32;
            if val > best_e {
                best_e = val;
                best_i = 1;
                best_j = j;
            }
        }
    }

    // Col 0 (i from 2 to q_len-1)
    for i in 2..q_len {
        let prev_bq = bq.get(i - 1, 0);
        if prev_bq != NA_VAL {
            bq.set(
                i,
                0,
                prev_bq + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32,
            );
        }
        if t_len >= 1 && prev_bq != NA_VAL {
            m.set(
                i,
                1,
                prev_bq + s_mat[q_char(i)][q_char(i - 1)][t_comp(1)][GAP_IDX] as i32,
            );
            let val = m.get(i, 1) + s_mat[GAP_IDX][q_char(i)][GAP_IDX][t_comp(1)] as i32;
            if val > best_e {
                best_e = val;
                best_i = i;
                best_j = 1;
            }
        }
    }

    // 2,2 Init (needs at least length 3 to have index 2)
    if q_len >= 3 && t_len >= 3 {
        let m11 = m.get(1, 1);
        if m11 != NA_VAL {
            bt.set(
                1,
                2,
                m11 + s_mat[GAP_IDX][q_char(1)][t_comp(2)][t_comp(1)] as i32,
            );
            bq.set(
                2,
                1,
                m11 + s_mat[q_char(2)][q_char(1)][GAP_IDX][t_comp(1)] as i32,
            );
            m.set(
                2,
                2,
                m11 + s_mat[q_char(2)][q_char(1)][t_comp(2)][t_comp(1)] as i32,
            );
            let val = m.get(2, 2) + s_mat[GAP_IDX][q_char(2)][GAP_IDX][t_comp(2)] as i32;
            if val > best_e {
                best_e = val;
                best_i = 2;
                best_j = 2;
            }
        }
        let m12 = m.get(1, 2);
        if m12 != NA_VAL {
            bq.set(
                2,
                2,
                m12 + s_mat[q_char(2)][q_char(1)][GAP_IDX][t_comp(2)] as i32,
            );
        }
        let m21 = m.get(2, 1);
        if m21 != NA_VAL {
            bt.set(
                2,
                2,
                m21 + s_mat[GAP_IDX][q_char(2)][t_comp(2)][t_comp(1)] as i32,
            );
        }
    }

    // Main DP
    for i in 2..q_len {
        for j in 2..t_len {
            if i == 2 && j == 2 {
                continue;
            } // Already done

            // Calc M[i,j]
            let m_diag = m.get(i - 1, j - 1);
            let bq_diag = bq.get(i - 1, j - 1);
            let bt_diag = bt.get(i - 1, j - 1);
            let mut val_m = NA_VAL;
            if m_diag != NA_VAL {
                val_m = val_m
                    .max(m_diag + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][t_comp(j - 1)] as i32);
            }
            if bq_diag != NA_VAL {
                val_m =
                    val_m.max(bq_diag + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][GAP_IDX] as i32);
            }
            if bt_diag != NA_VAL {
                val_m =
                    val_m.max(bt_diag + s_mat[q_char(i)][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32);
            }
            m.set(i, j, val_m);

            if val_m != NA_VAL {
                let curr_e = val_m + s_mat[GAP_IDX][q_char(i)][GAP_IDX][t_comp(j)] as i32;
                if curr_e > best_e {
                    best_e = curr_e;
                    best_i = i;
                    best_j = j;
                }
            }

            // Calc Bq[i,j]
            if i > 2 || (i == 2 && j > 2) {
                // check bounds logic
                let m_up = m.get(i - 1, j);
                let bq_up = bq.get(i - 1, j);
                let mut val_bq = NA_VAL;
                if m_up != NA_VAL {
                    val_bq = val_bq
                        .max(m_up + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][t_comp(j)] as i32);
                }
                if bq_up != NA_VAL {
                    val_bq = val_bq
                        .max(bq_up + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32);
                }
                bq.set(i, j, val_bq);
            }

            // Calc Bt[i,j]
            if j > 2 || (j == 2 && i > 2) {
                let m_left = m.get(i, j - 1);
                let bt_left = bt.get(i, j - 1);
                let mut val_bt = NA_VAL;
                if m_left != NA_VAL {
                    val_bt = val_bt
                        .max(m_left + s_mat[GAP_IDX][q_char(i)][t_comp(j)][t_comp(j - 1)] as i32);
                }
                if bt_left != NA_VAL {
                    val_bt = val_bt
                        .max(bt_left + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32);
                }
                bt.set(i, j, val_bt);
            }
        }
    }

    // Backtracking Logic for dp_left
    let mut i = best_i;
    let mut j = best_j;
    let mut fp = String::new();
    let mut state = DpState::Match;

    while i > 0 && j > 0 {
        let _qc = if i <= q_len { q_seq[q_start - i] } else { b'N' };

        let qc_byte = if i <= q_start {
            q_seq[q_start - i]
        } else {
            b'N'
        };
        let tc_byte = if t_start + j < t_seq.len() {
            t_seq[t_start + j]
        } else {
            b'N'
        };

        match state {
            DpState::Match => {
                let m_sc = m.get(i - 1, j - 1);
                let bq_sc = bq.get(i - 1, j - 1);
                let score = m.get(i, j);

                let val_m = if m_sc != NA_VAL {
                    m_sc + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][t_comp(j - 1)] as i32
                } else {
                    NA_FALLBACK
                };
                let val_bq = if bq_sc != NA_VAL {
                    bq_sc + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][GAP_IDX] as i32
                } else {
                    NA_FALLBACK
                };

                state = if score == val_m {
                    DpState::Match
                } else if score == val_bq {
                    DpState::GapQ
                } else {
                    DpState::GapT
                };
                fp.push(get_fingerprint_char(qc_byte, tc_byte));
                i -= 1;
                j -= 1;
            }
            DpState::GapQ => {
                let bq_up = bq.get(i - 1, j);
                let score = bq.get(i, j);
                let val_bq = if bq_up != NA_VAL {
                    bq_up + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32
                } else {
                    NA_VAL
                };

                fp.push('Q');
                state = if score == val_bq {
                    DpState::GapQ
                } else {
                    DpState::Match
                };
                i -= 1;
            }
            DpState::GapT => {
                let bt_left = bt.get(i, j - 1);
                let score = bt.get(i, j);
                let val_bt = if bt_left != NA_VAL {
                    bt_left + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32
                } else {
                    NA_VAL
                };

                fp.push('T');
                state = if score == val_bt {
                    DpState::GapT
                } else {
                    DpState::Match
                };
                j -= 1;
            }
        }
    }

    (best_e, best_i, best_j, fp)
}

fn dp_right(
    q_seq: &[u8],
    t_seq: &[u8],
    q_end: usize,
    t_end: usize,
    max_ext: usize,
) -> (i32, usize, usize, String) {
    // DP Right: Extend Query Right (towards 3'), Target Left (towards 5')
    // q_end: 3' end of seed (index incr)
    // t_end: 5' end of seed (index decr)

    // Query moves towards end.
    let q_len = (q_seq.len() - q_end).min(max_ext);
    // Target moves towards start (0).
    let t_len = (t_end + 1).min(max_ext);

    let s_mat = &DSM_T04_POS;

    // Q uses q_end + i (Moving 5' -> 3')
    // T uses t_end - j (Moving 3' -> 5' / Left)

    let q_char = |i: usize| {
        if q_end + i >= q_seq.len() {
            0
        } else {
            dsm_idx(q_seq[q_end + i])
        }
    };
    let t_comp = |j: usize| {
        if j > t_end {
            0
        } else {
            dsm_idx(t_seq[t_end - j])
        }
    };

    // Initial score from seed boundary
    let mut best_e = s_mat[q_char(0)][GAP_IDX][t_comp(0)][GAP_IDX] as i32;
    let mut best_i = 0;
    let mut best_j = 0;

    // C: if (lq <= 1 || lt <= 1) return best_e;
    if q_len <= 1 || t_len <= 1 {
        return (best_e, best_i, best_j, String::new());
    }

    let mut m = Grid::new(t_len + 1, q_len + 1, NA_VAL);
    let mut bq = Grid::new(t_len + 1, q_len + 1, NA_VAL);
    let mut bt = Grid::new(t_len + 1, q_len + 1, NA_VAL);

    m.set(0, 0, 0);

    // Init (0,1) Bt, (1,0) Bq, (1,1) M
    bt.set(0, 1, s_mat[q_char(0)][GAP_IDX][t_comp(0)][t_comp(1)] as i32);
    bq.set(1, 0, s_mat[q_char(0)][q_char(1)][t_comp(0)][GAP_IDX] as i32);
    m.set(
        1,
        1,
        s_mat[q_char(0)][q_char(1)][t_comp(0)][t_comp(1)] as i32,
    );

    let val = m.get(1, 1) + s_mat[q_char(1)][GAP_IDX][t_comp(1)][GAP_IDX] as i32;
    if val > best_e {
        best_e = val;
        best_i = 1;
        best_j = 1;
    }

    // Row 0 (j from 2 to t_len-1)
    for j in 2..t_len {
        let prev_bt = bt.get(0, j - 1);
        if prev_bt != NA_VAL {
            bt.set(
                0,
                j,
                prev_bt + s_mat[GAP_IDX][GAP_IDX][t_comp(j - 1)][t_comp(j)] as i32,
            );
        }
        if q_len >= 1 && prev_bt != NA_VAL {
            m.set(
                1,
                j,
                prev_bt + s_mat[GAP_IDX][q_char(1)][t_comp(j - 1)][t_comp(j)] as i32,
            );
            let val = m.get(1, j) + s_mat[q_char(1)][GAP_IDX][t_comp(j)][GAP_IDX] as i32;
            if val > best_e {
                best_e = val;
                best_i = 1;
                best_j = j;
            }
        }
    }

    // Col 0 (i from 2 to q_len-1)
    for i in 2..q_len {
        let prev_bq = bq.get(i - 1, 0);
        if prev_bq != NA_VAL {
            bq.set(
                i,
                0,
                prev_bq + s_mat[q_char(i - 1)][q_char(i)][GAP_IDX][GAP_IDX] as i32,
            );
        }
        if t_len >= 1 && prev_bq != NA_VAL {
            m.set(
                i,
                1,
                prev_bq + s_mat[q_char(i - 1)][q_char(i)][GAP_IDX][t_comp(1)] as i32,
            );
            let val = m.get(i, 1) + s_mat[q_char(i)][GAP_IDX][t_comp(1)][GAP_IDX] as i32;
            if val > best_e {
                best_e = val;
                best_i = i;
                best_j = 1;
            }
        }
    }

    // 2,2 Init (needs at least length 3 to have index 2)
    if q_len >= 3 && t_len >= 3 {
        let m11 = m.get(1, 1);
        if m11 != NA_VAL {
            bt.set(
                1,
                2,
                m11 + s_mat[q_char(1)][GAP_IDX][t_comp(1)][t_comp(2)] as i32,
            );
            bq.set(
                2,
                1,
                m11 + s_mat[q_char(1)][q_char(2)][t_comp(1)][GAP_IDX] as i32,
            );
            m.set(
                2,
                2,
                m11 + s_mat[q_char(1)][q_char(2)][t_comp(1)][t_comp(2)] as i32,
            );
            let val = m.get(2, 2) + s_mat[q_char(2)][GAP_IDX][t_comp(2)][GAP_IDX] as i32;
            if val > best_e {
                best_e = val;
                best_i = 2;
                best_j = 2;
            }
        }
        let m12 = m.get(1, 2);
        if m12 != NA_VAL {
            bq.set(
                2,
                2,
                m12 + s_mat[q_char(1)][q_char(2)][t_comp(2)][GAP_IDX] as i32,
            );
        }
        let m21 = m.get(2, 1);
        if m21 != NA_VAL {
            bt.set(
                2,
                2,
                m21 + s_mat[q_char(2)][GAP_IDX][t_comp(1)][t_comp(2)] as i32,
            );
        }
    }

    // Main DP
    for i in 2..q_len {
        for j in 2..t_len {
            if i == 2 && j == 2 {
                continue;
            }

            let m_diag = m.get(i - 1, j - 1);
            let bq_diag = bq.get(i - 1, j - 1);
            let bt_diag = bt.get(i - 1, j - 1);

            let mut val_m = NA_VAL;
            if m_diag != NA_VAL {
                val_m = val_m
                    .max(m_diag + s_mat[q_char(i - 1)][q_char(i)][t_comp(j - 1)][t_comp(j)] as i32);
            }
            if bq_diag != NA_VAL {
                val_m =
                    val_m.max(bq_diag + s_mat[q_char(i - 1)][q_char(i)][GAP_IDX][t_comp(j)] as i32);
            }
            if bt_diag != NA_VAL {
                val_m =
                    val_m.max(bt_diag + s_mat[GAP_IDX][q_char(i)][t_comp(j - 1)][t_comp(j)] as i32);
            }
            m.set(i, j, val_m);

            if val_m != NA_VAL {
                let curr_e = val_m + s_mat[q_char(i)][GAP_IDX][t_comp(j)][GAP_IDX] as i32;
                if curr_e > best_e {
                    best_e = curr_e;
                    best_i = i;
                    best_j = j;
                }
            }

            if i > 2 || (i == 2 && j > 2) {
                let m_up = m.get(i - 1, j);
                let bq_up = bq.get(i - 1, j);
                let mut val_bq = NA_VAL;
                if m_up != NA_VAL {
                    val_bq = val_bq
                        .max(m_up + s_mat[q_char(i - 1)][q_char(i)][t_comp(j)][GAP_IDX] as i32);
                }
                if bq_up != NA_VAL {
                    val_bq = val_bq
                        .max(bq_up + s_mat[q_char(i - 1)][q_char(i)][GAP_IDX][GAP_IDX] as i32);
                }
                bq.set(i, j, val_bq);
            }

            if j > 2 || (j == 2 && i > 2) {
                let m_left = m.get(i, j - 1);
                let bt_left = bt.get(i, j - 1);
                let mut val_bt = NA_VAL;
                if m_left != NA_VAL {
                    val_bt = val_bt
                        .max(m_left + s_mat[q_char(i)][GAP_IDX][t_comp(j - 1)][t_comp(j)] as i32);
                }
                if bt_left != NA_VAL {
                    val_bt = val_bt
                        .max(bt_left + s_mat[GAP_IDX][GAP_IDX][t_comp(j - 1)][t_comp(j)] as i32);
                }
                bt.set(i, j, val_bt);
            }
        }
    }

    // Backtracking Logic for dp_right
    let mut i = best_i;
    let mut j = best_j;
    let mut fp = String::new();
    let mut state = DpState::Match;

    while i > 0 && j > 0 {
        let qc_byte = if q_end + i < q_seq.len() {
            q_seq[q_end + i]
        } else {
            b'N'
        };
        let tc_byte = if t_end >= j { t_seq[t_end - j] } else { b'N' };

        match state {
            DpState::Match => {
                let m_sc = m.get(i - 1, j - 1);
                let bq_sc = bq.get(i - 1, j - 1);
                let score = m.get(i, j);

                let val_m = if m_sc != NA_VAL {
                    m_sc + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][t_comp(j - 1)] as i32
                } else {
                    NA_FALLBACK
                };
                let val_bq = if bq_sc != NA_VAL {
                    bq_sc + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][GAP_IDX] as i32
                } else {
                    NA_FALLBACK
                };

                state = if score == val_m {
                    DpState::Match
                } else if score == val_bq {
                    DpState::GapQ
                } else {
                    DpState::GapT
                };
                fp.push(get_fingerprint_char(qc_byte, tc_byte));
                i -= 1;
                j -= 1;
            }
            DpState::GapQ => {
                let bq_up = bq.get(i - 1, j);
                let score = bq.get(i, j);
                let val_bq = if bq_up != NA_VAL {
                    bq_up + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32
                } else {
                    NA_VAL
                };

                fp.push('Q');
                state = if score == val_bq {
                    DpState::GapQ
                } else {
                    DpState::Match
                };
                i -= 1;
            }
            DpState::GapT => {
                let bt_left = bt.get(i, j - 1);
                let score = bt.get(i, j);
                let val_bt = if bt_left != NA_VAL {
                    bt_left + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32
                } else {
                    NA_VAL
                };

                fp.push('T');
                state = if score == val_bt {
                    DpState::GapT
                } else {
                    DpState::Match
                };
                j -= 1;
            }
        }
    }

    (best_e, best_i, best_j, fp)
}

fn find_best_extent(
    q_seq: &[u8],
    t_seq: &[u8],
    seed_q: usize,
    seed_t: usize,
    seed_len: usize,
    opts: &SearchArgs,
) -> (
    usize,
    usize,
    usize,
    usize,
    usize, // l_q
    usize, // l_t
    usize, // r_q
    usize, // r_t
) {
    let max_ext = opts.max_extension as usize;
    let safe_ext = max_ext.min(MAX_DP_EXT);

    // Correct calls matching extend_seed
    let t_match_end = seed_t + seed_len - 1;

    // dp_left: Extend Q Left (5'), T Right (3')
    // Access T from t_match_end moving Right
    let (_, l_q, l_t, _) = dp_left(q_seq, t_seq, seed_q, t_match_end, safe_ext);

    // dp_right: Extend Q Right (3'), T Left (5')
    // Access T from seed_t moving Left
    let (_, r_q, r_t, _) = dp_right(q_seq, t_seq, seed_q + seed_len - 1, seed_t, safe_ext);

    // Coordinates:
    // Q Start: seed_q - l_q
    // Q End:   (seed_q + seed_len - 1) + r_q
    // T Start: seed_t - r_t (Extending Left means subtracting index)
    // T End:   (seed_t + seed_len - 1) + l_t (Extending Right means adding index)

    (
        seed_q - l_q,
        seed_q + seed_len - 1 + r_q,
        seed_t - r_t,
        seed_t + seed_len - 1 + l_t,
        l_q,
        l_t,
        r_q,
        r_t,
    )
}

/// RNA reverse complement helper
#[inline]
fn complement_rna(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'U',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'U' | b'u' | b'T' | b't' => b'A',
        _ => b'N',
    }
}

fn reverse_complement_rna(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_rna(b)).collect()
}

/// Determines the fingerprint character for a query-target nucleotide pair
/// - 'P' = Watson-Crick pair (A-U or G-C)
/// - 'W' = Wobble pair (G-U)
/// - 'U' = Unmatched/mismatch
#[inline]
fn get_fingerprint_char(qc: u8, tc: u8) -> char {
    const A: usize = 1;
    const G: usize = 2;
    const C: usize = 3;
    const U: usize = 4;

    match (dsm_idx(qc), dsm_idx(tc)) {
        (A, U) | (U, A) | (G, C) | (C, G) => 'P', // Watson-Crick pair
        (G, U) | (U, G) => 'W',                   // Wobble pair
        _ => 'U',                                 // Unmatched
    }
}

// Re-implement DP with traceback for printing
fn trace_left(
    q_seq: &[u8],
    t_seq: &[u8],
    q_start: usize,
    t_start: usize,
    best_i: usize,
    best_j: usize,
) -> (String, String, String) {
    let q_len = best_i;
    let t_len = best_j;

    let s_mat = &DSM_T04_POS;

    let q_char = |i: usize| {
        if i > q_start {
            0
        } else {
            dsm_idx(q_seq[q_start - i])
        }
    };
    let t_comp = |j: usize| {
        if t_start + j >= t_seq.len() {
            0
        } else {
            dsm_idx(t_seq[t_start + j])
        }
    };

    let rows = t_len + 1;
    let cols = q_len + 1;
    let mut m = Grid::new(rows, cols, NA_VAL);
    let mut bq = Grid::new(rows, cols, NA_VAL);
    let mut bt = Grid::new(rows, cols, NA_VAL);

    m.set(0, 0, 0);

    // Initial Scoring Logic
    if cols > 1 {
        bq.set(0, 1, s_mat[GAP_IDX][q_char(0)][t_comp(1)][t_comp(0)] as i32);
    }

    // Fill Matrix
    for j in 2..t_len {
        let prev_bt = bt.get(0, j - 1);
        if prev_bt != NA_VAL {
            bt.set(
                0,
                j,
                prev_bt + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32,
            );
        }
        if q_len >= 1 && prev_bt != NA_VAL {
            m.set(
                1,
                j,
                prev_bt + s_mat[q_char(1)][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32,
            );
        }
    }
    for i in 2..q_len {
        let prev_bq = bq.get(i - 1, 0);
        if prev_bq != NA_VAL {
            bq.set(
                i,
                0,
                prev_bq + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32,
            );
        }
        if t_len >= 1 && prev_bq != NA_VAL {
            m.set(
                i,
                1,
                prev_bq + s_mat[q_char(i)][q_char(i - 1)][t_comp(1)][GAP_IDX] as i32,
            );
        }
    }

    if q_len >= 2 && t_len >= 2 {
        let m11 = m.get(1, 1);
        if m11 != NA_VAL {
            bt.set(
                1,
                2,
                m11 + s_mat[GAP_IDX][q_char(1)][t_comp(2)][t_comp(1)] as i32,
            );
            bq.set(
                2,
                1,
                m11 + s_mat[q_char(2)][q_char(1)][GAP_IDX][t_comp(1)] as i32,
            );
            m.set(
                2,
                2,
                m11 + s_mat[q_char(2)][q_char(1)][t_comp(2)][t_comp(1)] as i32,
            );
        }
    }

    for i in 2..q_len {
        for j in 2..t_len {
            if i == 2 && j == 2 {
                continue;
            }
            let m_diag = m.get(i - 1, j - 1);
            let bq_diag = bq.get(i - 1, j - 1);
            let bt_diag = bt.get(i - 1, j - 1);
            let mut val_m = NA_VAL;
            if m_diag != NA_VAL {
                val_m = val_m
                    .max(m_diag + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][t_comp(j - 1)] as i32);
            }
            if bq_diag != NA_VAL {
                val_m =
                    val_m.max(bq_diag + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][GAP_IDX] as i32);
            }
            if bt_diag != NA_VAL {
                val_m =
                    val_m.max(bt_diag + s_mat[q_char(i)][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32);
            }
            m.set(i, j, val_m);

            if i > 2 || (i == 2 && j > 2) {
                let m_up = m.get(i - 1, j);
                let bq_up = bq.get(i - 1, j);
                let mut val_bq = NA_VAL;
                if m_up != NA_VAL {
                    val_bq = val_bq
                        .max(m_up + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][t_comp(j)] as i32);
                }
                if bq_up != NA_VAL {
                    val_bq = val_bq
                        .max(bq_up + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32);
                }
                bq.set(i, j, val_bq);
            }
            if j > 2 || (j == 2 && i > 2) {
                let m_left = m.get(i, j - 1);
                let bt_left = bt.get(i, j - 1);
                let mut val_bt = NA_VAL;
                if m_left != NA_VAL {
                    val_bt = val_bt
                        .max(m_left + s_mat[GAP_IDX][q_char(i)][t_comp(j)][t_comp(j - 1)] as i32);
                }
                if bt_left != NA_VAL {
                    val_bt = val_bt
                        .max(bt_left + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32);
                }
                bt.set(i, j, val_bt);
            }
        }
    }

    let mut i = best_i;
    let mut j = best_j;
    let mut fp = String::new();
    let mut qs = String::new();
    let mut ts = String::new();
    let mut state = DpState::Match;

    while i > 0 && j > 0 {
        let qc = if i <= q_start {
            q_seq[q_start - i]
        } else {
            b'N'
        };
        let tc = if t_start + j < t_seq.len() {
            t_seq[t_start + j]
        } else {
            b'N'
        };

        match state {
            DpState::Match => {
                let m_sc = m.get(i - 1, j - 1);
                let bq_sc = bq.get(i - 1, j - 1);
                let score = m.get(i, j);

                let from_m = if m_sc != NA_VAL {
                    m_sc + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][t_comp(j - 1)] as i32
                } else {
                    NA_FALLBACK
                };
                let from_bq = if bq_sc != NA_VAL {
                    bq_sc + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][GAP_IDX] as i32
                } else {
                    NA_FALLBACK
                };

                state = if score == from_m {
                    DpState::Match
                } else if score == from_bq {
                    DpState::GapQ
                } else {
                    DpState::GapT
                };
                fp.push(get_fingerprint_char(qc, tc));
                qs.push(qc as char);
                ts.push(tc as char);
                i -= 1;
                j -= 1;
            }
            DpState::GapQ => {
                let bq_up = bq.get(i - 1, j);
                let score = bq.get(i, j);
                let from_bq = if bq_up != NA_VAL {
                    bq_up + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32
                } else {
                    NA_VAL
                };

                fp.push('Q');
                qs.push(qc as char);
                ts.push('-');
                state = if score == from_bq {
                    DpState::GapQ
                } else {
                    DpState::Match
                };
                i -= 1;
            }
            DpState::GapT => {
                let bt_left = bt.get(i, j - 1);
                let score = bt.get(i, j);
                let from_bt = if bt_left != NA_VAL {
                    bt_left + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32
                } else {
                    NA_VAL
                };

                fp.push('T');
                qs.push('-');
                ts.push(tc as char);
                state = if score == from_bt {
                    DpState::GapT
                } else {
                    DpState::Match
                };
                j -= 1;
            }
        }
    }

    (fp, qs, ts)
}

fn trace_right(
    q_seq: &[u8],
    t_seq: &[u8],
    q_end: usize,
    t_end: usize,
    best_i: usize,
    best_j: usize,
) -> (String, String, String) {
    let q_len = best_i;
    let t_len = best_j;
    let s_mat = &DSM_T04_POS;
    let q_char = |i: usize| {
        if q_end + i >= q_seq.len() {
            0
        } else {
            dsm_idx(q_seq[q_end + i])
        }
    };
    let t_comp = |j: usize| {
        if j > t_end + 1 {
            0
        } else {
            dsm_idx(t_seq[t_end - j])
        }
    };

    let rows = t_len + 1;
    let cols = q_len + 1;
    let mut m = Grid::new(rows, cols, NA_VAL);
    let mut bq = Grid::new(rows, cols, NA_VAL);
    let mut bt = Grid::new(rows, cols, NA_VAL);

    m.set(0, 0, 0);

    if cols > 1 {
        bq.set(0, 1, s_mat[GAP_IDX][q_char(0)][t_comp(1)][t_comp(0)] as i32);
    }

    for j in 2..t_len {
        let prev_bt = bt.get(0, j - 1);
        if prev_bt != NA_VAL {
            bt.set(
                0,
                j,
                prev_bt + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32,
            );
        }
        if q_len >= 1 && prev_bt != NA_VAL {
            m.set(
                1,
                j,
                prev_bt + s_mat[q_char(1)][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32,
            );
        }
    }
    for i in 2..q_len {
        let prev_bq = bq.get(i - 1, 0);
        if prev_bq != NA_VAL {
            bq.set(
                i,
                0,
                prev_bq + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32,
            );
        }
        if t_len >= 1 && prev_bq != NA_VAL {
            m.set(
                i,
                1,
                prev_bq + s_mat[q_char(i)][q_char(i - 1)][t_comp(1)][GAP_IDX] as i32,
            );
        }
    }

    if q_len >= 2 && t_len >= 2 {
        let m11 = m.get(1, 1);
        if m11 != NA_VAL {
            bt.set(
                1,
                2,
                m11 + s_mat[GAP_IDX][q_char(1)][t_comp(2)][t_comp(1)] as i32,
            );
            bq.set(
                2,
                1,
                m11 + s_mat[q_char(2)][q_char(1)][GAP_IDX][t_comp(1)] as i32,
            );
            m.set(
                2,
                2,
                m11 + s_mat[q_char(2)][q_char(1)][t_comp(2)][t_comp(1)] as i32,
            );
        }
    }

    for i in 2..q_len {
        for j in 2..t_len {
            if i == 2 && j == 2 {
                continue;
            }
            let m_diag = m.get(i - 1, j - 1);
            let bq_diag = bq.get(i - 1, j - 1);
            let bt_diag = bt.get(i - 1, j - 1);
            let mut val_m = NA_VAL;
            if m_diag != NA_VAL {
                val_m = val_m
                    .max(m_diag + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][t_comp(j - 1)] as i32);
            }
            if bq_diag != NA_VAL {
                val_m =
                    val_m.max(bq_diag + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][GAP_IDX] as i32);
            }
            if bt_diag != NA_VAL {
                val_m =
                    val_m.max(bt_diag + s_mat[q_char(i)][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32);
            }
            m.set(i, j, val_m);

            if i > 2 || (i == 2 && j > 2) {
                let m_up = m.get(i - 1, j);
                let bq_up = bq.get(i - 1, j);
                let mut val_bq = NA_VAL;
                if m_up != NA_VAL {
                    val_bq = val_bq
                        .max(m_up + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][t_comp(j)] as i32);
                }
                if bq_up != NA_VAL {
                    val_bq = val_bq
                        .max(bq_up + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32);
                }
                bq.set(i, j, val_bq);
            }
            if j > 2 || (j == 2 && i > 2) {
                let m_left = m.get(i, j - 1);
                let bt_left = bt.get(i, j - 1);
                let mut val_bt = NA_VAL;
                if m_left != NA_VAL {
                    val_bt = val_bt
                        .max(m_left + s_mat[GAP_IDX][q_char(i)][t_comp(j)][t_comp(j - 1)] as i32);
                }
                if bt_left != NA_VAL {
                    val_bt = val_bt
                        .max(bt_left + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32);
                }
                bt.set(i, j, val_bt);
            }
        }
    }

    let mut i = best_i;
    let mut j = best_j;
    let mut fp = String::new();
    let mut qs = String::new();
    let mut ts = String::new();
    let mut state = DpState::Match;

    while i > 0 && j > 0 {
        let qc = if q_end + i < q_seq.len() {
            q_seq[q_end + i]
        } else {
            b'N'
        };
        let tc = if t_end >= j { t_seq[t_end - j] } else { b'N' };

        match state {
            DpState::Match => {
                let m_sc = m.get(i - 1, j - 1);
                let bq_sc = bq.get(i - 1, j - 1);
                let score = m.get(i, j);

                let from_m = if m_sc != NA_VAL {
                    m_sc + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][t_comp(j - 1)] as i32
                } else {
                    NA_FALLBACK
                };
                let from_bq = if bq_sc != NA_VAL {
                    bq_sc + s_mat[q_char(i)][q_char(i - 1)][t_comp(j)][GAP_IDX] as i32
                } else {
                    NA_FALLBACK
                };

                state = if score == from_m {
                    DpState::Match
                } else if score == from_bq {
                    DpState::GapQ
                } else {
                    DpState::GapT
                };
                fp.push(get_fingerprint_char(qc, tc));
                qs.push(qc as char);
                ts.push(tc as char);
                i -= 1;
                j -= 1;
            }
            DpState::GapQ => {
                let bq_up = bq.get(i - 1, j);
                let score = bq.get(i, j);
                let from_bq = if bq_up != NA_VAL {
                    bq_up + s_mat[q_char(i)][q_char(i - 1)][GAP_IDX][GAP_IDX] as i32
                } else {
                    NA_VAL
                };

                fp.push('Q');
                qs.push(qc as char);
                ts.push('-');
                state = if score == from_bq {
                    DpState::GapQ
                } else {
                    DpState::Match
                };
                i -= 1;
            }
            DpState::GapT => {
                let bt_left = bt.get(i, j - 1);
                let score = bt.get(i, j);
                let from_bt = if bt_left != NA_VAL {
                    bt_left + s_mat[GAP_IDX][GAP_IDX][t_comp(j)][t_comp(j - 1)] as i32
                } else {
                    NA_VAL
                };

                fp.push('T');
                qs.push('-');
                ts.push(tc as char);
                state = if score == from_bt {
                    DpState::GapT
                } else {
                    DpState::Match
                };
                j -= 1;
            }
        }
    }

    (fp, qs, ts)
}

#[allow(clippy::too_many_arguments)]
fn print_detailed_output(
    w: &mut dyn std::io::Write,
    q_id: &str,
    q_seq: &[u8],
    q_start: usize,
    q_end: usize,
    t_id: &str,
    t_seq: &[u8],
    t_start: usize,
    t_end: usize,
    out_ts: usize,
    out_te: usize,
    score: f64,
    strand: char,
    l_q: usize,
    l_t: usize,
    r_q: usize,
    r_t: usize,
) -> std::io::Result<()> {
    // Traceback logic to generate strings
    // Recover seed positions:
    let seed_q = q_start + l_q;
    let seed_t = t_start + r_t;
    let seed_len = (q_end - r_q) - seed_q + 1;
    let t_match_end = seed_t + seed_len - 1;

    let (fp_l, _qs_l, ts_l) = trace_left(q_seq, t_seq, seed_q, t_match_end, l_q, l_t);
    let (fp_r, _qs_r, ts_r) = trace_right(q_seq, t_seq, seed_q + seed_len - 1, seed_t, r_q, r_t);

    // Build Seed Strings (assumed exact/wobble match check done by seed index)
    // We construct them here.
    let mut fp_seed = String::with_capacity(seed_len);
    let mut qs_seed = String::with_capacity(seed_len);
    let mut ts_seed = String::with_capacity(seed_len);

    for k in 0..seed_len {
        let qc = q_seq[seed_q + k];
        let tc = t_seq[t_match_end - k]; // Antiparallel match
        qs_seed.push(qc as char);
        ts_seed.push(tc as char);
        fp_seed.push(get_fingerprint_char(qc, tc));
    }

    // Combine
    let full_fp = format!(
        "{}{}{}",
        fp_l,
        fp_seed,
        fp_r.chars().rev().collect::<String>()
    );
    // Note: C output doesn't seem to include Query String in the detailed (one-line) output?
    // User format request: QueryID \t QCompStart \t QCompEnd \t TargetID \t TStart \t TEnd \t Strand \t Energy \t IntString \t TargetString \t Flank5 \t Flank3
    // So we need IntString and TargetString.

    let full_ts = format!(
        "{}{}{}",
        ts_l,
        ts_seed,
        ts_r.chars().rev().collect::<String>()
    );

    // Normalize T -> U for output strings
    let norm_fp = full_fp; // P/W/U/y/x don't need normalization
    let norm_ts = full_ts.replace('T', "U").replace('t', "u");

    // Flanking - Normalize too
    // T starts at `t_start` (lowest index) and ends at `t_end` (highest index).
    // 5' Flank (Upstream relative to T sequence): t_start - 20 .. t_start
    let ctx_len = 20;
    let t_5_start = t_start.saturating_sub(ctx_len);
    let ctx_5 = String::from_utf8_lossy(&t_seq[t_5_start..t_start])
        .replace('T', "U")
        .replace('t', "u")
        .chars()
        .rev()
        .collect::<String>();

    // 3' Flank (Downstream): t_end + 1 .. t_end + 1 + 20
    let t_3_start = t_end + 1;
    let t_3_end = (t_3_start + ctx_len).min(t_seq.len());
    let ctx_3 = if t_3_start < t_seq.len() {
        String::from_utf8_lossy(&t_seq[t_3_start..t_3_end])
            .replace('T', "U")
            .replace('t', "u")
    } else {
        String::new()
    };

    // One-line TSV output to match C
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}",
        q_id.split_whitespace().next().unwrap_or(q_id), // Take first word of ID
        q_start + 1,                                    // 1-based (Query seems to use 1-based in C)
        q_end + 1,                                      // 1-based
        t_id.split_whitespace().next().unwrap_or(t_id), // Sanitize Target ID too
        out_ts,                                         // 0-based
        out_te,                                         // 0-based
        strand,
        score,
        norm_fp,
        norm_ts,
        ctx_5,
        ctx_3
    )?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]

    fn test_dp_left_returns_score() {
        // Verifies dp_left runs without panicking and returns a score tuple
        let q = b"AAAA";
        let t = b"UUUU";
        let (score, i, j, _) = dp_left(q, t, 3, 3, 10);
        // The actual score depends on the scoring matrix and sequence alignment
        // For now, just verify the function returns without panicking
        println!("dp_left score: {}, i: {}, j: {}", score, i, j);
    }

    #[test]
    fn test_dp_right_returns_score() {
        // Verifies dp_right runs without panicking and returns a score tuple
        let q = b"AAAA";
        let t = b"UUUU";
        let (score, i, j, _) = dp_right(q, t, 0, 0, 10);
        println!("dp_right score: {}, i: {}, j: {}", score, i, j);
    }
}
