//! Correctness tests for unlimited seed extension (`-l -1`).
//!
//! risearch2 (C) cannot serve as an oracle: it clamps `-l` with
//! `MAX(0, atoi(optarg))`, so `-l -1` silently becomes `0` (seed-only). The
//! unlimited window is a risearch3-only feature, so these tests use risearch3
//! itself as the oracle via metamorphic relations.
//!
//! Note on what is *not* tested: `-l -1` is **not** equivalent to a large fixed
//! window such as `-l 256`. Unlimited sizes each extension window to the query
//! bases available on that side (capped at the `dp::MAX_EXT` = 256 buffer
//! ceiling), whereas a fixed `-l k` permits up to `k` of extension including
//! large target-side bulges. The two therefore diverge on real data, so there
//! is no fixed-window oracle to compare against. Instead we assert the
//! feature's actual guarantees:
//!
//! 1. Reach — extension crosses past the exact-match seed to the query end
//!    (`unlimited_extends_past_the_exact_seed_to_the_query_end`).
//! 2. Rejection — a query longer than the ceiling is refused up front rather
//!    than silently clamped (`unlimited_rejects_a_query_longer_than_the_ceiling`).
//! 3. Determinism — repeated runs agree (`unlimited_is_deterministic`).
//!
//! Extension only does work when the optimal duplex extends past the exact
//! seed, so the reach test breaks exact complementarity with a single mismatch
//! near the 3' end: the seeder stops there, and only DP extension can bridge it.

use std::cmp::Ordering;
use std::fs;

use risearch::config::{
    ExtendConfig, FilterConfig, OutputCompression, OutputConfig, OutputFormat, ScoreConfig,
    SeedConfig,
};
use risearch::{
    run_search, DsmId, Energy, QueryRegistry, SearchConfig, SearchHit, Strand, TargetRegistry,
    VecSink,
};

// =============================================================================
// Harness
// =============================================================================

/// A search config in which only `max_extension` varies. The extension penalty
/// is zero so that added base pairs are favorable and extension actually runs
/// (a high per-nucleotide penalty would suppress it entirely).
fn config(max_extension: i32) -> SearchConfig {
    SearchConfig {
        seed: SeedConfig {
            seed_start: None,
            seed_end: None,
            seed_length: Some(7),
            seed_wobble: true,
            max_mismatches: 0,
            min_prefix_matches: 1,
            min_suffix_matches: 0,
        },
        score: ScoreConfig {
            dsm_id: DsmId::from("t04"),
            penalty: Energy::from_kcal(0.0),
            temperature: 37,
        },
        extend: ExtendConfig {
            max_extension,
            build_alignment: true,
        },
        filter: FilterConfig {
            delta_g: Energy::from_kcal(-8.0),
            seed_energy: Energy::from_kcal(0.0),
            no_max_prune: false,
            no_dedup: false,
        },
        output: OutputConfig {
            format: OutputFormat::Detailed,
            compress: OutputCompression::None,
            multifile: false,
        },
    }
}

/// Search a single raw `query` against a single raw `target` at the given
/// window, returning hits in canonical order (or the search error). Output
/// order is not stable (hits drain from a `HashMap`), so results are sorted
/// before comparison. Index build / query load are setup and unwrap; only the
/// search itself surfaces as `Err`.
fn search_or_err(query: &str, target: &str, max_extension: i32) -> Result<Vec<SearchHit>, String> {
    let cfg = config(max_extension);
    let tmp = tempfile::tempdir().unwrap();
    let query_fa = tmp.path().join("query.fa");
    let target_fa = tmp.path().join("target.fa");
    let idx = tmp.path().join("target.idx");
    fs::write(&query_fa, format!(">query\n{query}\n")).unwrap();
    fs::write(&target_fa, format!(">target\n{target}\n")).unwrap();

    TargetRegistry::build(&target_fa, &idx, None).unwrap();
    let store = TargetRegistry::open(&idx).unwrap();
    let queries = QueryRegistry::from_fasta(&query_fa, &cfg.seed).unwrap();

    let sink = VecSink::default();
    run_search(&queries, &store, &cfg, &sink).map_err(|e| e.to_string())?;

    let mut hits = sink.into_hits();
    hits.sort_by(cmp_hit);
    Ok(hits)
}

fn search_seqs(query: &str, target: &str, max_extension: i32) -> Vec<SearchHit> {
    search_or_err(query, target, max_extension).unwrap()
}

fn complement(base: char) -> char {
    match base {
        'A' => 'U',
        'U' => 'A',
        'G' => 'C',
        'C' => 'G',
        other => other,
    }
}

/// Reverse complement, so `revcomp(q)` forms a full antiparallel duplex with `q`.
fn revcomp(seq: &str) -> String {
    seq.chars().rev().map(complement).collect()
}

/// Deterministic pseudo-random RNA (LCG), for inputs too long to write by hand.
fn pseudo_random_rna(len: usize, seed: u64) -> String {
    let mut state = seed;
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            match (state >> 33) % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'U',
            }
        })
        .collect()
}

/// Observable projection of a hit. `Energy` is `Eq` (backed by an `i32`), so
/// equality here is exact, not float-fuzzy.
#[derive(Debug, PartialEq)]
struct Observed {
    query_idx: usize,
    target_idx: usize,
    q_start: usize,
    q_end: usize,
    t_start: usize,
    t_end: usize,
    strand: Strand,
    energy: Energy,
}

impl From<&SearchHit> for Observed {
    fn from(h: &SearchHit) -> Self {
        Self {
            query_idx: h.query_idx,
            target_idx: h.target_idx,
            q_start: h.q_start,
            q_end: h.q_end,
            t_start: h.t_start,
            t_end: h.t_end,
            strand: h.strand,
            energy: h.energy,
        }
    }
}

fn observed(hits: &[SearchHit]) -> Vec<Observed> {
    hits.iter().map(Observed::from).collect()
}

fn strand_rank(s: Strand) -> u8 {
    match s {
        Strand::Forward => 0,
        Strand::Reverse => 1,
    }
}

fn cmp_hit(a: &SearchHit, b: &SearchHit) -> Ordering {
    let coords = |h: &SearchHit| {
        (
            h.query_idx,
            h.target_idx,
            h.q_start,
            h.q_end,
            h.t_start,
            h.t_end,
            strand_rank(h.strand),
        )
    };
    coords(a)
        .cmp(&coords(b))
        .then_with(|| a.energy.to_kcal().total_cmp(&b.energy.to_kcal()))
}

// =============================================================================
// Tests
// =============================================================================

/// Reach: extension crosses a mismatch the exact-match seeder cannot, reaching
/// the query 3' end — something seed-only (`-l 0`) never does.
#[test]
fn unlimited_extends_past_the_exact_seed_to_the_query_end() {
    // Perfect 30-nt duplex, then one base near the 3' end mutated to break exact
    // matching there. The seeder stops before the mismatch; only DP extension
    // can bridge it and pair the final bases.
    const N: usize = 30;
    const BREAK: usize = 27; // mismatch position; indices 28,29 stay complementary

    let mut query: Vec<char> = pseudo_random_rna(N, 0x07).chars().collect();
    let target = revcomp(&query.iter().collect::<String>());
    // In an antiparallel duplex query index i pairs target index N-1-i. Setting
    // the query base equal to its opposing target base guarantees a mismatch
    // (identical bases never pair), without touching neighboring positions.
    let opposing = target.as_bytes()[N - 1 - BREAK] as char;
    query[BREAK] = opposing;
    let query: String = query.into_iter().collect();

    let seed_only = search_seqs(&query, &target, 0);
    let unlimited = search_seqs(&query, &target, -1);

    assert!(
        seed_only.iter().all(|h| h.q_end < N - 1),
        "seed-only (`-l 0`) must stop at the mismatch, not reach the 3' end; got {:#?}",
        observed(&seed_only)
    );
    assert!(
        unlimited
            .iter()
            .any(|h| h.q_start == 0 && h.q_end == N - 1 && h.energy.to_kcal() < 0.0),
        "unlimited extension must cross the mismatch and reach the query 3' end; got {:#?}",
        observed(&unlimited)
    );
}

/// Rejection: a query longer than the 256-nt buffer ceiling (`dp::MAX_EXT`) is
/// refused up front rather than silently clamped — `-l -1` cannot honor
/// "span the whole query" past the cap.
#[test]
fn unlimited_rejects_a_query_longer_than_the_ceiling() {
    const N: usize = 300; // > 256

    let query = pseudo_random_rna(N, 0x5151_2323);
    let target = revcomp(&query);

    let err = search_or_err(&query, &target, -1)
        .expect_err("`-l -1` must reject a query longer than the 256-nt cap");

    assert!(
        err.contains("cannot extend across it"),
        "unexpected error message: {err}"
    );
}

/// Determinism: `-l -1` yields identical results across repeated runs.
#[test]
fn unlimited_is_deterministic() {
    let query = pseudo_random_rna(60, 0x0B);
    let target = revcomp(&query);

    let first = search_seqs(&query, &target, -1);
    let second = search_seqs(&query, &target, -1);

    assert_eq!(
        observed(&first),
        observed(&second),
        "`-l -1` must be deterministic across runs"
    );
}
