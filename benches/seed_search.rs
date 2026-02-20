use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use risearch::config::{MismatchSpec, SeedConfig, SeedSpec};
use risearch::index::sa::SuffixArray;
use risearch::seed::searcher::{SeedMatch, SeedSearcher};
use risearch::seq::Sequence;
use risearch::types::Base;

// ============================================================================
// SIMPLE LINEAR CONGRUENTIAL GENERATOR (LCG)
// ============================================================================

struct SimpleLcg {
    state: u64,
}

impl SimpleLcg {
    const A: u64 = 6364136223846793005;
    const C: u64 = 1442695040888963407;

    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(Self::A).wrapping_add(Self::C);
        self.state
    }

    fn next_base(&mut self) -> Base {
        match (self.next() >> 32) % 4 {
            0 => Base::A,
            1 => Base::G,
            2 => Base::C,
            _ => Base::U,
        }
    }
}

// ============================================================================
// SEQUENCE GENERATION
// ============================================================================

fn generate_sequence(len: usize, seed: u64) -> Sequence {
    let mut rng = SimpleLcg::new(seed);
    let bases: Vec<Base> = (0..len).map(|_| rng.next_base()).collect();
    Sequence::from(bases)
}

// ============================================================================
// SETUP HELPERS
// ============================================================================

/// Build query SA data: RC the query, build SA on RC, append sentinel Gap.
fn build_query(q_len: usize, rng_seed: u64) -> (SuffixArray, Sequence) {
    let query = generate_sequence(q_len, rng_seed);
    let query_rc = query.reverse_complement();
    let reverse_sa = SuffixArray::try_from(&query_rc).expect("query SA");
    let mut rc_bases: Vec<Base> = query_rc.iter().copied().collect();
    rc_bases.push(Base::Gap);
    let padded_query_rc = Sequence::from(rc_bases);
    (reverse_sa, padded_query_rc)
}

/// Build target SA data: combined fwd ++ Gap ++ RC, build SA, append sentinel Gap.
fn build_target(t_len: usize, rng_seed: u64) -> (SuffixArray, Sequence) {
    let target = generate_sequence(t_len, rng_seed);
    let target_rc = target.reverse_complement();
    let mut combined_bases: Vec<Base> = Vec::with_capacity(2 * t_len + 2);
    combined_bases.extend_from_slice(&target);
    combined_bases.push(Base::Gap);
    combined_bases.extend_from_slice(&target_rc);
    let combined_for_sa = Sequence::from(combined_bases.clone());
    let combined_sa = SuffixArray::try_from(&combined_for_sa).expect("target SA");
    combined_bases.push(Base::Gap);
    let padded_combined = Sequence::from(combined_bases);
    (combined_sa, padded_combined)
}

/// Build a complete (query_sa, query_rc, target_sa, target_combined) pair.
fn build_search_pair(
    q_len: usize,
    t_len: usize,
    rng_seed: u64,
) -> (SuffixArray, Sequence, SuffixArray, Sequence) {
    let (q_sa, q_seq) = build_query(q_len, rng_seed);
    let (t_sa, t_seq) = build_target(t_len, rng_seed.wrapping_add(1));
    (q_sa, q_seq, t_sa, t_seq)
}

// ============================================================================
// BENCHMARK GROUPS
// ============================================================================

/// Partition scaling: fixed query (22nt), exact matching, varying target SA size.
fn bench_seed_exact(c: &mut Criterion) {
    let mut group = c.benchmark_group("seed_exact");

    let seed_config =
        SeedConfig::with_wobble(SeedSpec::LengthOnly(7), MismatchSpec::exact(), false);

    for t_len in [1_000, 10_000, 100_000] {
        group.bench_with_input(BenchmarkId::new("t_len", t_len), &t_len, |b, &t_len| {
            let (q_sa, q_seq, t_sa, t_seq) = build_search_pair(22, t_len, 42);
            let mut results: Vec<SeedMatch> = Vec::with_capacity(4096);

            b.iter(|| {
                results.clear();
                let searcher = SeedSearcher::new(
                    black_box(&q_sa),
                    black_box(&q_seq),
                    black_box(&t_sa),
                    black_box(&t_seq),
                    black_box(&seed_config),
                );
                searcher.search_length_range(7, 22, &mut results);
                black_box(results.len());
            });
        });
    }

    group.finish();
}

/// Mismatch branching overhead: fixed query/target, varying max mismatches.
fn bench_seed_mismatch(c: &mut Criterion) {
    let mut group = c.benchmark_group("seed_mismatch");

    for max_mm in [0, 1, 2] {
        let seed_config = SeedConfig::with_wobble(
            SeedSpec::LengthOnly(7),
            MismatchSpec::new(max_mm, 2, 2),
            false,
        );

        group.bench_with_input(BenchmarkId::new("max_mm", max_mm), &max_mm, |b, _| {
            let (q_sa, q_seq, t_sa, t_seq) = build_search_pair(22, 10_000, 42);
            let mut results: Vec<SeedMatch> = Vec::with_capacity(4096);

            b.iter(|| {
                results.clear();
                let searcher = SeedSearcher::new(
                    black_box(&q_sa),
                    black_box(&q_seq),
                    black_box(&t_sa),
                    black_box(&t_seq),
                    black_box(&seed_config),
                );
                searcher.search_length_range(7, 22, &mut results);
                black_box(results.len());
            });
        });
    }

    group.finish();
}

/// Representative batched workload: 10 queries (22nt) vs one 100K target,
/// with 1 mismatch (prefix=2, suffix=2) and wobble pairs enabled.
fn bench_seed_realistic(c: &mut Criterion) {
    let mut group = c.benchmark_group("seed_realistic");
    group.sample_size(10);

    let seed_config =
        SeedConfig::with_wobble(SeedSpec::LengthOnly(7), MismatchSpec::new(1, 2, 2), true);

    group.bench_function("10q_x_100k", |b| {
        let (t_sa, t_seq) = build_target(100_000, 99);
        let queries: Vec<_> = (0..10).map(|i| build_query(22, 100 + i as u64)).collect();
        let mut results: Vec<SeedMatch> = Vec::with_capacity(16384);

        b.iter(|| {
            for (q_sa, q_seq) in &queries {
                results.clear();
                let searcher = SeedSearcher::new(
                    black_box(q_sa),
                    black_box(q_seq),
                    black_box(&t_sa),
                    black_box(&t_seq),
                    black_box(&seed_config),
                );
                searcher.search_length_range(7, 22, &mut results);
                black_box(results.len());
            }
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_seed_exact,
    bench_seed_mismatch,
    bench_seed_realistic
);
criterion_main!(benches);
