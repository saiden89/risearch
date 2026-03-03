use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use risearch::config::{MismatchSpec, SeedConfig, SeedSpec};
use risearch::index::sa::SuffixArray;
use risearch::index::store::SA_CHAR_PADDING;
use risearch::seed::searcher::{SeedMatch, SeedSearcher};
use risearch::seq::Sequence;
use risearch::types::Base;

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

fn generate_sequence(len: usize, seed: u64) -> Sequence {
    let mut rng = SimpleLcg::new(seed);
    let bases: Vec<Base> = (0..len).map(|_| rng.next_base()).collect();
    Sequence::from(bases)
}

fn build_padded_query(q_len: usize, rng_seed: u64) -> (Vec<u64>, Vec<Base>, usize, usize) {
    let query = generate_sequence(q_len, rng_seed);
    let query_sa = SuffixArray::try_from(&query[..]).expect("query SA");

    let q_sa_len = query_sa.len();
    let mut padded_sa = query_sa.into_inner();
    padded_sa.resize(q_sa_len + SA_CHAR_PADDING, 0u64);

    let mut padded_seq: Vec<Base> = query.iter().copied().collect();
    padded_seq.resize(padded_seq.len() + SA_CHAR_PADDING, Base::Gap);

    (padded_sa, padded_seq, 0, q_sa_len)
}

fn build_padded_target(t_len: usize, rng_seed: u64) -> (Vec<u64>, Vec<Base>, usize) {
    let target = generate_sequence(t_len, rng_seed);
    let target_rc = target.reverse_complement();

    let mut combined_bases: Vec<Base> = Vec::with_capacity(2 * t_len + 2 + SA_CHAR_PADDING);
    combined_bases.extend_from_slice(target.as_ref());
    combined_bases.push(Base::Gap);
    combined_bases.extend_from_slice(target_rc.as_ref());

    let combined_sa = SuffixArray::try_from(combined_bases.as_slice()).expect("target SA");

    let t_sa_len = combined_sa.len();
    let mut padded_sa = combined_sa.into_inner();
    padded_sa.resize(t_sa_len + SA_CHAR_PADDING, 0u64);

    // +1 for trailing Gap separator after rc sequence, then SA_CHAR_PADDING sentinels
    combined_bases.resize(combined_bases.len() + 1 + SA_CHAR_PADDING, Base::Gap);
    (padded_sa, combined_bases, t_sa_len)
}

#[allow(clippy::type_complexity)]
fn build_search_pair(
    q_len: usize,
    t_len: usize,
    rng_seed: u64,
) -> (
    Vec<u64>,
    Vec<Base>,
    usize,
    usize,
    Vec<u64>,
    Vec<Base>,
    usize,
) {
    let (q_sa, q_seq, q_sa_start, q_sa_len) = build_padded_query(q_len, rng_seed);
    let (t_sa, t_seq, t_sa_len) = build_padded_target(t_len, rng_seed.wrapping_add(1));
    (q_sa, q_seq, q_sa_start, q_sa_len, t_sa, t_seq, t_sa_len)
}

fn bench_seed_exact(c: &mut Criterion) {
    let mut group = c.benchmark_group("seed_exact");
    let seed_config =
        SeedConfig::with_wobble(SeedSpec::LengthOnly(7), MismatchSpec::exact(), false);

    for t_len in [1_000, 10_000, 100_000] {
        group.bench_with_input(BenchmarkId::from_parameter(t_len), &t_len, |b, &t_len| {
            let (q_sa, q_seq, q_sa_start, q_sa_len, t_sa, t_seq, t_sa_len) =
                build_search_pair(22, t_len, 42);
            let mut results: Vec<SeedMatch> = Vec::with_capacity(4096);

            b.iter(|| {
                results.clear();
                let searcher = SeedSearcher::new(
                    black_box(q_sa.as_slice()),
                    black_box(q_seq.as_slice()),
                    black_box(q_sa_start),
                    black_box(q_sa_len),
                    black_box(t_sa.as_slice()),
                    black_box(t_seq.as_slice()),
                    black_box(t_sa_len),
                    black_box(&seed_config),
                );
                searcher.search_length_range(7, 22, &mut results);
                black_box(results.len());
            });
        });
    }

    group.finish();
}

fn bench_seed_mismatch(c: &mut Criterion) {
    let mut group = c.benchmark_group("seed_mismatch");

    for max_mm in [0, 1, 2] {
        let seed_config = SeedConfig::with_wobble(
            SeedSpec::LengthOnly(7),
            MismatchSpec::new(max_mm, 2, 2),
            false,
        );

        group.bench_with_input(BenchmarkId::new("max_mm", max_mm), &max_mm, |b, _| {
            let (q_sa, q_seq, q_sa_start, q_sa_len, t_sa, t_seq, t_sa_len) =
                build_search_pair(22, 10_000, 42);
            let mut results: Vec<SeedMatch> = Vec::with_capacity(4096);

            b.iter(|| {
                results.clear();
                let searcher = SeedSearcher::new(
                    black_box(q_sa.as_slice()),
                    black_box(q_seq.as_slice()),
                    black_box(q_sa_start),
                    black_box(q_sa_len),
                    black_box(t_sa.as_slice()),
                    black_box(t_seq.as_slice()),
                    black_box(t_sa_len),
                    black_box(&seed_config),
                );
                searcher.search_length_range(7, 22, &mut results);
                black_box(results.len());
            });
        });
    }

    group.finish();
}

fn bench_seed_realistic(c: &mut Criterion) {
    let mut group = c.benchmark_group("seed_realistic");
    group.sample_size(10);

    let seed_config =
        SeedConfig::with_wobble(SeedSpec::LengthOnly(7), MismatchSpec::new(1, 2, 2), true);

    group.bench_function("10q_x_100k", |b| {
        let (t_sa, t_seq, t_sa_len) = build_padded_target(100_000, 99);
        let queries: Vec<_> = (0..10)
            .map(|i| build_padded_query(22, 100 + i as u64))
            .collect();
        let mut results: Vec<SeedMatch> = Vec::with_capacity(16384);

        b.iter(|| {
            for (q_sa, q_seq, q_sa_start, q_sa_len) in &queries {
                results.clear();
                let searcher = SeedSearcher::new(
                    black_box(q_sa.as_slice()),
                    black_box(q_seq.as_slice()),
                    black_box(*q_sa_start),
                    black_box(*q_sa_len),
                    black_box(t_sa.as_slice()),
                    black_box(t_seq.as_slice()),
                    black_box(t_sa_len),
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
