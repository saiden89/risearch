use std::fs;
use std::path::Path;

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use risearch::config::{
    ExtendConfig, FilterConfig, Matrix, MismatchSpec, OutputConfig, OutputFormat, ScoreConfig,
    SearchArgs, SeedConfig, SeedSpec,
};
use risearch::index::sa::SuffixArray;
use risearch::index::store::SA_CHAR_PADDING;
use risearch::registry::{Query, QueryRegistry};
use risearch::search::run_search;
use risearch::seed::searcher::{SeedMatch, SeedSearcher};
use risearch::seq::Sequence;
use risearch::types::Base;
use risearch::TargetStore;
use tempfile::TempDir;

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

struct PreparedQuery {
    padded_q_sa: Vec<u64>,
    padded_q_seq: Vec<Base>,
    q_sa_start: usize,
    q_sa_len: usize,
    min_len: usize,
    max_len: usize,
}

struct ProductionSeedDataset {
    _tmpdir: TempDir,
    queries: QueryRegistry,
    store: TargetStore,
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

    // Mirror TargetStore layout:
    // fwd_comp[seq_len] + Gap + rc_comp[seq_len] + Gap
    let mut combined_bases: Vec<Base> = Vec::with_capacity(2 * t_len + 2 + SA_CHAR_PADDING);
    combined_bases.extend(target.iter().copied().map(Base::complement));
    combined_bases.push(Base::Gap);
    combined_bases.extend(target[..].iter().rev().copied());
    combined_bases.push(Base::Gap);

    let combined_sa = SuffixArray::try_from(combined_bases.as_slice()).expect("target SA");

    let t_sa_len = combined_sa.len();
    let mut padded_sa = combined_sa.into_inner();
    padded_sa.resize(t_sa_len + SA_CHAR_PADDING, 0u64);

    combined_bases.resize(combined_bases.len() + SA_CHAR_PADDING, Base::Gap);
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

fn base_to_byte(base: Base) -> u8 {
    match base {
        Base::A => b'A',
        Base::C => b'C',
        Base::G => b'G',
        Base::U => b'U',
        Base::N => b'N',
        Base::Gap => b'-',
    }
}

fn write_fasta(path: &Path, prefix: &str, seqs: &[Sequence]) {
    let mut buf = Vec::new();
    for (idx, seq) in seqs.iter().enumerate() {
        buf.extend_from_slice(format!(">{prefix}{idx}\n").as_bytes());
        for &base in seq.iter() {
            buf.push(base_to_byte(base));
        }
        buf.push(b'\n');
    }
    fs::write(path, buf).expect("write FASTA");
}

fn build_production_dataset(
    query_count: usize,
    query_len: usize,
    target_len: usize,
    seed_config: &SeedConfig,
) -> ProductionSeedDataset {
    let tmpdir = TempDir::new().expect("tempdir");
    let queries_path = tmpdir.path().join("queries.fa");
    let targets_path = tmpdir.path().join("targets.fa");
    let index_path = tmpdir.path().join("targets.rsidx");

    let queries: Vec<_> = (0..query_count)
        .map(|i| generate_sequence(query_len, 1_000 + i as u64))
        .collect();
    let targets = vec![generate_sequence(target_len, 9_999)];

    write_fasta(&queries_path, "q", &queries);
    write_fasta(&targets_path, "t", &targets);

    let queries = QueryRegistry::from_fasta(&queries_path, seed_config).expect("query registry");
    TargetStore::build_from_fasta(&targets_path, &index_path).expect("build target index");
    let store = TargetStore::open(&index_path).expect("open target index");

    ProductionSeedDataset {
        _tmpdir: tmpdir,
        queries,
        store,
    }
}

fn prepare_query_for_seed_search(
    query: &Query,
    seed_config: &SeedConfig,
) -> Option<PreparedQuery> {
    let interval = query.seed_interval();
    let min_len = seed_config
        .seed
        .normalize(query.sequence().len())
        .expect("prepared query must be used with a compatible seed config")
        .2;
    let max_len = interval.end.saturating_sub(interval.start);
    let q_sa_len = query.sa().len();
    if q_sa_len == 0 {
        return None;
    }

    let mut padded_q_sa = query.sa().to_vec();
    padded_q_sa.resize(q_sa_len + SA_CHAR_PADDING, 0u64);
    let mut padded_q_seq = query.seed_sequence().as_slice().to_vec();
    padded_q_seq.resize(query.seed_sequence().len() + SA_CHAR_PADDING, Base::Gap);

    Some(PreparedQuery {
        padded_q_sa,
        padded_q_seq,
        q_sa_start: 0,
        q_sa_len,
        min_len,
        max_len,
    })
}

fn make_search_args(seed_config: &SeedConfig) -> SearchArgs {
    SearchArgs {
        seed: seed_config.clone(),
        score: ScoreConfig {
            matrix: Matrix::T04,
            penalty: 0.0,
            matrix2: None,
            matpath: None,
            temperature: None,
            weights: None,
        },
        extend: ExtendConfig {
            max_extension: 0,
            band: None,
        },
        filter: FilterConfig {
            delta_g: f64::NEG_INFINITY,
            seed_energy: 0.0,
            no_max_prune: false,
        },
        output: OutputConfig {
            format: OutputFormat::Minimal,
            compress: None,
            level: None,
            multifile: false,
        },
        one_vs_one: false,
        three_prime_match: None,
        five_prime_match: None,
    }
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

    for max_mm in [1, 2, 3, 4, 5] {
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

fn bench_seed_prod_shaped_mismatch(c: &mut Criterion) {
    let mut group = c.benchmark_group("seed_prod_shaped_mismatch");
    group.sample_size(10);

    for max_mm in [1, 2, 3, 4, 5] {
        let seed_config = SeedConfig::with_wobble(
            SeedSpec::LengthOnly(7),
            MismatchSpec::new(max_mm, 2, 2),
            true,
        );
        let dataset = build_production_dataset(10, 22, 100_000, &seed_config);
        let global = dataset.store.global_view();
        let prepared_queries: Vec<_> = dataset
            .queries
            .entries()
            .iter()
            .filter_map(|query| prepare_query_for_seed_search(query, &seed_config))
            .collect();
        let mut results: Vec<SeedMatch> = Vec::with_capacity(16_384);

        group.bench_with_input(
            BenchmarkId::new("prepared_queries_10q_x_100k", max_mm),
            &max_mm,
            |b, _| {
                b.iter(|| {
                    for prepared in &prepared_queries {
                        results.clear();
                        let searcher = SeedSearcher::new(
                            black_box(prepared.padded_q_sa.as_slice()),
                            black_box(prepared.padded_q_seq.as_slice()),
                            black_box(prepared.q_sa_start),
                            black_box(prepared.q_sa_len),
                            black_box(global.combined_sa),
                            black_box(global.combined_seq),
                            black_box(global.sa_real_len),
                            black_box(&seed_config),
                        );
                        searcher.search_length_range(
                            prepared.min_len,
                            prepared.max_len,
                            &mut results,
                        );
                        black_box(results.len());
                    }
                });
            },
        );
    }

    group.finish();
}

fn bench_seed_realistic(c: &mut Criterion) {
    let mut group = c.benchmark_group("seed_realistic");
    group.sample_size(10);

    let seed_config =
        SeedConfig::with_wobble(SeedSpec::LengthOnly(7), MismatchSpec::new(1, 2, 2), true);

    group.bench_function("synthetic_searcher_10q_x_100k", |b| {
        let (t_sa, t_seq, t_sa_len) = build_padded_target(100_000, 99);
        let queries: Vec<_> = (0..10)
            .map(|i| build_padded_query(22, 100 + i as u64))
            .collect();
        let mut results: Vec<SeedMatch> = Vec::with_capacity(16_384);

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

fn bench_seed_prod_shaped_searcher(c: &mut Criterion) {
    let mut group = c.benchmark_group("seed_prod_shaped_searcher");
    group.sample_size(10);

    let seed_config =
        SeedConfig::with_wobble(SeedSpec::LengthOnly(7), MismatchSpec::new(1, 2, 2), true);
    let dataset = build_production_dataset(10, 22, 100_000, &seed_config);
    let global = dataset.store.global_view();
    let prepared_queries: Vec<_> = dataset
        .queries
        .entries()
        .iter()
        .filter_map(|query| prepare_query_for_seed_search(query, &seed_config))
        .collect();
    let mut results: Vec<SeedMatch> = Vec::with_capacity(16_384);

    group.bench_function("prepared_queries_10q_x_100k", |b| {
        b.iter(|| {
            for prepared in &prepared_queries {
                results.clear();
                let searcher = SeedSearcher::new(
                    black_box(prepared.padded_q_sa.as_slice()),
                    black_box(prepared.padded_q_seq.as_slice()),
                    black_box(prepared.q_sa_start),
                    black_box(prepared.q_sa_len),
                    black_box(global.combined_sa),
                    black_box(global.combined_seq),
                    black_box(global.sa_real_len),
                    black_box(&seed_config),
                );
                searcher.search_length_range(prepared.min_len, prepared.max_len, &mut results);
                black_box(results.len());
            }
        });
    });

    group.finish();
}

fn bench_seed_prod_shaped_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("seed_prod_shaped_search");
    group.sample_size(10);

    let seed_config =
        SeedConfig::with_wobble(SeedSpec::LengthOnly(7), MismatchSpec::new(1, 2, 2), true);
    let dataset = build_production_dataset(10, 22, 100_000, &seed_config);
    let args = make_search_args(&seed_config);
    let output_path = dataset._tmpdir.path().join("search.out");

    group.bench_function("run_search_10q_x_100k", |b| {
        b.iter(|| {
            let emitted = run_search(
                black_box(&dataset.queries),
                black_box(&dataset.store),
                black_box(&args),
                black_box(output_path.as_path()),
            )
            .expect("run search");
            black_box(emitted);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_seed_exact,
    bench_seed_mismatch,
    bench_seed_prod_shaped_mismatch,
    bench_seed_realistic,
    bench_seed_prod_shaped_searcher,
    bench_seed_prod_shaped_pipeline
);
criterion_main!(benches);
