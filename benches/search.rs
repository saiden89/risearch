use std::fs;
use std::path::Path;

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use risearch::config::{
    ExtendConfig, FilterConfig, OutputCompression, OutputConfig, OutputFormat,
    ScoreConfig, SearchConfig, SeedConfig,
};
use risearch::registry::QueryRegistry;
use risearch::search::{run_search, run_search_in_memory};
use risearch::Energy;
use risearch::seed::{SeedHit, SeedingEngine};
use risearch::seq::Sequence;
use risearch::types::{Base, DsmId};
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

struct ProductionSearchDataset {
    tmpdir: TempDir,
    queries: QueryRegistry,
    store: TargetStore,
}

fn generate_sequence(len: usize, seed: u64) -> Sequence {
    let mut rng = SimpleLcg::new(seed);
    let bases: Vec<Base> = (0..len).map(|_| rng.next_base()).collect();
    Sequence::from(bases)
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
) -> ProductionSearchDataset {
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
    TargetStore::build(&targets_path, &index_path).expect("build target index");
    let store = TargetStore::open(&index_path).expect("open target index");

    ProductionSearchDataset {
        tmpdir,
        queries,
        store,
    }
}

fn make_search_config(seed_config: &SeedConfig) -> SearchConfig {
    SearchConfig {
        seed: seed_config.clone(),
        score: ScoreConfig {
            dsm_id: DsmId::from("t04"),
            penalty: Energy::from_kcal(0.0),
            temperature: 37,
        },
        extend: ExtendConfig {
            max_extension: 20,
        },
        filter: FilterConfig {
            delta_g: Energy::from_kcal(f64::NEG_INFINITY),
            seed_energy: Energy::from_kcal(0.0),
            no_max_prune: false,
        },
        output: OutputConfig {
            format: OutputFormat::Minimal,
            compress: OutputCompression::None,
            multifile: false,
        },
    }
}

fn seed_count(groups: &[(usize, Vec<SeedHit>)]) -> usize {
    groups.iter().map(|(_, seeds)| seeds.len()).sum()
}

fn bench_search_prod_shaped_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("prod_shaped");
    group.sample_size(10);

    let seed_config = SeedConfig {
        seed_start: None,
        seed_end: None,
        seed_length: Some(7),
        seed_wobble: true,
        max_mismatches: 1,
        min_prefix_matches: 2,
        min_suffix_matches: 2,
    };
    let dataset = build_production_dataset(10, 22, 100_000, &seed_config);
    let args = make_search_config(&seed_config);
    let output_path = dataset.tmpdir.path().join("search.out");

    // Bench the seed stage alone on the same production-shaped dataset used
    // by the full search pipeline benches below.
    let case = "10q_x_100k";

    group.bench_with_input(BenchmarkId::new("collect", case), &case, |b, _| {
        b.iter(|| {
            let seeds = SeedingEngine::new(
                black_box(&dataset.queries),
                black_box(&dataset.store),
            ).run(black_box(&seed_config));
            black_box(seed_count(&seeds));
        });
    });

    // Bench the end-to-end search path when hits are collected in memory.
    group.bench_with_input(
        BenchmarkId::new("run_search_in_memory", case),
        &case,
        |b, _| {
            b.iter(|| {
                let hits = run_search_in_memory(
                    black_box(&dataset.queries),
                    black_box(&dataset.store),
                    black_box(&args),
                )
                .expect("run search in memory");
                black_box(hits.len());
            });
        },
    );

    // Bench the end-to-end search path including output writing.
    group.bench_with_input(BenchmarkId::new("run_search", case), &case, |b, _| {
        b.iter(|| {
            run_search(
                black_box(&dataset.queries),
                black_box(&dataset.store),
                black_box(&args),
                black_box(output_path.as_path()),
            )
            .expect("run search");
        });
    });

    group.finish();
}

criterion_group!(benches, bench_search_prod_shaped_pipeline);
criterion_main!(benches);
