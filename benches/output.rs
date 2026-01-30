//! Benchmark output formatting (where get_name is called)
//!
//! Run with: cargo bench --bench output

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use risearch::config::OutputFormat;
use risearch::registry::{QueryRegistry, Registry, SequenceIndex};
use risearch::search::SearchHit;
use risearch::seq::Sequence;
use risearch::types::{Base, Energy, Strand};

// =============================================================================
// TEST DATA GENERATION
// =============================================================================

fn generate_sequence(len: usize, seed: u64) -> Sequence {
    let mut rng = seed;
    let bases: Vec<Base> = (0..len)
        .map(|_| {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            match (rng >> 32) % 4 {
                0 => Base::A,
                1 => Base::G,
                2 => Base::C,
                _ => Base::U,
            }
        })
        .collect();
    Sequence::from(bases)
}

fn setup_registries(num_queries: usize, _num_targets: usize) -> QueryRegistry {
    // Create simple sequence indices for benchmarking
    let query_entries: Vec<_> = (0..num_queries)
        .map(|i| {
            let name = format!("query_{}", i);
            let seq = generate_sequence(100, i as u64);
            SequenceIndex::new(name, seq)
        })
        .collect();
    Registry::new(query_entries)
}

fn generate_hits(num_hits: usize, num_queries: u32, num_targets: u32) -> Vec<SearchHit> {
    (0..num_hits)
        .map(|i| SearchHit {
            query_idx: (i as u32) % num_queries,
            target_idx: (i as u32) % num_targets,
            q_start: 10,
            q_end: 30,
            t_start: 50,
            t_end: 70,
            output_q_start: 10,
            output_q_end: 30,
            output_t_start: 50,
            output_t_end: 70,
            strand: Strand::Forward,
            energy: Energy::from(-500),
            alignment: None,
            flank_5: vec![].into(),
            flank_3: vec![].into(),
        })
        .collect()
}

// =============================================================================
// BENCHMARKS
// =============================================================================

fn bench_output_formatting(c: &mut Criterion) {
    let mut group = c.benchmark_group("output_formatting");

    // Test with different numbers of hits
    for num_hits in [100, 1000, 10000] {
        let registry = setup_registries(100, 100);
        let hits = generate_hits(num_hits, 100, 100);

        group.bench_with_input(
            BenchmarkId::from_parameter(num_hits),
            &num_hits,
            |b, _| {
                b.iter(|| {
                    let mut output = Vec::with_capacity(num_hits * 128);
                    for hit in &hits {
                        hit.write_with_format(
                            &mut output,
                            OutputFormat::BindingSite,
                            &registry,
                            &registry, // Use same registry for both
                        )
                        .ok();
                    }
                    black_box(output);
                });
            },
        );
    }

    group.finish();
}

fn bench_get_name_only(c: &mut Criterion) {
    let mut group = c.benchmark_group("get_name_microbench");

    for num_entries in [10, 100, 1000] {
        let registry = setup_registries(num_entries, 1);

        group.bench_with_input(
            BenchmarkId::from_parameter(num_entries),
            &num_entries,
            |b, &n| {
                b.iter(|| {
                    // Simulate looking up names for many hits
                    let mut sum = 0usize;
                    for idx in 0..1000 {
                        let name = registry.get_name((idx % n) as u32);
                        sum += name.len();
                    }
                    black_box(sum);
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_output_formatting, bench_get_name_only);
criterion_main!(benches);
