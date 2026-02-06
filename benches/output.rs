//! Benchmark registry name lookups (where get_name optimization applies)
//!
//! Run with: cargo bench --bench output

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use risearch::registry::Registry;
use risearch::seq::Sequence;
use risearch::types::Base;

// =============================================================================
// TEST DATA GENERATION
// =============================================================================

fn generate_sequence(len: usize, seed: u64) -> Sequence {
    let mut rng = seed;
    let bases: Vec<Base> = (0..len)
        .map(|_| {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
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

// Simple struct that implements RegistryEntry for testing
struct TestEntry {
    name: String,
    _sequence: Sequence,
}

impl TestEntry {
    fn new(name: String, sequence: Sequence) -> Self {
        Self {
            name,
            _sequence: sequence,
        }
    }
}

impl risearch::registry::RegistryEntry for TestEntry {
    fn name(&self) -> &str {
        &self.name
    }
}

fn setup_registry(num_entries: usize) -> Registry<TestEntry> {
    let entries: Vec<_> = (0..num_entries)
        .map(|i| {
            let name = format!("entry_{:06}", i);
            let seq = generate_sequence(100, i as u64);
            TestEntry::new(name, seq)
        })
        .collect();
    Registry::new(entries)
}

// =============================================================================
// BENCHMARKS
// =============================================================================

/// Microbenchmark: Just the get_name call in a tight loop
fn bench_get_name_microbench(c: &mut Criterion) {
    let mut group = c.benchmark_group("get_name_microbench");

    for num_entries in [10, 100, 1000] {
        let registry = setup_registry(num_entries);

        group.bench_with_input(
            BenchmarkId::from_parameter(num_entries),
            &num_entries,
            |b, &n| {
                b.iter(|| {
                    // Simulate looking up names for many hits
                    let mut sum = 0usize;
                    for idx in 0..10000 {
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

/// Realistic benchmark: Simulates the output path where get_name is called
fn bench_simulated_output(c: &mut Criterion) {
    let mut group = c.benchmark_group("simulated_output");

    for num_hits in [100, 1000, 10000] {
        let query_registry = setup_registry(100);
        let target_registry = setup_registry(100);

        // Simulate hit indices
        let hit_indices: Vec<(u32, u32)> = (0..num_hits)
            .map(|i| ((i % 100) as u32, ((i * 7) % 100) as u32))
            .collect();

        group.bench_with_input(BenchmarkId::from_parameter(num_hits), &num_hits, |b, _| {
            b.iter(|| {
                // Simulate the output formatting path
                let mut output = Vec::with_capacity(num_hits * 64);
                for (q_idx, t_idx) in &hit_indices {
                    let q_name = query_registry.get_name(*q_idx);
                    let t_name = target_registry.get_name(*t_idx);
                    // Simulate formatting (concatenate names)
                    output.extend_from_slice(q_name.as_bytes());
                    output.push(b'\t');
                    output.extend_from_slice(t_name.as_bytes());
                    output.push(b'\n');
                }
                black_box(output);
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_get_name_microbench, bench_simulated_output);
criterion_main!(benches);
