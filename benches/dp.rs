use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use risearch::config::Matrix;
use risearch::dp::gotoh::Gotoh;
use risearch::dp::{DpGrid, DpView, ExtendDir};
use risearch::dsm::ScoringModel;
use risearch::seq::Sequence;
use risearch::types::Base;

// ============================================================================
// SIMPLE LINEAR CONGRUENTIAL GENERATOR (LCG)
// ============================================================================
// A deterministic PRNG without external dependencies. Provides good coverage
// of RNA bases for realistic DP workloads.

struct SimpleLcg {
    state: u64,
}

impl SimpleLcg {
    const A: u64 = 6364136223846793005;
    const C: u64 = 1442695040888963407;

    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Generate next pseudo-random u64
    fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(Self::A).wrapping_add(Self::C);
        self.state
    }

    /// Generate a random RNA base (A, G, C, U, or N).
    fn next_base(&mut self) -> Base {
        match (self.next() >> 32) % 5 {
            0 => Base::A,
            1 => Base::G,
            2 => Base::C,
            3 => Base::U,
            _ => Base::N,
        }
    }
}

// ============================================================================
// SEQUENCE GENERATION
// ============================================================================

/// Generate a pseudo-random RNA sequence of given length.
fn generate_sequence(len: usize, seed: u64) -> Sequence {
    let mut rng = SimpleLcg::new(seed);
    let bases: Vec<Base> = (0..len).map(|_| rng.next_base()).collect();
    Sequence::from(bases)
}

// ============================================================================
// BENCHMARK GROUPS
// ============================================================================

fn bench_extend_left(c: &mut Criterion) {
    let mut group = c.benchmark_group("extend_left");
    let model = ScoringModel::new(Matrix::T04, 0);
    let left_model = model.transpose();
    let gotoh = Gotoh::new(&left_model);

    for len in [10, 20, 30, 50].iter() {
        group.bench_with_input(BenchmarkId::from_parameter(len), len, |b, &len| {
            let query = generate_sequence(100, 12345);
            let target = generate_sequence(100, 54321);
            let q_start = 50;
            let t_start = 50;

            let mut grid = DpGrid::new(200);

            b.iter(|| {
                let view = DpView::new(
                    black_box(&query),
                    black_box(&target),
                    black_box(q_start),
                    black_box(t_start),
                    ExtendDir::Left,
                    black_box(len),
                );
                let _ = gotoh.extend(black_box(&view), &mut grid);
            });
        });
    }

    group.finish();
}

fn bench_extend_right(c: &mut Criterion) {
    let mut group = c.benchmark_group("extend_right");
    let model = ScoringModel::new(Matrix::T04, 0);
    let gotoh = Gotoh::new(&model);

    for len in [10, 20, 30, 50].iter() {
        group.bench_with_input(BenchmarkId::from_parameter(len), len, |b, &len| {
            let query = generate_sequence(100, 12345);
            let target = generate_sequence(100, 54321);
            let q_end = 49;
            let t_end = 49;

            let mut grid = DpGrid::new(200);

            b.iter(|| {
                let view = DpView::new(
                    black_box(&query),
                    black_box(&target),
                    black_box(q_end),
                    black_box(t_end),
                    ExtendDir::Right,
                    black_box(len),
                );
                let _ = gotoh.extend(black_box(&view), &mut grid);
            });
        });
    }

    group.finish();
}

/// Throughput benchmark: measure DP cells computed per second.
fn bench_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("throughput");
    group.sample_size(10); // Smaller sample size for realistic wall-clock time
    let model = ScoringModel::new(Matrix::T04, 0);
    let left_model = model.transpose();
    let gotoh_right = Gotoh::new(&model);
    let gotoh_left = Gotoh::new(&left_model);

    for len in [10, 20, 30, 50].iter() {
        // Left extension throughput: query_len × target_len DP cells
        group.bench_with_input(
            BenchmarkId::new("extend_left_cells", len),
            len,
            |b, &len| {
                let query = generate_sequence(100, 12345);
                let target = generate_sequence(100, 54321);
                let q_start = 50;
                let t_start = 50;

                let mut grid = DpGrid::new(200);

                b.iter_custom(|iters| {
                    let mut total_duration = std::time::Duration::ZERO;

                    for _ in 0..iters {
                        let start = std::time::Instant::now();
                        let view = DpView::new(
                            black_box(&query),
                            black_box(&target),
                            black_box(q_start),
                            black_box(t_start),
                            ExtendDir::Left,
                            black_box(len),
                        );
                        let result = gotoh_left.extend(black_box(&view), &mut grid);
                        total_duration += start.elapsed();

                        // Ensure result is not optimized away
                        black_box(result);
                    }

                    total_duration
                });
            },
        );

        // Right extension throughput: query_len × target_len DP cells
        group.bench_with_input(
            BenchmarkId::new("extend_right_cells", len),
            len,
            |b, &len| {
                let query = generate_sequence(100, 12345);
                let target = generate_sequence(100, 54321);
                let q_end = 49;
                let t_end = 49;

                let mut grid = DpGrid::new(200);

                b.iter_custom(|iters| {
                    let mut total_duration = std::time::Duration::ZERO;

                    for _ in 0..iters {
                        let start = std::time::Instant::now();
                        let view = DpView::new(
                            black_box(&query),
                            black_box(&target),
                            black_box(q_end),
                            black_box(t_end),
                            ExtendDir::Right,
                            black_box(len),
                        );
                        let result = gotoh_right.extend(black_box(&view), &mut grid);
                        total_duration += start.elapsed();

                        // Ensure result is not optimized away
                        black_box(result);
                    }

                    total_duration
                });
            },
        );
    }

    group.finish();
}

/// Many extensions benchmark: simulates a search workload with 100 extensions.
/// Each extension pair (left + right) is a typical search seed extension pattern.
fn bench_many_extensions(c: &mut Criterion) {
    let mut group = c.benchmark_group("many_extensions");
    group.sample_size(10);
    let model = ScoringModel::new(Matrix::T04, 0);
    let left_model = model.transpose();
    let gotoh_right = Gotoh::new(&model);
    let gotoh_left = Gotoh::new(&left_model);

    group.bench_function("100_extensions_len30", |b| {
        let queries: Vec<Sequence> = (0..10)
            .map(|i| generate_sequence(100, 1000 + i as u64))
            .collect();
        let targets: Vec<Sequence> = (0..10)
            .map(|i| generate_sequence(100, 5000 + i as u64))
            .collect();

        let mut grid = DpGrid::new(200);

        b.iter(|| {
            let mut total_score = 0i32;

            // Simulate 100 search hit extensions (10 query × 10 target combinations)
            for (q_idx, query) in queries.iter().enumerate() {
                for (t_idx, target) in targets.iter().enumerate() {
                    // Anchor positions based on indices for variety
                    let q_start = (50 + q_idx) % query.len();
                    let t_start = (50 + t_idx) % target.len();
                    let q_end = if q_start + 30 < query.len() {
                        q_start + 30
                    } else {
                        query.len() - 1
                    };
                    let t_end = if t_start + 30 < target.len() {
                        t_start + 30
                    } else {
                        target.len() - 1
                    };

                    // Left extension
                    let left_view = DpView::new(
                        black_box(query),
                        black_box(target),
                        black_box(q_start),
                        black_box(t_start),
                        ExtendDir::Left,
                        black_box(30),
                    );
                    let left_result = gotoh_left.extend(black_box(&left_view), &mut grid);
                    total_score = total_score.wrapping_add(left_result.energy);

                    // Right extension
                    let right_view = DpView::new(
                        black_box(query),
                        black_box(target),
                        black_box(q_end),
                        black_box(t_end),
                        ExtendDir::Right,
                        black_box(30),
                    );
                    let right_result = gotoh_right.extend(black_box(&right_view), &mut grid);
                    total_score = total_score.wrapping_add(right_result.energy);
                }
            }

            black_box(total_score);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_extend_left,
    bench_extend_right,
    bench_throughput,
    bench_many_extensions
);
criterion_main!(benches);
