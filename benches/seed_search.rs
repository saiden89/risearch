//! Rust-only performance benchmarks
//!
//! Run with: cargo bench --bench seed_search
//! Results in target/criterion/report/index.html

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use std::path::PathBuf;
use std::process::Command;
use std::time::Duration;

fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn query() -> PathBuf {
    project_root().join("legacy_c/RIsearch2/test_suite/mirnas.fa")
}

fn target() -> PathBuf {
    project_root().join("chr22.idx")
}

// =============================================================================
// RUST-ONLY SEARCH BENCHMARKS
// =============================================================================

fn run_rust_search(seed_len: usize, max_ext: usize) {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_risearch"));
    cmd.arg("search")
        .arg("-q")
        .arg(query())
        .arg("-i")
        .arg(target())
        .arg("-s")
        .arg(seed_len.to_string())
        .arg("-l")
        .arg(max_ext.to_string())
        .arg("-e")
        .arg("100.0")
        .arg("-o")
        .arg("/dev/null");

    let output = cmd.output().expect("Rust binary failed");
    black_box(output.stdout);
}

fn bench_rust_search(c: &mut Criterion) {
    let mut group = c.benchmark_group("rust_search");
    group.measurement_time(Duration::from_secs(10));

    for seed_len in [6, 8, 10] {
        for max_ext in [0, 10, 20] {
            let id = format!("s{}_l{}", seed_len, max_ext);
            group.bench_with_input(
                BenchmarkId::new("search", &id),
                &(seed_len, max_ext),
                |b, &(s, l)| b.iter(|| run_rust_search(s, l)),
            );
        }
    }

    group.finish();
}

fn bench_rust_extension(c: &mut Criterion) {
    let mut group = c.benchmark_group("rust_extension");
    group.measurement_time(Duration::from_secs(10));

    // Fixed seed, vary extension length
    let seed_len = 8;
    for max_ext in [0, 10, 20, 50] {
        let id = format!("ext_{}", max_ext);
        group.bench_with_input(BenchmarkId::new("extension", &id), &max_ext, |b, &l| {
            b.iter(|| run_rust_search(seed_len, l))
        });
    }

    group.finish();
}

criterion_group!(benches, bench_rust_search, bench_rust_extension,);
criterion_main!(benches);
