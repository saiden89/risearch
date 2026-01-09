//! Performance benchmarks: Rust vs C RIsearch
//!
//! Run with: cargo bench --bench search_bench
//! Results in target/criterion/report/index.html

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use std::path::PathBuf;
use std::process::Command;
use std::time::Duration;

// =============================================================================
// TEST DATA PATHS
// =============================================================================

fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn mirnas_fa() -> PathBuf {
    project_root().join("legacy_c/RIsearch2/test_suite/mirnas.fa")
}

fn rhoc_idx() -> PathBuf {
    project_root().join("legacy_c/RIsearch2/test_suite/RHOC.idx")
}

fn c_binary() -> PathBuf {
    // Try release first, fall back to debug
    let release = project_root().join("legacy_c/RIsearch2/bin/risearch2.x");
    if release.exists() {
        release
    } else {
        project_root().join("legacy_c/RIsearch2/bin/risearch2.dbg.x")
    }
}

// =============================================================================
// RUST SEARCH (via CLI to match C subprocess overhead)
// =============================================================================

fn run_rust_search(seed_len: usize, max_ext: usize, mismatch: Option<&str>) {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_risearch"));
    cmd.arg("search")
        .arg("-q")
        .arg(mirnas_fa())
        .arg("-i")
        .arg(rhoc_idx())
        .arg("-s")
        .arg(seed_len.to_string())
        .arg("-l")
        .arg(max_ext.to_string())
        .arg("-e")
        .arg("100.0");

    if let Some(m) = mismatch {
        cmd.arg("-m").arg(m);
    }

    let output = cmd.output().expect("Rust binary failed");
    black_box(output.stdout);
}

// =============================================================================
// C SEARCH (via subprocess)
// =============================================================================

fn run_c_search(seed_len: usize, max_ext: usize, mismatch: Option<&str>) {
    let mut cmd = Command::new(c_binary());
    cmd.arg("-q")
        .arg(mirnas_fa())
        .arg("-i")
        .arg(rhoc_idx())
        .arg("-s")
        .arg(seed_len.to_string())
        .arg("-l")
        .arg(max_ext.to_string())
        .arg("-e")
        .arg("100.0");

    if let Some(m) = mismatch {
        cmd.arg("-m").arg(m);
    }

    let output = cmd.output().expect("C binary failed");
    black_box(output.stdout);
}

// =============================================================================
// BENCHMARK GROUPS
// =============================================================================

fn bench_exact_match(c: &mut Criterion) {
    let mut group = c.benchmark_group("exact_match");
    group.measurement_time(Duration::from_secs(10));

    for seed_len in [6, 8, 10] {
        for max_ext in [0, 10, 20] {
            let id = format!("s{}_l{}", seed_len, max_ext);

            group.bench_with_input(
                BenchmarkId::new("rust", &id),
                &(seed_len, max_ext),
                |b, &(s, l)| b.iter(|| run_rust_search(s, l, None)),
            );

            group.bench_with_input(
                BenchmarkId::new("c", &id),
                &(seed_len, max_ext),
                |b, &(s, l)| b.iter(|| run_c_search(s, l, None)),
            );
        }
    }

    group.finish();
}

fn bench_mismatch_search(c: &mut Criterion) {
    let mut group = c.benchmark_group("mismatch");
    group.measurement_time(Duration::from_secs(15));
    group.sample_size(20); // Mismatch is slower

    for mismatch in ["1:0", "1:3", "2:2"] {
        for seed_len in [6, 8] {
            let max_ext = 10;
            let id = format!("m{}_s{}_l{}", mismatch.replace(':', "_"), seed_len, max_ext);

            group.bench_with_input(
                BenchmarkId::new("rust", &id),
                &(seed_len, max_ext, mismatch),
                |b, &(s, l, m)| b.iter(|| run_rust_search(s, l, Some(m))),
            );

            group.bench_with_input(
                BenchmarkId::new("c", &id),
                &(seed_len, max_ext, mismatch),
                |b, &(s, l, m)| b.iter(|| run_c_search(s, l, Some(m))),
            );
        }
    }

    group.finish();
}

fn bench_heavy_mismatch(c: &mut Criterion) {
    let mut group = c.benchmark_group("heavy_mismatch");
    group.measurement_time(Duration::from_secs(20));
    group.sample_size(10); // Very slow

    let mismatch = "3:0";
    let seed_len = 6;
    let max_ext = 10;
    let id = format!("m{}_s{}_l{}", mismatch.replace(':', "_"), seed_len, max_ext);

    group.bench_with_input(
        BenchmarkId::new("rust", &id),
        &(seed_len, max_ext, mismatch),
        |b, &(s, l, m)| b.iter(|| run_rust_search(s, l, Some(m))),
    );

    group.bench_with_input(
        BenchmarkId::new("c", &id),
        &(seed_len, max_ext, mismatch),
        |b, &(s, l, m)| b.iter(|| run_c_search(s, l, Some(m))),
    );

    group.finish();
}

// =============================================================================
// BENCHMARK GROUPS
// =============================================================================

criterion_group!(
    benches,
    bench_exact_match,
    bench_mismatch_search,
    bench_heavy_mismatch,
);
criterion_main!(benches);
