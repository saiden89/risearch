//! Mismatch parity tests: -m c:p (max mismatches : min consecutive matches).

mod support;

use rstest::rstest;
use std::path::PathBuf;
use support::{ParityRunner, query_fa, target_fa};

#[rstest]
fn mismatch(
    query_fa: PathBuf,
    target_fa: PathBuf,
    #[values("1:0", "1:3", "2:2")] spec: &str,
    #[values(6, 8, 10)] s: usize,
    #[values(0, 10)] l: usize,
) {

    let l_str = l.to_string();
    let s_str = s.to_string();
    let args = [
        "-l",
        &l_str,
        "-e",
        "100.0",
        "--seed-length",
        &s_str,
        "-m",
        spec,
        "-p3",
    ];

    let test_name = format!("{}_s{}_l{}", spec.replace(':', "_"), s, l);
    ParityRunner::new(&target_fa).assert_pass(&query_fa, &test_name, &args);
}
