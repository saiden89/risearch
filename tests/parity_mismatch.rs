//! Mismatch parity tests: -m c:p (max mismatches : min consecutive matches).

mod support;

use rstest::rstest;
use support::{workspace_root, ParityRunner};

#[rstest]
fn test_mismatch(
    #[values("1:0", "1:3", "2:2")] spec: &str,
    #[values(6, 8, 10)] s: usize,
    #[values(0, 10)] l: usize,
) {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");

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
    ParityRunner::new(&target).assert_pass(&query, &test_name, &args);
}
