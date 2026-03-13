//! Seed parity tests: seed length, interval, and spec configurations.

mod support;

use rstest::rstest;
use support::{workspace_root, ParityRunner};

/// Main seed length × extension limit matrix.
/// Tests all combinations on mirnas.fa vs RHOC.fa.
#[rstest]
fn test_length(
    #[values(6, 7, 8, 9, 10, 11, 12, 13, 14, 15)] s: usize,
    #[values(0, 5, 10, 15, 20, 25, 30, 35, 40)] l: usize,
    #[values(false, true)] strict: bool,
) {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");

    let l_str = l.to_string();
    let s_str = s.to_string();
    let mut args = vec!["-l", &l_str, "-e", "100.0", "--seed-length", &s_str, "-p3"];
    if strict {
        args.extend(["--seed-pairing", "strict"]);
    }

    let test_name = if strict {
        format!("s{}_l{}_strict", s, l)
    } else {
        format!("s{}_l{}", s, l)
    };
    ParityRunner::new(&target).assert_pass(&query, &test_name, &args);
}

/// Seed interval: --seed-start/--seed-end
#[rstest]
fn test_interval(
    #[values("1:8", "1:12", "2:10", "1:15")] spec: &str,
    #[values(0, 10, 20)] l: usize,
) {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");

    let l_str = l.to_string();
    let (start, end) = spec.split_once(':').unwrap();
    let args = [
        "-l",
        &l_str,
        "-e",
        "100.0",
        "--seed-start",
        start,
        "--seed-end",
        end,
        "-p3",
    ];

    let test_name = format!("interval_{}_l{}", spec.replace(':', "_"), l);
    ParityRunner::new(&target).assert_pass(&query, &test_name, &args);
}

/// Seed interval with min length: --seed-start/--seed-end/--seed-length
#[rstest]
fn test_interval_with_length(
    #[values("1:12/6", "2:10/5", "1:15/7")] spec: &str,
    #[values(0, 10, 20)] l: usize,
) {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");

    let l_str = l.to_string();
    let (start, rest) = spec.split_once(':').unwrap();
    let (end, len) = rest.split_once('/').unwrap();
    let args = [
        "-l",
        &l_str,
        "-e",
        "100.0",
        "--seed-start",
        start,
        "--seed-end",
        end,
        "--seed-length",
        len,
        "-p3",
    ];

    let test_name = format!("interval_{}_l{}", spec.replace([':', '/'], "_"), l);
    ParityRunner::new(&target).assert_pass(&query, &test_name, &args);
}
