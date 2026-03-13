//! Extension parity tests: DP extension logic (left, right, both).

mod support;

use support::SingleSeqRunner;

#[test]
fn test_seed_only() {
    // No extension: -l 0
    let args = ["-l", "0", "-e", "100.0", "--seed-length", "5", "-p3"];
    let query = "UGCUGCUGCUGCUGCUGCUG";
    let target = "CAGCAGCAGCAGCAGCAGCA";
    SingleSeqRunner::new(query, target).assert_pass("seed_only", &args);
}

#[test]
fn test_left_only() {
    // Seed at 3' end forces left extension only
    let args = ["-l", "20", "-e", "100.0", "--seed-length", "5", "-p3"];
    let query = "AAAAAUGCUG";
    let target = "CAGCAUUUUU";
    SingleSeqRunner::new(query, target).assert_pass("left_only", &args);
}

#[test]
fn test_right_only() {
    // Seed at 5' end forces right extension only
    let args = ["-l", "20", "-e", "10.0", "--seed-length", "5", "-p3"];
    let query = "UGCUGAAAAA";
    let target = "UUUUUCAGCA";
    SingleSeqRunner::new(query, target).assert_pass("right_only", &args);
}

#[test]
fn test_both_sides() {
    // Seed in middle, extension both sides
    let args = ["-l", "20", "-e", "100.0", "--seed-length", "5", "-p3"];
    let query = "AAAUGCUGAAA";
    let target = "UUUCAGCAUUU";
    SingleSeqRunner::new(query, target).assert_pass("both_sides", &args);
}

#[test]
fn test_internal_mismatch() {
    let args = ["-l", "20", "-e", "100.0", "--seed-length", "5", "-p3"];
    let query = "UGCUGCUGCCGCUGCUGCUG";
    let target = "GCAGCAGCAGCAGCAGCAGC";
    SingleSeqRunner::new(query, target).assert_pass("internal_mismatch", &args);
}
