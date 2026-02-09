//! Regression tests for seed-length pre-pruning logic in SA traversal.
//!
//! These cases are intentionally tiny so they run fast, but they exercise
//! root-interval behavior where suffix-length validity is not monotonic
//! in suffix-array order.

mod common;

use common::SingleSeqRunner;

#[test]
fn test_root_preprune_non_monotonic_validity() {
    // Query "AC" has SA [0,1]. For min_len=2 (offset=1), suffix validity in SA
    // order is [true,false], which is not partitioned for binary-search-based
    // validity checks. Rust and C must still agree on hits.
    let query = "AC";
    let target = "AC";
    let args = ["-l", "20", "-e", "10000", "--seed-length", "2", "-p3"];
    SingleSeqRunner::new(query, target).assert_pass("root_preprune_non_monotonic", &args);
}
