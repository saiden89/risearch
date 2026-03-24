//! Wobble parity tests: G-U wobble pairing in seed region.

mod support;

use support::SingleSeqRunner;

#[test]
fn wobble_seed() {
    let args = ["-l", "0", "-e", "100.0", "--seed-length", "5", "-p3"];
    let query = "UGUGUGUGUG"; // alternating U-G
    let target = "CGCGCGCGCG"; // complement with wobble
    SingleSeqRunner::new(query, target).assert_pass("wobble_seed", &args);
}
