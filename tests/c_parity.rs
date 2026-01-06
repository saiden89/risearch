//! C Parity Tests
//!
//! These tests compare Rust risearch output against the legacy C risearch2 implementation.
//! Run with RUST_LOG=debug for detailed diff analysis, or RUST_LOG=trace for raw data.

mod common;

use common::{run_file_parity, run_single_seq_parity, setup_common_test_files, workspace_root};

// =============================================================================
// FILE-BASED TESTS
// =============================================================================

#[test]
fn test_parity_full_pipeline() {
    let (_, query_path, target_path, _) = setup_common_test_files();
    let args = ["-l", "20", "-e", "-20", "-s", "6", "-p3"];
    run_file_parity(&query_path, &target_path, None, "default_config", &args);
}

#[test]
fn test_parity_long_seed() {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");
    let c_index = root.join("legacy_c/RIsearch2/test_suite/RHOC.pksuf");
    let args = ["-l", "0", "-e", "10000", "-s", "12", "-p3"];
    run_file_parity(&query, &target, Some(&c_index), "long_seed_no_ext", &args);
}

#[test]
fn test_parity_energy_threshold() {
    let (_, query_path, target_path, _) = setup_common_test_files();
    let args = ["-l", "10", "-e", "-10.0", "-s", "5", "-p3"];
    run_file_parity(&query_path, &target_path, None, "energy_only", &args);
}

#[test]
fn test_parity_alignment_repro() {
    let query = "uggcucaguucagcaggaacag";
    let target = "TGGCTCTGTGGGACACAGCAGG";
    let args = ["-l", "20", "-e", "-20", "-s", "6", "-p3"];
    run_single_seq_parity(query, target, "alignment_mismatch_repro", &args, false);
}

#[test]
fn test_parity_single_seq() {
    // Tests internal mismatch handling between Rust and C implementations.
    let args = ["-l", "20", "-e", "100.0", "-s", "5", "-p3"];

    let query = "UGCUGCUGCCGCUGCUGCUG"; // 20nt with internal variation
    let target = "GCAGCAGCAGCAGCAGCAGC"; // 20nt complement

    run_single_seq_parity(query, target, "internal_mismatch_20nt", &args, false);
}

// =============================================================================
// ISOLATED TESTS: Each tests a specific component to pinpoint differences
// =============================================================================

/// Tests seed matching only (no extension).
/// Use -l 0 to disable extension, so only seed pairing is tested.
#[test]
fn test_parity_seed_only() {
    // No extension: -l 0
    // This tests only the seed pairing and energy calculation
    let args = ["-l", "0", "-e", "100.0", "-s", "5", "-p3"];

    // Query and target are exact complements (20nt)
    // Should produce a single seed hit with no extensions
    let query = "UGCUGCUGCUGCUGCUGCUG"; // 20nt
    let target = "CAGCAGCAGCAGCAGCAGCA"; // Perfect complement, reversed

    run_single_seq_parity(query, target, "seed_only_no_extension", &args, false);
}

/// Tests left extension only (dp_left).
/// Design: seed at 3' end of query, extra bases only to the 5' side.
#[test]
fn test_parity_left_ext() {
    // Seed at 3' end of query forces only left extension
    let args = ["-l", "20", "-e", "100.0", "-s", "5", "-p3"];

    let query = "AAAAAUGCUG";
    let target = "CAGCAUUUUU";

    run_single_seq_parity(query, target, "left_extension_only", &args, false);
}

/// Tests right extension only (dp_right).
/// Design: seed at 5' end of query, extra bases only to the 3' side.
#[test]
fn test_parity_right_ext() {
    let args = ["-l", "20", "-e", "100.0", "-s", "5", "-p3"];

    let query = "UGCUGAAAAA";
    let target = "UUUUUCAGCA";

    run_single_seq_parity(query, target, "right_extension_only", &args, false);
}

/// Tests both left and right extension.
/// Design: seed in middle, extra bases on both sides.
#[test]
fn test_parity_both_ext() {
    let args = ["-l", "20", "-e", "100.0", "-s", "5", "-p3"];

    let query = "AAAUGCUGAAA";
    let target = "UUUCAGCAUUU";

    run_single_seq_parity(query, target, "both_extensions", &args, false);
}

/// Tests with wobble pairs (G-U) in the seed region.
#[test]
fn test_parity_wobble() {
    let args = ["-l", "0", "-e", "100.0", "-s", "5", "-p3"];

    let query = "UGUGUGUGUG"; // 10 alternating U-G pattern
    let target = "CGCGCGCGCG"; // Complement with wobble

    run_single_seq_parity(query, target, "wobble_seed_pairs", &args, false);
}

// =============================================================================
// ISOLATED FAILURE REPRODUCTION: Uses actual failing sequences from parity_default_config
// =============================================================================

/// Isolated repro of parity_default_config failure.
/// Uses actual hsa-miR-24-3p query and a short RHOC target region.
#[test]
fn test_parity_mir24_isolated() {
    let query = "uggcucaguucagcaggaacag"; // 22nt

    let target =
        "GGAAGACCTGCCTCCTCATCGTCTTCAGCAAGGATCAGTTTCCGGAGGTCTACGTCCCTACTGTCTTTGAGAACTATATTG";

    let args = [
        "-l",
        "20",
        "-e",
        "100.0",
        "-s",
        "6",
        "-p3",
        "--no-max-prune",
    ];

    run_single_seq_parity(query, target, "miR24_isolated_single", &args, true);
}

/// More targeted single-hit test with explicit segment expectations.
#[test]
fn test_parity_segments_debug() {
    let query = "AAAAAUGCUGUAAAAA"; // 5 + 6 + 5 = 16nt
    let target = "UUUUUACAGCAUUUUU";
    let args = ["-l", "20", "-e", "100.0", "-s", "6", "-p3"];

    run_single_seq_parity(query, target, "segment_debug", &args, false);
}

/// Isolated test for rust-worse energy case.
/// From full_pipeline: q=[1,22] t=[359,384] with E=-21.74 (C: -22.38)
/// The FP diff shows gap placement difference in 3' EXT:
///   C:    UQPPP (gap-then-pair)
///   Rust: UUQPP (unpaired-then-gap)
#[test]
fn test_parity_rust_worse_debug() {
    // Query: hsa-miR-24-3p
    let query = "uggcucaguucagcaggaacag"; // 22nt

    // Target: ENST00000414971 positions 320-420 (contains t=[359,384])
    // Extracted directly from RHOC.fa
    let target = "ATATTGCGGACATTGAGGTGGACGGCAAGCAGGTGGAGCTGGCTCTGTGGGACACAGCAGGGCAGGAAGACTATGATCGACTGCGGCCTCTCTCCTACCCG";

    let args = ["-l", "20", "-e", "100.0", "-s", "6", "-p3"];

    run_single_seq_parity(query, target, "rust_worse_debug", &args, false);
}

/// Minimal reproducible test for identical-alignment-different-energy issue.
/// q[1,17] t[20,33]: Rust=-15.04, C=-15.35 (31 centiunit difference)
/// Both produce identical alignment: Seed=PPWPPPPP, 3'EXT=UQQQPPPWP
/// but the energy scores differ, indicating DSM lookup or initialization mismatch.
#[test]
fn test_parity_energy_discrepancy_minimal() {
    // Query: hsa-miR-24-3p (positions 1-17 used in this hit)
    let query = "uggcucaguucagcaggaacag"; // 22nt

    // Target: ENST00000534717 (positions 20-33 used in this hit)
    // Only first 60nt needed to reproduce the issue
    let target = "AAGCCCGGGAAGCTGACTCCTTGCCCTGAGTCACAGGGAGGGGUGGGCAGGGCATGCGGC";

    // Args matching energy_threshold test: -l 10 -e -10 -s 5 -p3
    let args = ["-l", "10", "-e", "-10.0", "-s", "5", "-p3"];

    run_single_seq_parity(query, target, "energy_discrepancy_minimal", &args, false);
}
