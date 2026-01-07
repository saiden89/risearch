//! C Parity Tests
//!
//! These tests compare Rust risearch output against the legacy C risearch2 implementation.
//! Run with RUST_LOG=debug for detailed diff analysis, or RUST_LOG=trace for raw data.

mod common;

use common::{ParityRunner, SingleSeqRunner, workspace_root};
use rstest::rstest;

// =============================================================================
// FILE-BASED TESTS
// =============================================================================

#[test]
fn test_parity_full_pipeline() {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");
    let args = ["-l", "20", "-e", "-20", "-s", "6", "-p3"];

    ParityRunner::new(&target).assert_pass(&query, "default_config", &args);
}

#[test]
fn test_parity_long_seed() {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");
    let args = ["-l", "0", "-e", "10000", "-s", "12", "-p3"];

    ParityRunner::new(&target).assert_pass(&query, "long_seed_no_ext", &args);
}

#[test]
fn test_parity_energy_threshold() {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");
    let args = ["-l", "10", "-e", "-10.0", "-s", "5", "-p3"];

    ParityRunner::new(&target).assert_pass(&query, "energy_only", &args);
}

#[test]
fn test_parity_alignment_repro() {
    let query = "uggcucaguucagcaggaacag";
    let target = "TGGCTCTGTGGGACACAGCAGG";
    let args = ["-l", "20", "-e", "-20", "-s", "6", "-p3"];

    SingleSeqRunner::new(query, target).assert_pass("alignment_mismatch_repro", &args);
}

#[test]
fn test_parity_single_seq() {
    // Tests internal mismatch handling between Rust and C implementations.
    let args = ["-l", "20", "-e", "100.0", "-s", "5", "-p3"];

    let query = "UGCUGCUGCCGCUGCUGCUG"; // 20nt with internal variation
    let target = "GCAGCAGCAGCAGCAGCAGC"; // 20nt complement

    SingleSeqRunner::new(query, target).assert_pass("internal_mismatch_20nt", &args);
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

    SingleSeqRunner::new(query, target).assert_pass("seed_only_no_extension", &args);
}

/// Minimal test: 2bp seed = exactly ONE stack energy lookup.
/// CG paired with GC gives the strongest stack (-3.30 kcal/mol in Turner 04).
/// Expected energy: (-330 - 559) / -100 = 8.89 kcal/mol (wait, that's positive)
/// Actually: just seed_energy / -100 without extension = -3.30 kcal/mol
/// But with the -559 offset and no extension: (seed_energy + 0 + 0 - 559) / -100
/// For CG/GC: DSM value should be around -330 (3.30 kcal/mol stabilizing).
#[test]
fn test_parity_single_stack() {
    // 2bp seed, no extension
    let args = ["-l", "0", "-e", "10000.0", "-s", "6", "-p3"];

    // Query: CG (5'->3')
    // Target needs to be the complement in reverse: GC reading 3'->5' = CG reading 5'->3'
    // Wait, for antiparallel: Query 5'-CG-3' pairs with Target 3'-GC-5'
    // Target stored 5'->3' = CG, reverse complement to pair = GC
    let query = "GC";
    let target = "GC"; // This gives target bases GC when read for pairing

    SingleSeqRunner::new(query, target).assert_pass("single_stack_cg", &args);
}

/// Tests left extension only (dp_left).
/// Design: seed at 3' end of query, extra bases only to the 5' side.
#[test]
fn test_parity_left_ext() {
    // Seed at 3' end of query forces only left extension
    let args = ["-l", "20", "-e", "100.0", "-s", "5", "-p3"];

    let query = "AAAAAUGCUG";
    let target = "CAGCAUUUUU";

    SingleSeqRunner::new(query, target).assert_pass("left_extension_only", &args);
}

/// Tests right extension only (dp_right).
/// Design: seed at 5' end of query, extra bases only to the 3' side.
#[test]
fn test_parity_right_ext() {
    let args = ["-l", "20", "-e", "10.0", "-s", "5", "-p3"];

    let query = "UGCUGAAAAA";
    let target = "UUUUUCAGCA";

    SingleSeqRunner::new(query, target).assert_pass("right_extension_only", &args);
}

/// Tests both left and right extension.
/// Design: seed in middle, extra bases on both sides.
#[test]
fn test_parity_both_ext() {
    let args = ["-l", "20", "-e", "100.0", "-s", "5", "-p3"];

    let query = "AAAUGCUGAAA";
    let target = "UUUCAGCAUUU";

    SingleSeqRunner::new(query, target).assert_pass("both_extensions", &args);
}

/// Tests with wobble pairs (G-U) in the seed region.
#[test]
fn test_parity_wobble() {
    let args = ["-l", "0", "-e", "100.0", "-s", "5", "-p3"];

    let query = "UGUGUGUGUG"; // 10 alternating U-G pattern
    let target = "CGCGCGCGCG"; // Complement with wobble

    SingleSeqRunner::new(query, target).assert_pass("wobble_seed_pairs", &args);
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

    SingleSeqRunner::new(query, target).assert_pass_detailed("miR24_isolated_single", &args);
}

/// More targeted single-hit test with explicit segment expectations.
#[test]
fn test_parity_segments_debug() {
    let query = "AAAAAUGCUGUAAAAA"; // 5 + 6 + 5 = 16nt
    let target = "UUUUUACAGCAUUUUU";
    let args = ["-l", "20", "-e", "100.0", "-s", "6", "-p3"];

    SingleSeqRunner::new(query, target).assert_pass("segment_debug", &args);
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

    SingleSeqRunner::new(query, target).assert_pass("rust_worse_debug", &args);
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

    SingleSeqRunner::new(query, target).assert_pass("energy_discrepancy_minimal", &args);
}

// =============================================================================
// DINUCLEOTIDE STACK PARITY: Tests all 16 possible dinucleotide stack energies
// =============================================================================

/// Watson-Crick complement for RNA (returns reverse complement for antiparallel pairing).
fn wc_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'U',
            'U' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => c,
        })
        .collect()
}

/// Parameterized test for all 16 dinucleotide stack combinations.
///
/// Each (b1, b2) pair generates a 2bp query and its Watson-Crick complement target,
/// testing exactly one stacking energy lookup for parity between Rust and C.
#[rstest]
fn test_parity_dinucleotide_stack(
    #[values('A', 'C', 'G', 'U')] b1: char,
    #[values('A', 'C', 'G', 'U')] b2: char,
) {
    // 2bp seed, no extension, lenient energy threshold to catch all hits
    let args = ["-l", "0", "-e", "10000.0", "-s", "2", "-p3"];

    let query = format!("{}{}", b1, b2);
    let target = wc_complement(&query);

    SingleSeqRunner::new(&query, &target).assert_pass(&format!("dinuc_stack_{}_{}", b1, b2), &args);
}

/// Parameterized test for all 64 trinucleotide combinations.
///
/// Each (b1, b2, b3) triple generates a 3bp query testing TWO stacking
/// interactions for parity between Rust and C.
#[rstest]
fn test_parity_trinucleotide_stack(
    #[values('A', 'C', 'G', 'U')] b1: char,
    #[values('A', 'C', 'G', 'U')] b2: char,
    #[values('A', 'C', 'G', 'U')] b3: char,
) {
    // 3bp seed, no extension
    let args = ["-l", "0", "-e", "10000.0", "-s", "3", "-p3"];

    let query = format!("{}{}{}", b1, b2, b3);
    let target = wc_complement(&query);

    SingleSeqRunner::new(&query, &target)
        .assert_pass(&format!("trinuc_stack_{}_{}_{}", b1, b2, b3), &args);
}

/// Parameterized test for minimal right extension: 2bp seed + 1bp extension.
///
/// Query is 3bp: first 2bp form the seed, last 1bp triggers dp_right extension.
/// This tests the simplest possible DP extension path.
#[rstest]
fn test_parity_minimal_ext_right(
    #[values('A', 'C', 'G', 'U')] s1: char, // seed base 1
    #[values('A', 'C', 'G', 'U')] s2: char, // seed base 2
    #[values('A', 'C', 'G', 'U')] e1: char, // extension base
) {
    // 2bp seed, with extension enabled
    let args = ["-l", "10", "-e", "10000.0", "-s", "2", "-p3"];

    let query = format!("{}{}{}", s1, s2, e1); // seed at 5', ext at 3'
    let target = wc_complement(&query);

    SingleSeqRunner::new(&query, &target)
        .assert_pass(&format!("min_ext_r_{}_{}_{}", s1, s2, e1), &args);
}

/// Parameterized test for minimal left extension: 1bp extension + 2bp seed.
///
/// Query is 3bp: first 1bp triggers dp_left extension, last 2bp form the seed.
#[rstest]
fn test_parity_minimal_ext_left(
    #[values('A', 'C', 'G', 'U')] e1: char, // extension base
    #[values('A', 'C', 'G', 'U')] s1: char, // seed base 1
    #[values('A', 'C', 'G', 'U')] s2: char, // seed base 2
) {
    // 2bp seed, with extension enabled
    let args = ["-l", "10", "-e", "10000.0", "-s", "2", "-p3"];

    let query = format!("{}{}{}", e1, s1, s2); // ext at 5', seed at 3'
    let target = wc_complement(&query);

    SingleSeqRunner::new(&query, &target)
        .assert_pass(&format!("min_ext_l_{}_{}_{}", e1, s1, s2), &args);
}

/// Parameterized test for extension on both sides: 1bp left + 2bp seed + 1bp right.
///
/// Query is 4bp: first 1bp triggers dp_left, middle 2bp is seed, last 1bp triggers dp_right.
/// This tests both DP extension paths simultaneously.
#[rstest]
fn test_parity_ext_both_sides(
    #[values('A', 'C', 'G', 'U')] el: char, // left extension
    #[values('A', 'C', 'G', 'U')] s1: char, // seed base 1
    #[values('A', 'C', 'G', 'U')] s2: char, // seed base 2
    #[values('A', 'C', 'G', 'U')] er: char, // right extension
) {
    // 2bp seed, with extension enabled on both sides
    let args = ["-l", "10", "-e", "10000.0", "-s", "2", "-p3"];

    let query = format!("{}{}{}{}", el, s1, s2, er);
    let target = wc_complement(&query);

    SingleSeqRunner::new(&query, &target)
        .assert_pass(&format!("ext_both_{}_{}{}_{}", el, s1, s2, er), &args);
}
