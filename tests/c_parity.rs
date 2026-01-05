//! C Parity Tests
//!
//! These tests compare Rust risearch output against the legacy C risearch2 implementation.
//! Run with RUST_LOG=debug for detailed diff analysis, or RUST_LOG=trace for raw data.

mod common;

use common::{
    LogTag, c_binary_path, compare_results, create_c_index, index_and_search_rust,
    run_single_seq_parity, search_c, setup_common_test_files, workspace_root,
};
use log::{debug, info, trace};
use std::fs;

// =============================================================================
// TESTS
// =============================================================================

#[test]
fn test_parity_full_pipeline() {
    let (tmpdir, query_path, target_path, c_bin) = setup_common_test_files();
    let c_index = tmpdir.path().join("c_target.pksuf");
    let rust_idx = tmpdir.path().join("rust_target.idx");

    create_c_index(&target_path, &c_index, &c_bin);

    let args = ["-l", "20", "-e", "-20", "-s", "6", "-p3"];

    let rust_output = index_and_search_rust(&query_path, &target_path, &rust_idx, &args);
    let rust_out = String::from_utf8_lossy(&rust_output.stdout).to_string();

    // Log C_DEBUG traces at trace level
    for line in rust_out.lines() {
        if line.contains("C_DEBUG:") {
            if line.contains("final=-20.49") || line.contains("final=-21.50") {
                trace!("{} {}", LogTag::Rust, line);
            }
        }
    }

    let c_out = search_c(&query_path, &c_index, &c_bin, &args);
    for line in c_out.lines() {
        if line.contains("C_DEBUG:") {
            if line.contains("final=-20.49") || line.contains("final=-21.50") {
                trace!("{} {}", LogTag::C, line);
            }
        }
    }

    compare_results(&rust_out, &c_out, "default_config");
}

#[test]
fn test_parity_long_seed() {
    // User requested: -l 0 -e 10000 -s 12 -p 3
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");
    let c_index = root.join("legacy_c/RIsearch2/test_suite/RHOC.pksuf");
    let c_bin = c_binary_path(&root);
    let tmpdir = tempfile::tempdir().expect("tempdir");
    let rust_idx = tmpdir.path().join("RHOC.idx");

    let args = ["-l", "0", "-e", "10000", "-s", "12", "-p3"];

    let rust_output = index_and_search_rust(&query, &target, &rust_idx, &args);
    let rust_out = String::from_utf8_lossy(&rust_output.stdout).to_string();
    let c_out = search_c(&query, &c_index, &c_bin, &args);

    compare_results(&rust_out, &c_out, "long_seed_no_ext");
}

#[test]
fn test_parity_energy_threshold() {
    let (tmpdir, query_path, target_path, c_bin) = setup_common_test_files();
    let c_index = tmpdir.path().join("c_target.pksuf");
    let rust_idx = tmpdir.path().join("rust_target.idx");

    create_c_index(&target_path, &c_index, &c_bin);

    let args = ["-l", "10", "-e", "-10.0", "-s", "5", "-p3"];

    let rust_output = index_and_search_rust(&query_path, &target_path, &rust_idx, &args);
    let rust_out = String::from_utf8_lossy(&rust_output.stdout).to_string();
    let c_out = search_c(&query_path, &c_index, &c_bin, &args);

    compare_results(&rust_out, &c_out, "energy_only");
}

#[test]
fn test_parity_alignment_repro() {
    let root = workspace_root();
    let tmpdir = tempfile::tempdir().expect("tempdir");

    // Create inputs
    let query_path = tmpdir.path().join("query.fa");
    let target_path = tmpdir.path().join("target.fa");

    fs::write(&query_path, ">query\nuggcucaguucagcaggaacag\n").expect("write query");
    fs::write(&target_path, ">target\nTGGCTCTGTGGGACACAGCAGG\n").expect("write target");

    let c_bin = c_binary_path(&root);
    let c_index = tmpdir.path().join("target.pksuf");

    let c_index_cmd = std::process::Command::new(&c_bin)
        .arg("-c")
        .arg(target_path.to_str().unwrap())
        .arg("-o")
        .arg(c_index.to_str().unwrap())
        .output()
        .expect("c index creation");

    if !c_index_cmd.status.success() {
        panic!("C indexing failed");
    }

    let rust_idx = tmpdir.path().join("target.idx");

    let args = ["-l", "20", "-e", "-20", "-s", "6", "-p3"];

    let rust_output = index_and_search_rust(&query_path, &target_path, &rust_idx, &args);
    let rust_out = String::from_utf8_lossy(&rust_output.stdout).to_string();
    let c_out = search_c(&query_path, &c_index, &c_bin, &args);

    debug!("{} Rust Output:\n{}", LogTag::Parity, rust_out);
    debug!("{} C Output:\n{}", LogTag::Parity, c_out);

    if rust_out.trim() == c_out.trim() {
        info!(
            "{} alignment_mismatch_repro - Outputs identical",
            LogTag::Parity
        );
    } else {
        debug!("{} Outputs differ", LogTag::Parity);
    }

    compare_results(&rust_out, &c_out, "alignment_mismatch_repro");
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

    info!("{} mir24_isolated_single test", LogTag::Parity);
    debug!("{} Query:  {}", LogTag::Parity, query);
    debug!("{} Target: {}", LogTag::Parity, target);
    debug!("{} Args:   {:?}", LogTag::Parity, args);

    run_single_seq_parity(query, target, "miR24_isolated_single", &args, true);
}

/// More targeted single-hit test with explicit segment expectations.
#[test]
fn test_parity_segments_debug() {
    let query = "AAAAAUGCUGUAAAAA"; // 5 + 6 + 5 = 16nt
    let target = "UUUUUACAGCAUUUUU";

    let args = ["-l", "20", "-e", "100.0", "-s", "6", "-p3"];

    info!("{} segment_debug test", LogTag::Parity);
    debug!("{} Query:  {} (len={})", LogTag::Parity, query, query.len());
    debug!(
        "{} Target: {} (len={})",
        LogTag::Parity,
        target,
        target.len()
    );
    debug!("{} Expected structure:", LogTag::Parity);
    debug!(
        "{}   LEFT_EXT  (dp_left):  AAAAA -> should extend with UUUUU",
        LogTag::Parity
    );
    debug!(
        "{}   SEED      (6nt):      UGCUGU -> matches ACAGCA",
        LogTag::Parity
    );
    debug!(
        "{}   RIGHT_EXT (dp_right): AAAAA -> should extend with UUUUU",
        LogTag::Parity
    );

    run_single_seq_parity(query, target, "segment_debug", &args, false);
}
