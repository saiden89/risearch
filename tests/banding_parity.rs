//! Experimental banded DP parity checks.
//!
//! These tests compare banded Rust output against the legacy C implementation.
//! They are ignored by default because banding is experimental.

mod common;

use common::{ParityRunner, workspace_root};
use rstest::rstest;

#[rstest]
#[ignore]
fn test_parity_banded_soft(
    #[values(6, 8)] s: usize,
    #[values(10, 20)] l: usize,
    #[values(8, 12)] band: usize,
) {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");

    let l_str = l.to_string();
    let s_str = s.to_string();
    let band_str = band.to_string();

    let args = [
        "-l",
        &l_str,
        "-e",
        "100.0",
        "--seed-length",
        &s_str,
        "-p3",
        "--experimental",
        "--dp-band",
        &band_str,
        "--dp-band-mode",
        "soft",
    ];

    let test_name = format!("band_soft_s{}_l{}_b{}", s, l, band);
    ParityRunner::new(&target).assert_pass(&query, &test_name, &args);
}

#[rstest]
#[ignore]
fn test_parity_banded_hard(
    #[values(6, 8)] s: usize,
    #[values(10, 20)] l: usize,
    #[values(8, 12)] band: usize,
) {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");

    let l_str = l.to_string();
    let s_str = s.to_string();
    let band_str = band.to_string();

    let args = [
        "-l",
        &l_str,
        "-e",
        "100.0",
        "--seed-length",
        &s_str,
        "-p3",
        "--experimental",
        "--dp-band",
        &band_str,
        "--dp-band-mode",
        "hard",
    ];

    let test_name = format!("band_hard_s{}_l{}_b{}", s, l, band);
    ParityRunner::new(&target).assert_pass(&query, &test_name, &args);
}
