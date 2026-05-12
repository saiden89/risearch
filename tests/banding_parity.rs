//! Experimental banded DP parity checks.
//!
//! These tests compare banded Rust output against the legacy C implementation.
//! They are ignored by default because banding is experimental.

mod support;

use rstest::rstest;
use std::path::PathBuf;
use support::{ParityRunner, query_fa, target_fa};

#[rstest]
#[ignore]
fn banded_soft(
    query_fa: PathBuf,
    target_fa: PathBuf,
    #[values(6, 8)] s: usize,
    #[values(10, 20)] l: usize,
    #[values(8, 12)] band: usize,
) {

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
    ParityRunner::new(&target_fa).assert_pass(&query_fa, &test_name, &args);
}

#[rstest]
#[ignore]
fn banded_hard(
    query_fa: PathBuf,
    target_fa: PathBuf,
    #[values(6, 8)] s: usize,
    #[values(10, 20)] l: usize,
    #[values(8, 12)] band: usize,
) {

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
    ParityRunner::new(&target_fa).assert_pass(&query_fa, &test_name, &args);
}
