//! DSM parity tests: dinucleotide/trinucleotide stacking energy calculations.

mod common;

use common::SingleSeqRunner;
use rstest::rstest;

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

/// All 16 dinucleotide stacks (one stacking interaction).
#[rstest]
fn test_dinucleotide(
    #[values('A', 'C', 'G', 'U')] b1: char,
    #[values('A', 'C', 'G', 'U')] b2: char,
) {
    let args = ["-l", "0", "-e", "10000.0", "--seed-length", "2", "-p3"];
    let query = format!("{}{}", b1, b2);
    let target = wc_complement(&query);
    SingleSeqRunner::new(&query, &target).assert_pass(&format!("{}_{}", b1, b2), &args);
}

/// All 64 trinucleotide combinations (two stacking interactions).
#[rstest]
fn test_trinucleotide(
    #[values('A', 'C', 'G', 'U')] b1: char,
    #[values('A', 'C', 'G', 'U')] b2: char,
    #[values('A', 'C', 'G', 'U')] b3: char,
) {
    let args = ["-l", "0", "-e", "10000.0", "--seed-length", "3", "-p3"];
    let query = format!("{}{}{}", b1, b2, b3);
    let target = wc_complement(&query);
    SingleSeqRunner::new(&query, &target).assert_pass(&format!("{}_{}_{}", b1, b2, b3), &args);
}

/// Minimal right extension: 2bp seed + 1bp extension.
#[rstest]
fn test_min_ext_right(
    #[values('A', 'C', 'G', 'U')] s1: char,
    #[values('A', 'C', 'G', 'U')] s2: char,
    #[values('A', 'C', 'G', 'U')] e1: char,
) {
    let args = ["-l", "10", "-e", "10000.0", "--seed-length", "2", "-p3"];
    let query = format!("{}{}{}", s1, s2, e1);
    let target = wc_complement(&query);
    SingleSeqRunner::new(&query, &target).assert_pass(&format!("r_{}_{}_{}", s1, s2, e1), &args);
}

/// Minimal left extension: 1bp extension + 2bp seed.
#[rstest]
fn test_min_ext_left(
    #[values('A', 'C', 'G', 'U')] e1: char,
    #[values('A', 'C', 'G', 'U')] s1: char,
    #[values('A', 'C', 'G', 'U')] s2: char,
) {
    let args = ["-l", "10", "-e", "10000.0", "--seed-length", "2", "-p3"];
    let query = format!("{}{}{}", e1, s1, s2);
    let target = wc_complement(&query);
    SingleSeqRunner::new(&query, &target).assert_pass(&format!("l_{}_{}_{}", e1, s1, s2), &args);
}

/// Both sides extension: 1bp left + 2bp seed + 1bp right.
#[rstest]
fn test_ext_both(
    #[values('A', 'C', 'G', 'U')] el: char,
    #[values('A', 'C', 'G', 'U')] s1: char,
    #[values('A', 'C', 'G', 'U')] s2: char,
    #[values('A', 'C', 'G', 'U')] er: char,
) {
    let args = ["-l", "10", "-e", "10000.0", "--seed-length", "2", "-p3"];
    let query = format!("{}{}{}{}", el, s1, s2, er);
    let target = wc_complement(&query);
    SingleSeqRunner::new(&query, &target)
        .assert_pass(&format!("b_{}_{}{}_{}", el, s1, s2, er), &args);
}
