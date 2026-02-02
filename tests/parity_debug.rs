//! Debug/regression parity tests: isolated reproductions for specific issues.

mod common;

use common::SingleSeqRunner;

#[test]
fn test_alignment_repro() {
    let query = "ucaguucagcaggaacag";
    let target = "TGGCTCTGTGGGACACAGCAGG";
    let args = ["-l", "20", "-e", "-20", "--seed-length", "6", "-p3"];
    SingleSeqRunner::new(query, target).assert_pass("alignment_repro", &args);
}

#[test]
fn test_custom_seq() {
    let args = ["-l", "20", "-e", "10000", "--seed-length", "7", "-p3"];
    let query = "uggcucaguucagcaggaacag";
    let target = "ggaagaccgacuaggagacgacuugcugcuacuccuccgucccugcaucuggaggccuuu";
    SingleSeqRunner::new(query, target).assert_pass("custom_seq", &args);
}

#[test]
fn test_mir24_isolated() {
    let query = "uggcucaguucagcaggaacag";
    let target = "GGAAGACCTGCCTCCTCATCGTCTTCAGCAAGGATCAGTTTCCGGAGGTCTACGTCCCTACTGTCTTTGAGAACTATATTG";
    let args = ["-l", "20", "-e", "100.0", "--seed-length", "6", "-p3", "--no-max-prune"];
    SingleSeqRunner::new(query, target).assert_pass("mir24_isolated", &args);
}

#[test]
fn test_segments_debug() {
    let query = "AAAAAUGCUGUAAAAA";
    let target = "UUUUUACAGCAUUUUU";
    let args = ["-l", "20", "-e", "100.0", "--seed-length", "6", "-p3"];
    SingleSeqRunner::new(query, target).assert_pass("segments", &args);
}

#[test]
fn test_rust_worse_debug() {
    let query = "uggcucaguucagcaggaacag";
    let target = "ATATTGCGGACATTGAGGTGGACGGCAAGCAGGTGGAGCTGGCTCTGTGGGACACAGCAGGGCAGGAAGACTATGATCGACTGCGGCCTCTCTCCTACCCG";
    let args = ["-l", "20", "-e", "100.0", "--seed-length", "6", "-p3"];
    SingleSeqRunner::new(query, target).assert_pass("rust_worse", &args);
}

#[test]
fn test_missing_544() {
    let query = "uggcucaguucagcaggaacag";
    let target = "TGATTGCCCTCCATCAACACTGCCCACCCCAGGTTGGGGCTACCCCAGCCCATCTTTACAAAACAGGGCAAGGTGAACTAATGGAGTGGGTGGAGGAGTTGGAAGA";
    let args = ["-l", "20", "-e", "100.0", "--seed-length", "6", "-p3"];
    SingleSeqRunner::new(query, target).assert_pass("missing_544", &args);
}

#[test]
fn test_energy_discrepancy() {
    let query = "cucaguucagcaggaacag";
    let target = "CTGACTCCTTGCCCTGAGT";
    let args = ["-l", "10", "-e", "-10.0", "--seed-length", "5", "-p3"];
    SingleSeqRunner::new(query, target).assert_pass("energy_discrepancy", &args);
}
