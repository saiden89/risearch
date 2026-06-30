//! Minimal reproducer for strict mismatch parity divergence seen in benchmark.

mod support;

use flate2::read::GzDecoder;
use std::io::Read;
use std::io::Write;
use std::process::Command;
use support::{workspace_root, SingleSeqRunner};
use tempfile::tempdir;

#[test]
fn test_release_c_strict_m1_seed17_first_pos_mismatch() -> Result<(), Box<dyn std::error::Error>> {
    let root = workspace_root();
    let c_bin = root.join("legacy_c/RIsearch2/bin/risearch2.x");
    assert!(
        c_bin.exists(),
        "Missing C release binary at {}",
        c_bin.display()
    );

    let tmp = tempdir()?;
    let query_path = tmp.path().join("q.fa");
    let target_path = tmp.path().join("t.fa");
    let rust_index = tmp.path().join("t.idx");
    let c_index = tmp.path().join("t.cidx");

    // Exact minimal pair extracted from the benchmark delta.
    let mut qf = std::fs::File::create(&query_path)?;
    writeln!(qf, ">q\nCTTCCAGCAACCCCAGG")?;
    let mut tf = std::fs::File::create(&target_path)?;
    writeln!(tf, ">t\nCCTGGGGTTGCTGGAAA")?;

    let rust_bin = assert_cmd::cargo::cargo_bin!("risearch");
    let index_status = Command::new(rust_bin)
        .arg("index")
        .arg(&target_path)
        .arg(&rust_index)
        .status()?;
    assert!(index_status.success(), "Rust index build failed");

    let c_index_status = Command::new(&c_bin)
        .arg("-c")
        .arg(&target_path)
        .arg("-o")
        .arg(&c_index)
        .status()?;
    assert!(c_index_status.success(), "C index build failed");

    let rust_out = Command::new(rust_bin)
        .arg("search")
        .arg("-q")
        .arg(&query_path)
        .arg("-t")
        .arg(&rust_index)
        .arg("--seed-length")
        .arg("17")
        .arg("-l")
        .arg("0")
        .arg("--mismatch-max")
        .arg("1")
        .arg("--mismatch-prefix")
        .arg("0")
        .arg("--mismatch-suffix")
        .arg("0")
        .arg("--no-seed-wobble")
        .arg("--no-dedup")
        .arg("-j")
        .arg("1")
        .arg("-e")
        .arg("10000")
        .arg("--format")
        .arg("minimal")
        .arg("-o")
        .arg("-")
        .output()?;
    assert!(rust_out.status.success(), "Rust search failed");

    let c_out = Command::new(&c_bin)
        .current_dir(tmp.path())
        .arg("-q")
        .arg(&query_path)
        .arg("-i")
        .arg(&c_index)
        .arg("-s")
        .arg("17")
        .arg("-l")
        .arg("0")
        .arg("-m")
        .arg("1:0")
        .arg("-e")
        .arg("10000")
        .arg("--noGUseed")
        .arg("-t")
        .arg("1")
        .arg("-p4")
        .output()?;
    assert!(
        c_out.status.success(),
        "C search failed: status={:?}\nstderr={}\nstdout={}",
        c_out.status,
        String::from_utf8_lossy(&c_out.stderr),
        String::from_utf8_lossy(&c_out.stdout)
    );

    let rust_count = String::from_utf8_lossy(&rust_out.stdout)
        .lines()
        .filter(|l| !l.trim().is_empty())
        .count();
    let mut c_files: Vec<_> = std::fs::read_dir(tmp.path())?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with("risearch_") && n.ends_with(".out.gz"))
        })
        .collect();
    c_files.sort();

    let mut c_count = 0usize;
    for path in c_files {
        let bytes = std::fs::read(path)?;
        let mut decoder = GzDecoder::new(&bytes[..]);
        let mut s = String::new();
        decoder.read_to_string(&mut s)?;
        c_count += s.lines().filter(|l| !l.trim().is_empty()).count();
    }

    assert_eq!(
        rust_count, c_count,
        "Mismatch on minimal repro (query=CTTCCAGCAACCCCAGG target=CCTGGGGTTGCTGGAAA): rust={} c={}",
        rust_count, c_count
    );

    Ok(())
}

#[test]
#[ignore = "Run with PARITY_C_BIN=release RUST_LOG=debug to get parity diff tables on release-C mismatch"]
fn test_release_c_strict_m1_seed17_table_repro() {
    let query = "CTTCCAGCAACCCCAGG";
    let target = "CCTGGGGTTGCTGGAAA";
    let args = [
        "-l",
        "0",
        "-e",
        "10000.0",
        "--seed-length",
        "17",
        "-m",
        "1:0",
        "--no-seed-wobble",
        "-p3",
    ];
    SingleSeqRunner::new(query, target).assert_pass("release_c_strict_m1_seed17_table", &args);
}
