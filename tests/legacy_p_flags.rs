//! Deprecation warnings for the legacy `-p` and `-i` flags.

use assert_cmd::cargo::cargo_bin_cmd;
use predicates::prelude::*;
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn legacy_p_flags_emit_deprecation_warning() -> Result<(), Box<dyn std::error::Error>> {
    // 1. Setup
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q\nAAAAUA")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t\nUUUUAU")?;

    let index_file = NamedTempFile::new()?;

    // 2. Build Index
    cargo_bin_cmd!("risearch")
        .arg("index")
        .arg(target_file.path())
        .arg(index_file.path())
        .assert()
        .success();

    // 3. Test -p2 (Cigar) warning
    cargo_bin_cmd!("risearch")
        .arg("search")
        .arg("-q")
        .arg(query_file.path())
        .arg("-t")
        .arg(index_file.path())
        .arg("-o")
        .arg("-")
        .arg("-p2") // Legacy flag
        .assert()
        .success()
        .stderr(predicate::str::contains(
            "Legacy -p/--report-alignment is deprecated; use --format cigar.",
        ));

    // 4. Test -p (Detailed) warning
    cargo_bin_cmd!("risearch")
        .arg("search")
        .arg("-q")
        .arg(query_file.path())
        .arg("-t")
        .arg(index_file.path())
        .arg("-o")
        .arg("-")
        .arg("-p") // Legacy flag default
        .assert()
        .success()
        .stderr(predicate::str::contains(
            "Legacy -p/--report-alignment is deprecated; use --format detailed.",
        ));

    Ok(())
}

#[test]
fn legacy_i_flag_emits_deprecation_warning() -> Result<(), Box<dyn std::error::Error>> {
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q\nAAAAUA")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t\nUUUUAU")?;

    let index_file = NamedTempFile::new()?;

    cargo_bin_cmd!("risearch")
        .arg("index")
        .arg(target_file.path())
        .arg(index_file.path())
        .assert()
        .success();

    cargo_bin_cmd!("risearch")
        .arg("search")
        .arg("-q")
        .arg(query_file.path())
        .arg("-i")
        .arg(index_file.path())
        .arg("-o")
        .arg("-")
        .arg("--seed-length")
        .arg("2")
        .arg("-l")
        .arg("0")
        .arg("-e")
        .arg("100.0")
        .assert()
        .success()
        .stderr(predicate::str::contains(
            "'-i' is deprecated; use -t/--target instead.",
        ));

    Ok(())
}
