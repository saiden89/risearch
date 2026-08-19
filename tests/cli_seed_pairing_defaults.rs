use assert_cmd::cargo::cargo_bin_cmd;
use predicates::prelude::*;
use std::io::Write;
use std::path::Path;
use tempfile::NamedTempFile;

fn build_index(target_path: &Path, index_path: &Path) {
    let mut cmd_index = cargo_bin_cmd!("risearch");
    cmd_index
        .arg("index")
        .arg(target_path)
        .arg(index_path)
        .assert()
        .success();
}

fn run_search_and_read(
    query: &Path,
    index: &Path,
    output: &Path,
    extra_args: &[&str],
) -> Result<String, Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("risearch");
    cmd.arg("search")
        .arg("-q")
        .arg(query)
        .arg("-t")
        .arg(index)
        .arg("-o")
        .arg(output)
        .arg("--seed-length")
        .arg("2")
        .arg("-l")
        .arg("0")
        .arg("-e")
        .arg("10000.0")
        .arg("--format")
        .arg("minimal")
        .args(extra_args)
        .assert()
        .success();
    Ok(std::fs::read_to_string(output)?)
}

#[test]
fn seed_wobble_opt_in_enables_wobble_only_hits() -> Result<(), Box<dyn std::error::Error>> {
    // This fixture has only G-U seed matches, so hits appear only when wobble is opted in.
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q\nGGGG")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t\nUUUUUUUU")?;

    let index_file = NamedTempFile::new()?;
    build_index(target_file.path(), index_file.path());

    let default_out = NamedTempFile::new()?;
    let wobble_out = NamedTempFile::new()?;

    let default_content = run_search_and_read(
        query_file.path(),
        index_file.path(),
        default_out.path(),
        &[],
    )?;
    let wobble_content = run_search_and_read(
        query_file.path(),
        index_file.path(),
        wobble_out.path(),
        &["--seed-wobble"],
    )?;

    let default_count = default_content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .count();
    let wobble_count = wobble_content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .count();

    assert!(
        wobble_count > 0,
        "--seed-wobble should permit wobble seed hits"
    );
    assert!(
        default_count < wobble_count,
        "default mode should reduce hits for wobble-only seeds"
    );

    Ok(())
}

#[test]
fn legacy_no_guseed_warns_it_is_now_a_noop() -> Result<(), Box<dyn std::error::Error>> {
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q\nGGGG")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t\nUUUUUUUU")?;

    let index_file = NamedTempFile::new()?;
    build_index(target_file.path(), index_file.path());

    cargo_bin_cmd!("risearch")
        .arg("search")
        .arg("-q")
        .arg(query_file.path())
        .arg("-t")
        .arg(index_file.path())
        .arg("-o")
        .arg("-")
        .arg("--seed-length")
        .arg("2")
        .arg("-l")
        .arg("0")
        .arg("-e")
        .arg("10000.0")
        .arg("--format")
        .arg("minimal")
        .arg("--noGUseed")
        .assert()
        .success()
        .stderr(predicate::str::contains(
            "'--noGUseed' is deprecated and now a no-op; G-U wobble is off by default, enable it with --seed-wobble.",
        ));

    Ok(())
}
