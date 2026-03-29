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
fn default_seed_wobble_matches_enabled_behavior() -> Result<(), Box<dyn std::error::Error>> {
    // This fixture has only G-U seed matches, so disabling seed wobble should produce no hits.
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q\nGGGG")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t\nUUUUUUUU")?;

    let index_file = NamedTempFile::new()?;
    build_index(target_file.path(), index_file.path());

    let default_out = NamedTempFile::new()?;
    let strict_out = NamedTempFile::new()?;

    let default_content = run_search_and_read(
        query_file.path(),
        index_file.path(),
        default_out.path(),
        &[],
    )?;
    let strict_content = run_search_and_read(
        query_file.path(),
        index_file.path(),
        strict_out.path(),
        &["--no-seed-wobble"],
    )?;

    let default_count = default_content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .count();
    let strict_count = strict_content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .count();

    assert!(
        default_count > 0,
        "default mode should permit wobble seed hits"
    );
    assert!(
        strict_count < default_count,
        "strict mode should reduce hits for wobble-only seeds"
    );

    Ok(())
}

#[test]
fn legacy_no_guseed_warning_points_to_no_seed_wobble() -> Result<(), Box<dyn std::error::Error>> {
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
            "'--noGUseed' is deprecated; use --no-seed-wobble instead.",
        ));

    Ok(())
}
