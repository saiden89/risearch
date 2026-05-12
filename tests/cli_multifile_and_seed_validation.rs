use assert_cmd::cargo::cargo_bin_cmd;
use predicates::prelude::*;
use std::io::Write;
use std::path::Path;
use tempfile::{tempdir, NamedTempFile};

fn build_index(target_path: &Path, index_path: &Path) {
    let mut cmd_index = cargo_bin_cmd!("risearch");
    cmd_index
        .arg("index")
        .arg(target_path)
        .arg(index_path)
        .assert()
        .success();
}

#[test]
fn seed_length_zero_uses_full_interval() -> Result<(), Box<dyn std::error::Error>> {
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q\nAAAA")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t\nUUUUUUUU")?;

    let index_file = NamedTempFile::new()?;
    build_index(target_file.path(), index_file.path());

    let out_file = NamedTempFile::new()?;
    let mut cmd_search = cargo_bin_cmd!("risearch");
    cmd_search
        .arg("search")
        .arg("-q")
        .arg(query_file.path())
        .arg("-t")
        .arg(index_file.path())
        .arg("-o")
        .arg(out_file.path())
        .arg("--seed-start")
        .arg("1")
        .arg("--seed-end")
        .arg("4")
        .arg("--seed-length")
        .arg("0")
        .arg("-l")
        .arg("0")
        .arg("-e")
        .arg("100.0")
        .arg("--format")
        .arg("minimal")
        .assert()
        .success();

    Ok(())
}

#[test]
fn multifile_output_avoids_sanitized_name_collisions() -> Result<(), Box<dyn std::error::Error>> {
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">a/b\nAAAA\n>a:b\nCCCC")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t1\nUUUUUUUUUU\n>t2\nGGGGGGGGGG")?;

    let index_file = NamedTempFile::new()?;
    build_index(target_file.path(), index_file.path());

    let out_dir = tempdir()?;
    let mut cmd_search = cargo_bin_cmd!("risearch");
    cmd_search
        .arg("search")
        .arg("-q")
        .arg(query_file.path())
        .arg("-t")
        .arg(index_file.path())
        .arg("-o")
        .arg(out_dir.path())
        .arg("--multifile")
        .arg("--seed-length")
        .arg("1")
        .arg("-l")
        .arg("0")
        .arg("-e")
        .arg("10000.0")
        .arg("--format")
        .arg("minimal")
        .assert()
        .success();

    let mut files = std::fs::read_dir(out_dir.path())?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .collect::<Vec<_>>();
    files.sort();

    assert_eq!(files.len(), 2, "expected one output file per query");
    assert_ne!(files[0].file_name(), files[1].file_name());
    assert!(files.iter().any(|p| {
        p.file_name()
            .is_some_and(|n| n.to_string_lossy().starts_with("a_b.tsv"))
    }));
    assert!(files.iter().all(|p| {
        std::fs::read_to_string(p)
            .map(|s| !s.trim().is_empty())
            .unwrap_or(false)
    }));

    Ok(())
}

#[test]
fn multifile_rejects_stdout_output() -> Result<(), Box<dyn std::error::Error>> {
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q\nAAAA")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t\nUUUUUUUU")?;

    let index_file = NamedTempFile::new()?;
    build_index(target_file.path(), index_file.path());

    let mut cmd_search = cargo_bin_cmd!("risearch");
    cmd_search
        .arg("search")
        .arg("-q")
        .arg(query_file.path())
        .arg("-t")
        .arg(index_file.path())
        .arg("-o")
        .arg("-")
        .arg("--multifile")
        .arg("--seed-length")
        .arg("1")
        .arg("-l")
        .arg("0")
        .arg("-e")
        .arg("10000.0")
        .arg("--format")
        .arg("minimal")
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "--multifile requires -o/--output to be a directory path",
        ));

    Ok(())
}

#[test]
fn bindingsite_output_includes_non_empty_flanks() -> Result<(), Box<dyn std::error::Error>> {
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q\nAAAA")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t\nUUUUUUUUUUUUUUUUUUUUUUUUUUUU")?;

    let index_file = NamedTempFile::new()?;
    build_index(target_file.path(), index_file.path());

    let out_file = NamedTempFile::new()?;
    let mut cmd_search = cargo_bin_cmd!("risearch");
    cmd_search
        .arg("search")
        .arg("-q")
        .arg(query_file.path())
        .arg("-t")
        .arg(index_file.path())
        .arg("-o")
        .arg(out_file.path())
        .arg("--seed-length")
        .arg("2")
        .arg("-l")
        .arg("0")
        .arg("-e")
        .arg("10000.0")
        .arg("--format")
        .arg("bindingsite")
        .assert()
        .success();

    let content = std::fs::read_to_string(out_file.path())?;
    assert!(!content.trim().is_empty(), "expected at least one hit");
    let has_non_empty_flank = content.lines().any(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        fields.len() >= 12 && (!fields[10].is_empty() || !fields[11].is_empty())
    });
    assert!(
        has_non_empty_flank,
        "expected bindingsite output to include non-empty target flanks"
    );

    Ok(())
}
