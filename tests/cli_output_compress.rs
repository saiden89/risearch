//! Compressed output selected by `--compress` or inferred from the file extension.

use assert_cmd::cargo::cargo_bin_cmd;
use flate2::read::MultiGzDecoder;
use std::io::{Read, Write};
use tempfile::{Builder, NamedTempFile};

fn build_index(target_path: &std::path::Path, index_path: &std::path::Path) {
    let mut cmd_index = cargo_bin_cmd!("risearch");
    cmd_index
        .arg("index")
        .arg(target_path)
        .arg(index_path)
        .assert()
        .success();
}

fn run_search(
    query_path: &std::path::Path,
    target_path: &std::path::Path,
    out_path: &std::path::Path,
    extra_args: &[&str],
) {
    let mut cmd_search = cargo_bin_cmd!("risearch");
    let mut cmd = cmd_search
        .arg("search")
        .arg("-q")
        .arg(query_path)
        .arg("-t")
        .arg(target_path)
        .arg("-o")
        .arg(out_path)
        .arg("-e")
        .arg("100.0")
        .arg("--seed-length")
        .arg("2")
        .arg("-l")
        .arg("0");

    for arg in extra_args {
        cmd = cmd.arg(arg);
    }

    cmd.assert().success();
}

fn assert_gzip_file(path: &std::path::Path) -> Result<(), Box<dyn std::error::Error>> {
    let bytes = std::fs::read(path)?;
    assert!(bytes.len() >= 2, "gzip file too short");
    assert_eq!(bytes[0], 0x1f);
    assert_eq!(bytes[1], 0x8b);

    let mut decoder = MultiGzDecoder::new(&bytes[..]);
    let mut s = String::new();
    decoder.read_to_string(&mut s)?;
    Ok(())
}

#[test]
fn gzip_output_by_extension() -> Result<(), Box<dyn std::error::Error>> {
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q\nAAAAUA")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t\nUUUUAU")?;

    let index_file = NamedTempFile::new()?;
    build_index(target_file.path(), index_file.path());

    let out_file = Builder::new().suffix(".gz").tempfile()?;
    run_search(query_file.path(), index_file.path(), out_file.path(), &[]);

    assert_gzip_file(out_file.path())?;
    Ok(())
}

#[test]
fn gzip_output_by_flag() -> Result<(), Box<dyn std::error::Error>> {
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q\nAAAAUA")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t\nUUUUAU")?;

    let index_file = NamedTempFile::new()?;
    build_index(target_file.path(), index_file.path());

    let out_file = NamedTempFile::new()?;
    run_search(
        query_file.path(),
        index_file.path(),
        out_file.path(),
        &["--compress", "gzip"],
    );

    assert_gzip_file(out_file.path())?;
    Ok(())
}
