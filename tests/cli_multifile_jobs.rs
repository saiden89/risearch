//! Multi-file output stays consistent across `-j/--jobs` values.

use assert_cmd::cargo::cargo_bin_cmd;
use std::collections::BTreeMap;
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

fn run_search_multifile(query_path: &Path, target_index_path: &Path, out_dir: &Path, jobs: usize) {
    let mut cmd_search = cargo_bin_cmd!("risearch");
    cmd_search
        .arg("-j")
        .arg(jobs.to_string())
        .arg("search")
        .arg("-q")
        .arg(query_path)
        .arg("-t")
        .arg(target_index_path)
        .arg("-o")
        .arg(out_dir)
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
}

fn sorted_lines(text: &str) -> String {
    let mut lines: Vec<&str> = text.lines().collect();
    lines.sort_unstable();
    lines.join("\n")
}

fn read_multifile_dir(dir: &Path) -> Result<BTreeMap<String, String>, Box<dyn std::error::Error>> {
    let mut out = BTreeMap::new();
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        if !entry.file_type()?.is_file() {
            continue;
        }
        let name = entry.file_name().to_string_lossy().to_string();
        let data = std::fs::read_to_string(entry.path())?;
        out.insert(name, sorted_lines(&data));
    }
    Ok(out)
}

#[test]
fn multifile_output_consistent_across_jobs() -> Result<(), Box<dyn std::error::Error>> {
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q1\nAAAA\n>q2\nCCCC")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t1\nTTTT\n>t2\nGGGG")?;

    let index_file = NamedTempFile::new()?;
    build_index(target_file.path(), index_file.path());

    let out_j1 = tempdir()?;
    let out_j4 = tempdir()?;
    run_search_multifile(query_file.path(), index_file.path(), out_j1.path(), 1);
    run_search_multifile(query_file.path(), index_file.path(), out_j4.path(), 4);

    let files_j1 = read_multifile_dir(out_j1.path())?;
    let files_j4 = read_multifile_dir(out_j4.path())?;

    assert_eq!(
        files_j1, files_j4,
        "multifile output changed across --jobs values"
    );
    Ok(())
}
