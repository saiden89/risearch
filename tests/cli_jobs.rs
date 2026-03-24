use assert_cmd::cargo::cargo_bin_cmd;
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

fn run_search(
    query_path: &Path,
    target_index_path: &Path,
    jobs: usize,
) -> Result<String, Box<dyn std::error::Error>> {
    let out_file = NamedTempFile::new()?;

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
        .arg(out_file.path())
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

    Ok(std::fs::read_to_string(out_file.path())?)
}

#[test]
fn jobs_flag_produces_consistent_streamed_output() -> Result<(), Box<dyn std::error::Error>> {
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q1\nAAAA\n>q2\nCCCC")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t1\nTTTT\n>t2\nGGGG")?;

    let index_file = NamedTempFile::new()?;
    build_index(target_file.path(), index_file.path());

    let out_j1 = run_search(query_file.path(), index_file.path(), 1)?;
    let out_j4 = run_search(query_file.path(), index_file.path(), 4)?;

    fn sort_lines(s: String) -> String {
        let mut lines: Vec<_> = s.lines().collect();
        lines.sort_unstable();
        lines.join("\n")
    }

    assert_eq!(
        sort_lines(out_j1),
        sort_lines(out_j4),
        "search output changed across --jobs values"
    );
    Ok(())
}
