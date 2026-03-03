use assert_cmd::cargo::cargo_bin_cmd;
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

fn sorted_non_empty_lines(s: &str) -> Vec<&str> {
    let mut lines: Vec<_> = s.lines().map(str::trim).filter(|l| !l.is_empty()).collect();
    lines.sort_unstable();
    lines
}

fn assert_minimal_hit_line(line: &str) {
    let fields: Vec<&str> = line.split('\t').collect();
    assert_eq!(fields.len(), 8, "expected 8 TSV fields in minimal output");
    assert_eq!(fields[0], "q", "unexpected query name");
    assert_eq!(fields[3], "t", "unexpected target name");
    assert_eq!(fields[6], "+", "unexpected strand");
    assert!(
        fields[1].parse::<u32>().is_ok()
            && fields[2].parse::<u32>().is_ok()
            && fields[4].parse::<u32>().is_ok()
            && fields[5].parse::<u32>().is_ok(),
        "coordinate fields should be integers: {line}"
    );
    assert!(
        fields[7].parse::<f32>().is_ok(),
        "energy field should be numeric: {line}"
    );
}

#[test]
fn test_search_stdout_output() -> Result<(), Box<dyn std::error::Error>> {
    // Deterministic fixture that yields non-empty minimal-format hits.
    let mut query_file = NamedTempFile::new()?;
    writeln!(query_file, ">q\nAAAA")?;

    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, ">t\nUUUUUUUU")?;

    let index_file = NamedTempFile::new()?;
    build_index(target_file.path(), index_file.path());

    // Use a temp working directory so we can assert that no file literally named "-"
    // is created as a side effect.
    let cwd = tempdir()?;
    let dash_path = cwd.path().join("-");

    // Run search with output explicitly directed to stdout.
    let mut cmd_search = cargo_bin_cmd!("risearch");
    let assert = cmd_search
        .current_dir(cwd.path())
        .arg("-vvv")
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
        .arg("--format")
        .arg("minimal")
        .arg("-e")
        .arg("10000.0")
        .assert()
        .success();

    let output = assert.get_output();
    let stdout = String::from_utf8(output.stdout.clone())?;
    let stderr = String::from_utf8(output.stderr.clone())?;

    let stdout_lines = sorted_non_empty_lines(&stdout);
    assert!(!stdout_lines.is_empty(), "expected hits on stdout");
    for line in &stdout_lines {
        assert_minimal_hit_line(line);
    }

    // Ensure stderr contains logs only and no hit-like TSV rows.
    assert!(
        stderr
            .lines()
            .all(|l| !l.starts_with("q\t") && !l.starts_with("t\t")),
        "stderr should not contain hit rows; got:\n{stderr}"
    );

    // Ensure stdout mode does not create a file named "-".
    assert!(
        !dash_path.exists(),
        "search with -o - created an unexpected file named '-'"
    );

    // Cross-check: same args with file output should match stdout content.
    let out_file = NamedTempFile::new()?;
    let mut cmd_file = cargo_bin_cmd!("risearch");
    cmd_file
        .current_dir(cwd.path())
        .arg("-vvv")
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
        .arg("--format")
        .arg("minimal")
        .arg("-e")
        .arg("10000.0")
        .assert()
        .success();

    let file_out = std::fs::read_to_string(out_file.path())?;
    assert_eq!(
        stdout_lines,
        sorted_non_empty_lines(&file_out),
        "stdout output and file output diverged for equivalent search settings"
    );

    Ok(())
}
