use assert_cmd::cargo::cargo_bin_cmd;
use flate2::read::GzDecoder;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Key {
    q_id: String,
    q_start: u32,
    q_end: u32,
    t_id: String,
    t_start: u32,
    t_end: u32,
    strand: String,
}

#[derive(Debug, Clone)]
struct Rec {
    key: Key,
    energy: String,
}

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn parse_first_8_fields(output: &str) -> Vec<Rec> {
    let mut out = Vec::new();

    for (line_no, line) in output.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 8 {
            panic!("Expected >= 8 fields at line {}: {:?}", line_no + 1, line);
        }

        let rec = Rec {
            key: Key {
                q_id: fields[0].to_string(),
                q_start: fields[1].parse().expect("q_start int"),
                q_end: fields[2].parse().expect("q_end int"),
                t_id: fields[3].to_string(),
                t_start: fields[4].parse().expect("t_start int"),
                t_end: fields[5].parse().expect("t_end int"),
                strand: fields[6].to_string(),
            },
            energy: fields[7].to_string(),
        };

        out.push(rec);
    }

    out
}

fn index_and_search_rust(query: &Path, target: &Path, tmp_idx: &Path) -> String {
    // Build index
    let mut cmd = cargo_bin_cmd!("risearch");
    cmd.args(["index", target.to_str().unwrap(), tmp_idx.to_str().unwrap()])
        .assert()
        .success();

    // Create a temp file for output
    let tmp_output = tempfile::NamedTempFile::new().expect("temp output file");
    let output_path = tmp_output.path();

    // Search
    let mut cmd = cargo_bin_cmd!("risearch");
    cmd.args([
        "search",
        "-q",
        query.to_str().unwrap(),
        "-i",
        tmp_idx.to_str().unwrap(),
        "-o",
        output_path.to_str().unwrap(),
        "-l",
        "20",
        "-e",
        "-20",
        "-s",
        "6",
    ])
    .assert()
    .success();

    // Read the output file
    std::fs::read_to_string(output_path).expect("read rust output file")
}

fn search_c(query: &Path, c_index: &Path, c_bin: &Path) -> String {
    // Legacy RIsearch2 writes results to files (e.g. risearch_<query-id>.out.gz)
    // in the current working directory by default, not to stdout.
    let tmpdir = tempfile::tempdir().expect("tempdir");

    let out = std::process::Command::new(c_bin)
        .current_dir(tmpdir.path())
        .args([
            "-q",
            query.to_str().unwrap(),
            "-i",
            c_index.to_str().unwrap(),
            "-l",
            "20",
            "-e",
            "-20",
            "-s",
            "6",
        ])
        .output()
        .expect("run legacy C risearch2");

    if !out.status.success() {
        panic!(
            "Legacy C risearch2 failed: status={:?}\nstdout=\n{}\nstderr=\n{}\n",
            out.status,
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
    }

    let mut out_files: Vec<PathBuf> = fs::read_dir(tmpdir.path())
        .expect("read legacy C output dir")
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with("risearch_") && n.ends_with(".out.gz"))
        })
        .collect();

    out_files.sort();

    if out_files.is_empty() {
        let mut files: Vec<String> = fs::read_dir(tmpdir.path())
            .ok()
            .into_iter()
            .flat_map(|it| it)
            .filter_map(|e| e.ok())
            .filter_map(|e| e.file_name().into_string().ok())
            .collect();
        files.sort();

        panic!(
            "Legacy C risearch2 produced no risearch_*.out.gz files.\n\
status={:?}\n\
stdout=\n{}\n\
stderr=\n{}\n\
cwd={:?}\n\
files={:?}\n",
            out.status,
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr),
            tmpdir.path(),
            files,
        );
    }

    let mut combined = String::new();
    for path in out_files {
        let bytes = fs::read(&path).expect("read legacy output file");
        let mut decoder = GzDecoder::new(&bytes[..]);
        let mut s = String::new();
        decoder
            .read_to_string(&mut s)
            .expect("decode legacy .out.gz as utf8");
        combined.push_str(&s);
        if !combined.ends_with('\n') {
            combined.push('\n');
        }
    }

    combined
}

#[test]
fn parity_against_legacy_c_on_test_suite() {
    let root = workspace_root();

    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");

    let c_index = root.join("legacy_c/RIsearch2/test_suite/RHOC.pksuf");
    let c_bin = root.join("legacy_c/RIsearch2/bin/risearch2.dbg.x");

    assert!(query.exists(), "missing query file: {query:?}");
    assert!(target.exists(), "missing target file: {target:?}");
    assert!(c_index.exists(), "missing C index file: {c_index:?}");
    assert!(c_bin.exists(), "missing C binary: {c_bin:?}");

    let tmpdir = tempfile::tempdir().expect("tempdir");
    let rust_idx = tmpdir.path().join("RHOC.idx");

    let rust_out = index_and_search_rust(&query, &target, &rust_idx);
    let c_out = search_c(&query, &c_index, &c_bin);

    let rust_recs = parse_first_8_fields(&rust_out);
    let c_recs = parse_first_8_fields(&c_out);

    assert!(
        !rust_recs.is_empty(),
        "Rust search returned 0 records; stdout was empty or had only whitespace.\n\
query={query:?}\n\
target={target:?}\n\
rust_idx={rust_idx:?}\n\
stdout=\n{rust_out}\n"
    );
    assert!(
        !c_recs.is_empty(),
        "Legacy C search returned 0 records; stdout was empty or had only whitespace.\n\
query={query:?}\n\
c_index={c_index:?}\n\
c_bin={c_bin:?}\n\
stdout=\n{c_out}\n"
    );

    let rust_map: BTreeMap<Key, String> =
        rust_recs.into_iter().map(|r| (r.key, r.energy)).collect();
    let c_map: BTreeMap<Key, String> = c_recs.into_iter().map(|r| (r.key, r.energy)).collect();

    let rust_keys: BTreeSet<Key> = rust_map.keys().cloned().collect();
    let c_keys: BTreeSet<Key> = c_map.keys().cloned().collect();

    let only_in_rust: Vec<_> = rust_keys.difference(&c_keys).cloned().collect();
    let only_in_c: Vec<_> = c_keys.difference(&rust_keys).cloned().collect();

    let mut energy_mismatches = Vec::new();
    for k in rust_keys.intersection(&c_keys) {
        let r = rust_map.get(k).unwrap();
        let c = c_map.get(k).unwrap();
        if r != c {
            energy_mismatches.push((k.clone(), c.clone(), r.clone()));
        }
    }

    if !only_in_rust.is_empty() || !only_in_c.is_empty() || !energy_mismatches.is_empty() {
        let mut msg = String::new();

        if !only_in_c.is_empty() {
            msg.push_str("Keys only in C (first 20):\n");
            for k in only_in_c.iter().take(20) {
                msg.push_str(&format!("  {k:?}\n"));
            }
        }

        if !only_in_rust.is_empty() {
            msg.push_str("Keys only in Rust (first 20):\n");
            for k in only_in_rust.iter().take(20) {
                msg.push_str(&format!("  {k:?}\n"));
            }
        }

        if !energy_mismatches.is_empty() {
            msg.push_str("Energy mismatches (first 30):\n");
            for (k, c_e, r_e) in energy_mismatches.iter().take(30) {
                msg.push_str(&format!("  {k:?}: C={c_e} Rust={r_e}\n"));
            }
        }

        panic!("Rust/C parity failed:\n{msg}");
    }
}
