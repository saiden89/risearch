use assert_cmd::cargo::cargo_bin_cmd;
use flate2::read::GzDecoder;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Rec {
    // Identity fields for grouping
    q_id: String,
    t_id: String,

    // Coordinate fields
    q_start: usize,
    q_end: usize,
    t_start: usize,
    t_end: usize,
    strand: String,

    // Data fields
    energy: String, // Keep as string for exact comparison, parse for fuzzy
    interaction: String,
    target_seq: String,
    flank_5: String,
    flank_3: String,
}

impl Rec {
    fn from_line(line: &str) -> Option<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return None;
        }

        Some(Rec {
            q_id: fields[0].to_string(),
            q_start: fields[1].parse().unwrap_or(0),
            q_end: fields[2].parse().unwrap_or(0),
            t_id: fields[3].to_string(),
            t_start: fields[4].parse().unwrap_or(0),
            t_end: fields[5].parse().unwrap_or(0),
            strand: fields[6].to_string(),
            energy: fields[7].to_string(),
            interaction: fields[8].to_string(),
            target_seq: fields[9].to_string(),
            flank_5: fields[10].to_string(),
            flank_3: fields[11].to_string(),
        })
    }

    /// Returns true if coordinates and strand match
    fn coords_match(&self, other: &Self) -> bool {
        self.q_start == other.q_start
            && self.q_end == other.q_end
            && self.t_start == other.t_start
            && self.t_end == other.t_end
            && self.strand == other.strand
    }
}

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn parse_output(output: &str) -> Vec<Rec> {
    output
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty())
        .filter_map(Rec::from_line)
        .collect()
}

fn index_and_search_rust(query: &Path, target: &Path, tmp_idx: &Path, args: &[&str]) -> String {
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
    let mut final_args = vec![
        "search",
        "-q",
        query.to_str().unwrap(),
        "-i",
        tmp_idx.to_str().unwrap(),
        "-o",
        output_path.to_str().unwrap(),
    ];
    final_args.extend_from_slice(args);

    cmd.args(&final_args).assert().success();

    // Read the output file
    std::fs::read_to_string(output_path).expect("read rust output file")
}

fn search_c(query: &Path, c_index: &Path, c_bin: &Path, args: &[&str]) -> String {
    // Legacy RIsearch2 writes results to files (e.g. risearch_<query-id>.out.gz)
    // in the current working directory by default.
    let tmpdir = tempfile::tempdir().expect("tempdir");

    let mut final_args = vec![
        "-q",
        query.to_str().unwrap(),
        "-i",
        c_index.to_str().unwrap(),
    ];
    final_args.extend_from_slice(args);

    let out = std::process::Command::new(c_bin)
        .current_dir(tmpdir.path())
        .args(&final_args)
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

    if out_files.is_empty() {
        panic!(
            "Legacy C risearch2 produced no output files!\nArgs: {:?}\nStderr:\n{}\nStdout:\n{}\n",
            final_args,
            String::from_utf8_lossy(&out.stderr),
            String::from_utf8_lossy(&out.stdout)
        );
    }

    out_files.sort();

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

    // If no output files, check if stdout was used (legacy behavior varies, but typically files).
    // In our tests, we use -q and -i, which typically generate files.
    if combined.trim().is_empty() {
        panic!(
            "Legacy C risearch2 produced empty output!\nArgs: {:?}\nStderr:\n{}\nStdout:\n{}\n",
            final_args,
            String::from_utf8_lossy(&out.stderr),
            String::from_utf8_lossy(&out.stdout)
        );
    }

    if combined.is_empty() {
        String::from_utf8_lossy(&out.stdout).to_string()
    } else {
        combined
    }
}

fn compare_results(rust_out: &str, c_out: &str, test_name: &str) {
    // 1. Strict Byte-Level Parity Check
    if rust_out.trim() == c_out.trim() {
        return;
    }

    // 2. Granular Reporting
    let rust_recs = parse_output(rust_out);
    let c_recs = parse_output(c_out);

    let mut report = String::new();
    report.push_str(&format!("\n=== PARITY FAILURE REPORT: {} ===\n", test_name));
    report.push_str(&format!("Rust Count: {}\n", rust_recs.len()));
    report.push_str(&format!("C Count:    {}\n\n", c_recs.len()));

    // Group by (q_id, t_id)
    let mut keys = HashSet::new();
    for r in &rust_recs {
        keys.insert((r.q_id.clone(), r.t_id.clone()));
    }
    for r in &c_recs {
        keys.insert((r.q_id.clone(), r.t_id.clone()));
    }
    let mut sorted_keys: Vec<_> = keys.into_iter().collect();
    sorted_keys.sort();

    let mut missing_count = 0;
    let mut extra_count = 0;
    let mut mismatch_count = 0;

    for (q, t) in sorted_keys {
        let mut r_group: Vec<&Rec> = rust_recs
            .iter()
            .filter(|r| r.q_id == q && r.t_id == t)
            .collect();
        r_group.sort_by(|a, b| {
            a.q_start
                .cmp(&b.q_start)
                .then(a.q_end.cmp(&b.q_end))
                .then(a.t_start.cmp(&b.t_start))
                .then(a.t_end.cmp(&b.t_end))
        });
        let mut c_group: Vec<&Rec> = c_recs
            .iter()
            .filter(|r| r.q_id == q && r.t_id == t)
            .collect();
        c_group.sort_by(|a, b| {
            a.q_start
                .cmp(&b.q_start)
                .then(a.q_end.cmp(&b.q_end))
                .then(a.t_start.cmp(&b.t_start))
                .then(a.t_end.cmp(&b.t_end))
        });

        // Pass 1: Remove Exact Matches
        let mut c_matched = vec![false; c_group.len()];
        let mut r_remaining = Vec::new();

        for r in r_group {
            let mut found = false;
            for (i, c) in c_group.iter().enumerate() {
                if !c_matched[i] && r == *c {
                    c_matched[i] = true;
                    found = true;
                    break;
                }
            }
            if !found {
                r_remaining.push(r);
            }
        }

        let c_remaining: Vec<&Rec> = c_group
            .into_iter()
            .enumerate()
            .filter(|(i, _)| !c_matched[*i])
            .map(|(_, r)| r)
            .collect();

        if r_remaining.is_empty() && c_remaining.is_empty() {
            continue;
        }

        report.push_str(&format!("[{}:{}]\n", q, t));

        let mut c_rem_matched = vec![false; c_remaining.len()];

        for r in &r_remaining {
            let mut matched_kind = "EXTRA (No matching coordinates in C)";
            let mut best_match_idx = None;

            for (i, c) in c_remaining.iter().enumerate() {
                if c_rem_matched[i] {
                    continue;
                }
                if r.coords_match(c) {
                    best_match_idx = Some(i);
                    matched_kind = "MISMATCH (Coordinates match, content differs)";
                    break;
                }
            }

            if let Some(idx) = best_match_idx {
                c_rem_matched[idx] = true;
                mismatch_count += 1;
                let c = c_remaining[idx];

                report.push_str(&format!("  {}\n", matched_kind));
                report.push_str(&format!(
                    "    Rust: coords={}-{}:{}-{} E={} FP={}\n",
                    r.q_start, r.q_end, r.t_start, r.t_end, r.energy, r.interaction
                ));
                report.push_str(&format!(
                    "    C:    coords={}-{}:{}-{} E={} FP={}\n",
                    c.q_start, c.q_end, c.t_start, c.t_end, c.energy, c.interaction
                ));

                if r.energy != c.energy {
                    report.push_str(&format!(
                        "    -> Energy Diff: Rust='{}' vs C='{}'\n",
                        r.energy, c.energy
                    ));
                }
                if r.interaction != c.interaction {
                    report.push_str(&format!(
                        "    -> Interaction Diff: Rust='{}' vs C='{}'\n",
                        r.interaction, c.interaction
                    ));
                }
                if r.target_seq != c.target_seq {
                    report.push_str(&format!(
                        "    -> TargetSeq Diff: Rust='{}' vs C='{}'\n",
                        r.target_seq, c.target_seq
                    ));
                }
                if r.flank_5 != c.flank_5 {
                    report.push_str(&format!(
                        "    -> Flank5 Diff: Rust='{}' vs C='{}'\n",
                        r.flank_5, c.flank_5
                    ));
                }
                if r.flank_3 != c.flank_3 {
                    report.push_str(&format!(
                        "    -> Flank3 Diff: Rust='{}' vs C='{}'\n",
                        r.flank_3, c.flank_3
                    ));
                }
            } else {
                extra_count += 1;
                report.push_str(&format!(
                    "  EXTRA IN RUST: coords={}-{}:{}-{} E={} FP={}\n",
                    r.q_start, r.q_end, r.t_start, r.t_end, r.energy, r.interaction
                ));
            }
        }

        for (i, c) in c_remaining.iter().enumerate() {
            if !c_rem_matched[i] {
                missing_count += 1;
                report.push_str(&format!(
                    "  MISSING IN RUST: coords={}-{}:{}-{} E={} FP={}\n",
                    c.q_start, c.q_end, c.t_start, c.t_end, c.energy, c.interaction
                ));
            }
        }
    }

    report.push_str("\nSUMMARY:\n");
    report.push_str(&format!(
        "  Total Mismatches (Coords matched): {}\n",
        mismatch_count
    ));
    report.push_str(&format!(
        "  Total Missing (In C, not Rust):    {}\n",
        missing_count
    ));
    report.push_str(&format!(
        "  Total Extra (In Rust, not C):      {}\n",
        extra_count
    ));

    if mismatch_count > 0 || missing_count > 0 || extra_count > 0 {
        panic!("{}", report);
    }
}

#[test]
fn parity_default_config() {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");
    let c_index = root.join("legacy_c/RIsearch2/test_suite/RHOC.pksuf");
    let c_bin = root.join("legacy_c/RIsearch2/bin/risearch2.x");
    let tmpdir = tempfile::tempdir().expect("tempdir");
    let rust_idx = tmpdir.path().join("RHOC.idx");

    let args = ["-l", "20", "-e", "-20", "-s", "6", "-p3"];

    let rust_out = index_and_search_rust(&query, &target, &rust_idx, &args);
    let c_out = search_c(&query, &c_index, &c_bin, &args);

    compare_results(&rust_out, &c_out, "default_config");
}

#[test]
fn parity_long_seed_no_ext() {
    // User requested: -l 0 -e 10000 -s 12 -p 3
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");
    let c_index = root.join("legacy_c/RIsearch2/test_suite/RHOC.pksuf");
    let c_bin = root.join("legacy_c/RIsearch2/bin/risearch2.x");
    let tmpdir = tempfile::tempdir().expect("tempdir");
    let rust_idx = tmpdir.path().join("RHOC.idx");

    let args = ["-l", "0", "-e", "10000", "-s", "12", "-p3"];

    let rust_out = index_and_search_rust(&query, &target, &rust_idx, &args);
    let c_out = search_c(&query, &c_index, &c_bin, &args);

    compare_results(&rust_out, &c_out, "long_seed_no_ext");
}
