use assert_cmd::cargo::cargo_bin_cmd;
use flate2::read::GzDecoder;
use std::collections::HashSet;
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
    query_seq: String,
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
            query_seq: if fields.len() > 12 {
                fields[12].to_string()
            } else {
                String::new()
            },
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

fn index_and_search_rust(
    query: &Path,
    target: &Path,
    index_path: &Path,
    args: &[&str],
) -> std::process::Output {
    // Build index
    let mut cmd = cargo_bin_cmd!("risearch");
    cmd.args([
        "index",
        target.to_str().unwrap(),
        index_path.to_str().unwrap(),
    ])
    .assert()
    .success();

    let mut cmd = cargo_bin_cmd!("risearch");
    // Search
    let mut final_args = vec![
        "search",
        "-i",
        index_path.to_str().unwrap(),
        "-q",
        query.to_str().unwrap(),
        "-o",
        "-", // Output to stdout
    ];
    final_args.extend_from_slice(args);

    let output = cmd.args(&final_args).output().expect("run risearch");

    if !output.status.success() {
        panic!("Rust search failed: {:?}", output.status);
    }
    output
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

fn setup_common_test_files() -> (tempfile::TempDir, PathBuf, PathBuf, PathBuf) {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");
    let c_bin = root.join("legacy_c/RIsearch2/bin/risearch2.x");
    let tmpdir = tempfile::tempdir().expect("tempdir");
    (tmpdir, query, target, c_bin)
}

fn create_c_index(target: &Path, c_index: &Path, c_bin: &Path) {
    let output = std::process::Command::new(c_bin)
        .arg("-c")
        .arg(target)
        .arg("-o")
        .arg(c_index)
        .output()
        .expect("create c index");
    if !output.status.success() {
        panic!("C index creation failed");
    }
}

fn compare_results(rust_out: &str, c_out: &str, test_name: &str) {
    // 1. Strict Byte-Level Parity Check (ignoring order and duplicates)
    let normalize = |s: &str| -> String {
        let mut lines: Vec<&str> = s.trim().split('\n').filter(|l| !l.is_empty()).collect();
        lines.sort();
        lines.dedup();
        lines.join("\n")
    };

    let r_norm = normalize(rust_out);
    let c_norm = normalize(c_out);

    if r_norm == c_norm {
        return;
    }

    // 2. Granular Reporting
    let rust_recs = parse_output(rust_out);
    let c_recs = parse_output(c_out);

    println!(
        "Compare Results Debug: Rust Recs: {}, C Recs: {}",
        rust_recs.len(),
        c_recs.len()
    );

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
                    "    Rust: coords={}-{}:{}-{} S={} E={}\n",
                    r.q_start, r.q_end, r.t_start, r.t_end, r.strand, r.energy
                ));
                if !r.query_seq.is_empty() {
                    report.push_str(&format!("          Query: {}\n", r.query_seq));
                }
                report.push_str(&format!("          FP:    {}\n", r.interaction));
                report.push_str(&format!("          Tgt:   {}\n", r.target_seq));
                report.push_str(&format!(
                    "    C:    coords={}-{}:{}-{} S={} E={}\n",
                    c.q_start, c.q_end, c.t_start, c.t_end, c.strand, c.energy
                ));
                if !c.query_seq.is_empty() {
                    report.push_str(&format!("          Query: {}\n", c.query_seq));
                }
                report.push_str(&format!("          FP:    {}\n", c.interaction));
                report.push_str(&format!("          Tgt:   {}\n", c.target_seq));

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
                    "  EXTRA IN RUST: coords={}-{}:{}-{} S={} E={}\n",
                    r.q_start, r.q_end, r.t_start, r.t_end, r.strand, r.energy
                ));
                if !r.query_seq.is_empty() {
                    report.push_str(&format!("                 Query: {}\n", r.query_seq));
                }
                report.push_str(&format!("                 FP:    {}\n", r.interaction));
                report.push_str(&format!("                 Tgt:   {}\n", r.target_seq));
            }
        }

        for (i, c) in c_remaining.iter().enumerate() {
            if !c_rem_matched[i] {
                missing_count += 1;
                report.push_str(&format!(
                    "  MISSING IN RUST: coords={}-{}:{}-{} S={} E={}\n",
                    c.q_start, c.q_end, c.t_start, c.t_end, c.strand, c.energy
                ));
                if !c.query_seq.is_empty() {
                    report.push_str(&format!("                   Query: {}\n", c.query_seq));
                }
                report.push_str(&format!("                   FP:    {}\n", c.interaction));
                report.push_str(&format!("                   Tgt:   {}\n", c.target_seq));
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
    let (tmpdir, query_path, target_path, c_bin) = setup_common_test_files();
    let c_index = tmpdir.path().join("c_target.pksuf");
    let rust_idx = tmpdir.path().join("rust_target.idx");

    create_c_index(&target_path, &c_index, &c_bin);

    let args = ["-l", "20", "-e", "-20", "-s", "6", "-p3"];

    let rust_output = index_and_search_rust(&query_path, &target_path, &rust_idx, &args);
    let rust_out = String::from_utf8_lossy(&rust_output.stdout).to_string();

    // Debug: Print extend_seed traces captured from stdout
    for line in rust_out.lines() {
        if line.contains("C_DEBUG:") {
            if line.contains("final=-20.49") || line.contains("final=-21.50") {
                println!("RUST: {}", line);
            }
        }
    }

    let c_out = search_c(&query_path, &c_index, &c_bin, &args);
    for line in c_out.lines() {
        if line.contains("C_DEBUG:") {
            if line.contains("final=-20.49") || line.contains("final=-21.50") {
                println!("C   : {}", line);
            }
        }
    }

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

    let rust_output = index_and_search_rust(&query, &target, &rust_idx, &args);
    let rust_out = String::from_utf8_lossy(&rust_output.stdout).to_string();
    let c_out = search_c(&query, &c_index, &c_bin, &args);

    compare_results(&rust_out, &c_out, "long_seed_no_ext");
}

#[test]
fn parity_energy_only() {
    let (tmpdir, query_path, target_path, c_bin) = setup_common_test_files();
    let c_index = tmpdir.path().join("c_target.pksuf");
    let rust_idx = tmpdir.path().join("rust_target.idx");

    create_c_index(&target_path, &c_index, &c_bin);

    // Test with just energy threshold, no detailed output (default output format)
    // Note: C implementation output format might differ slightly for default output
    // We'll use -p3 for consistent parsing but change other parameters
    let args = ["-l", "10", "-e", "-10.0", "-s", "5", "-p3"];

    let rust_output = index_and_search_rust(&query_path, &target_path, &rust_idx, &args);
    let rust_out = String::from_utf8_lossy(&rust_output.stdout).to_string();
    let c_out = search_c(&query_path, &c_index, &c_bin, &args);

    compare_results(&rust_out, &c_out, "energy_only");
}

#[test]
fn parity_reproduce_alignment_mismatch() {
    // User requested reproduction parameters:
    // target: TGGCTCTGTGGGACACAGCAGG
    // query: uggcucaguucagcaggaacag
    // args: -l 20 -e -20 -s 6 -p3

    let root = workspace_root();
    let tmpdir = tempfile::tempdir().expect("tempdir");

    // Create inputs
    let query_path = tmpdir.path().join("query.fa");
    let target_path = tmpdir.path().join("target.fa");

    fs::write(&query_path, ">query\nuggcucaguucagcaggaacag\n").expect("write query");
    fs::write(&target_path, ">target\nTGGCTCTGTGGGACACAGCAGG\n").expect("write target");

    let c_bin = root.join("legacy_c/RIsearch2/bin/risearch2.x");
    // We need to build a C index first
    let c_index = tmpdir.path().join("target.pksuf");

    // Index for C
    // Usage: risearch2.x -c <target_fasta> -o <output_index>
    let c_index_cmd = std::process::Command::new(&c_bin)
        .arg("-c")
        .arg(target_path.to_str().unwrap())
        .arg("-o")
        .arg(c_index.to_str().unwrap())
        .output()
        .expect("c index creation");

    if !c_index_cmd.status.success() {
        panic!("C indexing failed");
    }

    let rust_idx = tmpdir.path().join("target.idx");

    let args = ["-l", "20", "-e", "-20", "-s", "6", "-p3"];

    let rust_output = index_and_search_rust(&query_path, &target_path, &rust_idx, &args);
    let rust_out = String::from_utf8_lossy(&rust_output.stdout).to_string();
    let c_out = search_c(&query_path, &c_index, &c_bin, &args);

    println!("Rust Output:\n{}", rust_out);
    println!("C Output:\n{}", c_out);

    if rust_out != c_out {
        println!("Outputs differ!");
    } else {
        println!("Outputs are identical!");
    }

    // Force failure if inputs differ
    // We use compare_results to handle normalization (deduplication)
    // if rust_out.trim() != c_out.trim() { ... } -- REMOVED

    compare_results(&rust_out, &c_out, "alignment_mismatch_repro");
}
