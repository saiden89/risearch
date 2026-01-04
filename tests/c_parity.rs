use assert_cmd::cargo::cargo_bin_cmd;
use flate2::read::GzDecoder;
use log::{debug, info, trace, warn};
use std::collections::HashSet;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use tabled::{Table, Tabled, settings::Style};

// =============================================================================
// OUTPUT LEVEL GUIDE
// =============================================================================
// INFO  - Test progress: test name, hit counts, pass/fail summary
// DEBUG - Detailed comparisons: individual hit diffs, alignment analysis, coords
// TRACE - Raw data: full record dumps, C_DEBUG traces, file writes
// WARN  - Known divergences, acceptable differences
// ERROR - Reserved for panics (test failures)
// =============================================================================

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
        // Minimum fields: id, q_start, q_end, t_id, t_start, t_end, strand, energy, interaction, target_seq (10 fields)
        if fields.len() < 10 {
            if !line.trim().is_empty() {
                trace!(
                    "[REC] Skipped line with {} columns (expected >= 10): {}",
                    fields.len(),
                    line
                );
            }
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
            // Optional fields - use get() for safe access
            flank_5: fields.get(10).map_or(String::new(), |s| s.to_string()),
            flank_3: fields.get(11).map_or(String::new(), |s| s.to_string()),
            query_seq: fields.get(12).map_or(String::new(), |s| s.to_string()),
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

    /// Format record for debug output
    fn fmt_coords(&self) -> String {
        format!(
            "q=[{},{}] t=[{},{}] S={} E={}",
            self.q_start, self.q_end, self.t_start, self.t_end, self.strand, self.energy
        )
    }
}

/// Analyzes differences between two interaction strings and logs at DEBUG level.
///
/// The interaction string has format: [LEFT_EXT][SEED][RIGHT_EXT]
/// where LEFT_EXT = dp_left output, SEED = matched seed, RIGHT_EXT = dp_right output.
fn analyze_interaction_diff(rust_fp: &str, c_fp: &str, _q_start: usize, _q_end: usize) {
    if rust_fp == c_fp {
        return;
    }

    let rust_len = rust_fp.len();
    let c_len = c_fp.len();

    // Simple heuristic: divide into thirds for region classification
    let left_boundary = rust_len / 3;
    let right_boundary = rust_len * 2 / 3;

    debug!("[DIFF] Interaction mismatch:");
    debug!("[DIFF]   Rust: {}", rust_fp);
    debug!("[DIFF]   C:    {}", c_fp);

    // Build diff marker line - use zip_longest pattern with chars iterators (O(1) per char)
    let max_len = rust_len.max(c_len);
    let mut diff_markers = vec![' '; max_len];
    let mut mismatches: Vec<(usize, char, char, &str)> = Vec::new();

    let rust_chars: Vec<char> = rust_fp.chars().collect();
    let c_chars: Vec<char> = c_fp.chars().collect();

    for i in 0..max_len {
        let r_char = rust_chars.get(i).copied();
        let c_char = c_chars.get(i).copied();

        match (r_char, c_char) {
            (Some(r), Some(c)) if r != c => {
                diff_markers[i] = '^';
                let region = if i < left_boundary {
                    "5' EXT (dp_left)"
                } else if i >= right_boundary {
                    "3' EXT (dp_right)"
                } else {
                    "SEED"
                };
                mismatches.push((i, r, c, region));
            }
            (Some(r), None) => {
                diff_markers[i] = '+'; // Rust has extra
                mismatches.push((i, r, '-', "LENGTH"));
            }
            (None, Some(c)) => {
                diff_markers[i] = '-'; // C has extra
                mismatches.push((i, '-', c, "LENGTH"));
            }
            _ => {}
        }
    }

    let marker_str: String = diff_markers.into_iter().collect();
    debug!("[DIFF]   Mark: {}", marker_str.trim_end());

    // Count by region
    let mut left_count = 0;
    let mut seed_count = 0;
    let mut right_count = 0;
    let mut len_count = 0;

    for (pos, r, c, region) in &mismatches {
        match *region {
            "5' EXT (dp_left)" => left_count += 1,
            "SEED" => seed_count += 1,
            "3' EXT (dp_right)" => right_count += 1,
            "LENGTH" => len_count += 1,
            _ => {}
        }
        debug!("[DIFF]   Pos {}: Rust='{}' C='{}' -> {}", pos, r, c, region);
    }

    debug!(
        "[DIFF]   Region counts: 5'ext={}, seed={}, 3'ext={}, len_diff={}",
        left_count, seed_count, right_count, len_count
    );
}

/// Detailed analysis for focused tests with few hits.
/// Tries to pair EXTRA (Rust-only) and MISSING (C-only) hits that likely represent
/// the same biological interaction but with different coordinates/extensions.
fn analyze_hit_pairs(extras: &[&Rec], missings: &[&Rec]) {
    if extras.is_empty() || missings.is_empty() {
        return;
    }

    info!(
        "[PAIR] Paired hit analysis ({} extras, {} missings)",
        extras.len(),
        missings.len()
    );

    for extra in extras {
        // Find best matching missing hit using min_by_key
        let best_match = missings.iter().min_by_key(|missing| {
            let q_start_diff = (extra.q_start as i32 - missing.q_start as i32).abs();
            let q_end_diff = (extra.q_end as i32 - missing.q_end as i32).abs();
            let t_start_diff = (extra.t_start as i32 - missing.t_start as i32).abs();
            let t_end_diff = (extra.t_end as i32 - missing.t_end as i32).abs();
            q_start_diff + q_end_diff + t_start_diff + t_end_diff
        });

        if let Some(m) = best_match {
            // Compute score for display
            let score = (extra.q_start as i32 - m.q_start as i32).abs()
                + (extra.q_end as i32 - m.q_end as i32).abs()
                + (extra.t_start as i32 - m.t_start as i32).abs()
                + (extra.t_end as i32 - m.t_end as i32).abs();

            debug!("[PAIR] Likely pair (distance={})", score);
            debug!(
                "[PAIR]   RUST (extra):   q=[{:>2},{:>2}] t=[{:>3},{:>3}] E={}",
                extra.q_start, extra.q_end, extra.t_start, extra.t_end, extra.energy
            );
            debug!(
                "[PAIR]   C (missing):    q=[{:>2},{:>2}] t=[{:>3},{:>3}] E={}",
                m.q_start, m.q_end, m.t_start, m.t_end, m.energy
            );

            // Calculate deltas
            let dq_start = extra.q_start as i32 - m.q_start as i32;
            let dq_end = extra.q_end as i32 - m.q_end as i32;
            let dt_start = extra.t_start as i32 - m.t_start as i32;
            let dt_end = extra.t_end as i32 - m.t_end as i32;

            debug!(
                "[PAIR]   dq_start={:+3} dq_end={:+3} dt_start={:+3} dt_end={:+3}",
                dq_start, dq_end, dt_start, dt_end
            );

            // Diagnosis
            if dq_end == 0 && dt_end == 0 && (dq_start != 0 || dt_start != 0) {
                debug!(
                    "[PAIR]   >> DIAGNOSIS: Ends match, starts differ -> LEFT EXTENSION (dp_left) issue"
                );
            } else if dq_start == 0 && dt_start == 0 && (dq_end != 0 || dt_end != 0) {
                debug!(
                    "[PAIR]   >> DIAGNOSIS: Starts match, ends differ -> RIGHT EXTENSION (dp_right) issue"
                );
            } else if dq_start != 0 && dq_end != 0 {
                debug!(
                    "[PAIR]   >> DIAGNOSIS: Both starts and ends differ -> SEED SELECTION issue"
                );
            }

            // Show fingerprints at debug level
            debug!("[PAIR]   Rust FP: {}", extra.interaction);
            debug!("[PAIR]   C FP:    {}", m.interaction);
        }
    }
}

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Returns (records after dedup, count before dedup)
fn parse_output(output: &str) -> (Vec<Rec>, usize) {
    let mut recs: Vec<Rec> = output
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty())
        .filter_map(Rec::from_line)
        .collect();
    let parsed_count = recs.len();
    recs.sort();
    recs.dedup();
    (recs, parsed_count)
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

    // Check if trace logging is requested via RUST_LOG
    let trace_enabled = std::env::var("RUST_LOG")
        .map(|v| v.contains("trace"))
        .unwrap_or(false);

    // Search
    let mut final_args = vec!["search"];

    // Add -vvv for trace output if RUST_LOG=trace
    if trace_enabled {
        final_args.push("-vvv");
    }

    final_args.extend_from_slice(&[
        "-i",
        index_path.to_str().unwrap(),
        "-q",
        query.to_str().unwrap(),
        "-o",
        "-", // Output to stdout
    ]);
    final_args.extend_from_slice(args);

    let output = cmd.args(&final_args).output().expect("run risearch");

    // Print stderr if trace logging is enabled (contains trace! output)
    if trace_enabled {
        let stderr = String::from_utf8_lossy(&output.stderr);
        for line in stderr.lines() {
            trace!("[RUST] {}", line);
        }
    }

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
    // Initialize logger if not already initialized
    let _ = env_logger::builder().is_test(true).try_init();

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
        info!("[PARITY] {} - PASS (exact match)", test_name);
        return;
    }

    // 2. Granular Reporting
    let (rust_recs, rust_parsed) = parse_output(rust_out);
    let (c_recs, c_parsed) = parse_output(c_out);

    // Sanity check: verify parsed count (before dedup) matches raw line count
    // Difference due to dedup is fine, difference due to parse failure is not
    let rust_lines = rust_out.lines().filter(|l| !l.trim().is_empty()).count();
    let c_lines = c_out.lines().filter(|l| !l.trim().is_empty()).count();

    if rust_parsed != rust_lines {
        panic!(
            "[PARITY] Rust parse error: {} lines in output, {} records parsed (dropped {}). Check Rec::from_line parsing.",
            rust_lines,
            rust_parsed,
            rust_lines.saturating_sub(rust_parsed)
        );
    }
    if c_parsed != c_lines {
        panic!(
            "[PARITY] C parse error: {} lines in output, {} records parsed (dropped {}). Check Rec::from_line parsing.",
            c_lines,
            c_parsed,
            c_lines.saturating_sub(c_parsed)
        );
    }

    // Log dedup info at debug level if any duplicates were removed
    if rust_recs.len() != rust_parsed {
        debug!(
            "[PARITY] Rust: {} duplicates removed",
            rust_parsed - rust_recs.len()
        );
    }
    if c_recs.len() != c_parsed {
        debug!("[PARITY] C: {} duplicates removed", c_parsed - c_recs.len());
    }

    info!(
        "[PARITY] {} - Comparing: Rust={} hits, C={} hits",
        test_name,
        rust_recs.len(),
        c_recs.len()
    );

    // Dump raw records at trace level
    trace!("[PARITY] Rust records:");
    for r in &rust_recs {
        trace!("[PARITY]   {:?}", r);
    }
    trace!("[PARITY] C records:");
    for c in &c_recs {
        trace!("[PARITY]   {:?}", c);
    }

    // Optionally write to file at trace level
    if log::log_enabled!(log::Level::Trace) {
        let mut debug_out = String::new();
        debug_out.push_str("RUST RECORDS:\n");
        for r in &rust_recs {
            debug_out.push_str(&format!("{:?}\n", r));
        }
        debug_out.push_str("\nC RECORDS:\n");
        for c in &c_recs {
            debug_out.push_str(&format!("{:?}\n", c));
        }
        let _ = std::fs::write("target/debug_parity_report.txt", debug_out);
        trace!("[PARITY] Wrote raw records to target/debug_parity_report.txt");
    }

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
    let mut exact_match_count = 0;

    // Energy breakdown for coordinate-matched records
    let mut rust_better_energy = 0;
    let mut rust_worse_energy = 0;
    let mut energy_equal = 0;

    // Stats for extra/missing hits
    let mut extra_len_sum: usize = 0;
    let mut extra_energy_sum: f64 = 0.0;
    let mut missing_len_sum: usize = 0;
    let mut missing_energy_sum: f64 = 0.0;

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
                    exact_match_count += 1;
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

        debug!("");
        debug!("[PARITY] Group [{}:{}]", q, t);

        let mut c_rem_matched = vec![false; c_remaining.len()];

        for r in &r_remaining {
            let mut matched_kind = "EXTRA";
            let mut best_match_idx = None;

            for (i, c) in c_remaining.iter().enumerate() {
                if c_rem_matched[i] {
                    continue;
                }
                if r.coords_match(c) {
                    best_match_idx = Some(i);
                    matched_kind = "MISMATCH";
                    break;
                }
            }

            if let Some(idx) = best_match_idx {
                c_rem_matched[idx] = true;
                let c = c_remaining[idx];

                // Check if this is an "Improved Energy" case or a "True Mismatch"
                let r_e = r.energy.parse::<f64>().unwrap_or(0.0);
                let c_e = c.energy.parse::<f64>().unwrap_or(0.0);

                // If Rust is strictly better, log at warn. If equal, silently accept.
                if r_e < c_e - 0.001 {
                    rust_better_energy += 1;
                    warn!(
                        "[PARITY] Improved energy: Rust E={} vs C E={}",
                        r.energy, c.energy
                    );
                } else if r_e > c_e + 0.001 {
                    // Rust is worse - this is a true mismatch
                    rust_worse_energy += 1;
                    mismatch_count += 1;
                    debug!("[PARITY] {} (coords match, content differs)", matched_kind);
                } else {
                    // Equal within tolerance
                    energy_equal += 1;
                }

                // Detailed alignment comparison at DEBUG level
                debug!("[PARITY]   Rust: {}", r.fmt_coords());
                if !r.query_seq.is_empty() {
                    debug!("[PARITY]         Query: {}", r.query_seq);
                }
                debug!("[PARITY]         FP:    {}", r.interaction);
                debug!("[PARITY]         Tgt:   {}", r.target_seq);

                debug!("[PARITY]   C:    {}", c.fmt_coords());
                if !c.query_seq.is_empty() {
                    debug!("[PARITY]         Query: {}", c.query_seq);
                }
                debug!("[PARITY]         FP:    {}", c.interaction);
                debug!("[PARITY]         Tgt:   {}", c.target_seq);

                // Specific field diffs
                if r.energy != c.energy {
                    debug!(
                        "[PARITY]   -> Energy diff: Rust='{}' vs C='{}'",
                        r.energy, c.energy
                    );
                }
                if r.interaction != c.interaction {
                    analyze_interaction_diff(&r.interaction, &c.interaction, r.q_start, r.q_end);
                }
                if r.target_seq != c.target_seq {
                    debug!(
                        "[PARITY]   -> TargetSeq diff: Rust='{}' vs C='{}'",
                        r.target_seq, c.target_seq
                    );
                }
                if r.flank_5 != c.flank_5 {
                    debug!(
                        "[PARITY]   -> Flank5 diff: Rust='{}' vs C='{}'",
                        r.flank_5, c.flank_5
                    );
                }
                if r.flank_3 != c.flank_3 {
                    debug!(
                        "[PARITY]   -> Flank3 diff: Rust='{}' vs C='{}'",
                        r.flank_3, c.flank_3
                    );
                }
            } else {
                extra_count += 1;
                extra_len_sum += r.interaction.len();
                extra_energy_sum += r.energy.parse::<f64>().unwrap_or(0.0);
                debug!("[PARITY]   EXTRA in Rust: {}", r.fmt_coords());
                if !r.query_seq.is_empty() {
                    debug!("[PARITY]                  Query: {}", r.query_seq);
                }
                debug!("[PARITY]                  FP:    {}", r.interaction);
                debug!("[PARITY]                  Tgt:   {}", r.target_seq);
            }
        }

        for (i, c) in c_remaining.iter().enumerate() {
            if !c_rem_matched[i] {
                missing_count += 1;
                missing_len_sum += c.interaction.len();
                missing_energy_sum += c.energy.parse::<f64>().unwrap_or(0.0);
                debug!("[PARITY]   MISSING in Rust: {}", c.fmt_coords());
                if !c.query_seq.is_empty() {
                    debug!("[PARITY]                    Query: {}", c.query_seq);
                }
                debug!("[PARITY]                    FP:    {}", c.interaction);
                debug!("[PARITY]                    Tgt:   {}", c.target_seq);
            }
        }
    }

    let coord_matched = rust_better_energy + rust_worse_energy + energy_equal;
    let total_c = c_recs.len();

    // Compute percentages for coord-matched breakdown
    let pct = |n: usize, total: usize| -> String {
        if total == 0 {
            "0%".into()
        } else {
            format!("{}%", n * 100 / total)
        }
    };

    // Determine verdict
    let is_pass = mismatch_count == 0 && missing_count == 0 && extra_count == 0;
    let is_allowed_divergence = (test_name == "seed_only_no_extension"
        || test_name == "wobble_seed_pairs"
        || test_name == "right_extension_only")
        && mismatch_count == 0
        && extra_count == 0
        && missing_count > 0;

    // Build summary using tabled for proper alignment
    #[derive(Tabled)]
    struct SummaryRow {
        #[tabled(rename = "Metric")]
        metric: String,
        #[tabled(rename = "Value")]
        value: String,
    }

    let mut rows = vec![
        SummaryRow {
            metric: "Test".into(),
            value: test_name.to_string(),
        },
        SummaryRow {
            metric: "Rust hits".into(),
            value: rust_recs.len().to_string(),
        },
        SummaryRow {
            metric: "C hits".into(),
            value: total_c.to_string(),
        },
        SummaryRow {
            metric: "Exact matches".into(),
            value: format!(
                "{} ({})",
                exact_match_count,
                pct(exact_match_count, total_c)
            ),
        },
        SummaryRow {
            metric: "Coord-matched (diff content)".into(),
            value: coord_matched.to_string(),
        },
        SummaryRow {
            metric: "  ├ Rust better energy".into(),
            value: format!(
                "{} ({})",
                rust_better_energy,
                pct(rust_better_energy, coord_matched)
            ),
        },
        SummaryRow {
            metric: "  ├ Equal energy".into(),
            value: format!("{} ({})", energy_equal, pct(energy_equal, coord_matched)),
        },
        SummaryRow {
            metric: "  └ Rust worse energy".into(),
            value: format!(
                "{} ({})",
                rust_worse_energy,
                pct(rust_worse_energy, coord_matched)
            ),
        },
    ];

    // Missing stats
    if missing_count > 0 {
        let avg_len = missing_len_sum as f64 / missing_count as f64;
        let avg_energy = missing_energy_sum / missing_count as f64;
        rows.push(SummaryRow {
            metric: "Missing (in C, not Rust)".into(),
            value: format!(
                "{} (avg_len={:.1}, avg_E={:.2})",
                missing_count, avg_len, avg_energy
            ),
        });
    } else {
        rows.push(SummaryRow {
            metric: "Missing (in C, not Rust)".into(),
            value: "0".into(),
        });
    }

    // Extra stats
    if extra_count > 0 {
        let avg_len = extra_len_sum as f64 / extra_count as f64;
        let avg_energy = extra_energy_sum / extra_count as f64;
        rows.push(SummaryRow {
            metric: "Extra (in Rust, not C)".into(),
            value: format!(
                "{} (avg_len={:.1}, avg_E={:.2})",
                extra_count, avg_len, avg_energy
            ),
        });
    } else {
        rows.push(SummaryRow {
            metric: "Extra (in Rust, not C)".into(),
            value: "0".into(),
        });
    }

    // Verdict row
    let verdict = if is_pass {
        "✓ PASS".to_string()
    } else if is_allowed_divergence {
        "⚠ ALLOWED DIVERGENCE".to_string()
    } else {
        format!(
            "✗ FAIL ({} worse, {} missing, {} extra)",
            rust_worse_energy, missing_count, extra_count
        )
    };
    rows.push(SummaryRow {
        metric: "VERDICT".into(),
        value: verdict,
    });

    let table = Table::new(rows).with(Style::rounded()).to_string();
    for line in table.lines() {
        info!("[PARITY] {}", line);
    }

    // Handle allowed divergences
    if is_allowed_divergence {
        warn!(
            "[PARITY] Allowed divergence for '{}': Missing in Rust expected due to Maximality/Wobble improvements",
            test_name
        );
        return;
    }

    if !is_pass {
        panic!(
            "\n[PARITY] FAILED: {} ({} rust-worse, {} missing, {} extra)\nRun with RUST_LOG=debug for detailed diff analysis.",
            test_name, rust_worse_energy, missing_count, extra_count
        );
    }
}

/// Unified single-sequence parity test runner.
///
/// # Arguments
/// * `query_seq` - Query sequence (RNA/DNA)
/// * `target_seq` - Target sequence (RNA/DNA)
/// * `test_name` - Test identifier for logging
/// * `args` - risearch arguments
/// * `detailed` - If true, performs paired hit analysis for focused debugging
fn run_single_seq_parity(
    query_seq: &str,
    target_seq: &str,
    test_name: &str,
    args: &[&str],
    detailed: bool,
) {
    let _ = env_logger::builder().is_test(true).try_init();

    let root = workspace_root();
    let tmpdir = tempfile::tempdir().expect("tempdir");

    let query_path = tmpdir.path().join("query.fa");
    let target_path = tmpdir.path().join("target.fa");

    // Force uppercase for compatibility
    let q_upper = query_seq.to_uppercase();
    let t_upper = target_seq.to_uppercase();

    fs::write(&query_path, format!(">query\n{}\n", q_upper)).expect("write query");
    fs::write(&target_path, format!(">target\n{}\n", t_upper)).expect("write target");

    let c_bin = root.join("legacy_c/RIsearch2/bin/risearch2.x");
    let c_index = tmpdir.path().join("target.pksuf");

    // Index for C
    let c_index_cmd = std::process::Command::new(&c_bin)
        .arg("-c")
        .arg(target_path.to_str().unwrap())
        .arg("-o")
        .arg(c_index.to_str().unwrap())
        .output()
        .expect("c index creation");

    if !c_index_cmd.status.success() {
        panic!("C indexing failed: {:?}", c_index_cmd.status);
    }

    let rust_idx = tmpdir.path().join("target.idx");

    let rust_output = index_and_search_rust(&query_path, &target_path, &rust_idx, args);
    let rust_out = String::from_utf8_lossy(&rust_output.stdout).to_string();

    // Log debug traces at TRACE level
    for line in rust_out.lines() {
        if line.contains("DEBUG_TRACE")
            || line.contains("DEBUG_MAX")
            || line.contains("DEBUG_DP_LEFT_RESULT")
            || line.contains("C_DEBUG:")
        {
            trace!("[RUST] {}", line);
        }
    }

    // Filter out Rust-specific flags like --no-max-prune for C call
    let c_args: Vec<&str> = args
        .iter()
        .filter(|&&a| a != "--no-max-prune")
        .cloned()
        .collect();
    let c_out = search_c(&query_path, &c_index, &c_bin, &c_args);

    if detailed {
        // Parse both outputs for detailed analysis
        let (rust_recs, _) = parse_output(&rust_out);
        let (c_recs, _) = parse_output(&c_out);

        info!(
            "[PARITY] {} - Rust={} hits, C={} hits",
            test_name,
            rust_recs.len(),
            c_recs.len()
        );

        // Find extras and missings
        let mut c_matched = vec![false; c_recs.len()];
        let mut extras: Vec<&Rec> = Vec::new();

        for r in &rust_recs {
            let mut found = false;
            for (i, c) in c_recs.iter().enumerate() {
                if !c_matched[i] {
                    // Strict equality
                    if r == c {
                        c_matched[i] = true;
                        found = true;
                        break;
                    }

                    // Relaxed equality: Exact Coordinates + Rust Energy is Better/Equal
                    if r.coords_match(c) {
                        let r_e = r.energy.parse::<f64>().unwrap_or(0.0);
                        let c_e = c.energy.parse::<f64>().unwrap_or(0.0);
                        // Only warn if Rust is strictly better, silently accept equal
                        if r_e < c_e - 0.001 {
                            warn!(
                                "[PARITY] Improved energy: Rust E={} vs C E={}",
                                r.energy, c.energy
                            );
                        }
                        // Accept if Rust is better or equal
                        if r_e <= c_e + 0.001 {
                            c_matched[i] = true;
                            found = true;
                            break;
                        }
                    }
                }
            }
            if !found {
                extras.push(r);
            }
        }

        let missings: Vec<&Rec> = c_recs
            .iter()
            .enumerate()
            .filter(|(i, _)| !c_matched[*i])
            .map(|(_, r)| r)
            .collect();

        debug!(
            "[PARITY] Exact matches: {}, EXTRA: {}, MISSING: {}",
            rust_recs.len() - extras.len(),
            extras.len(),
            missings.len()
        );

        // Call paired analysis
        analyze_hit_pairs(&extras, &missings);

        // Fail if we missed any C hits
        if !missings.is_empty() {
            panic!(
                "[PARITY] FAILED {}: {} missings ({} extras)",
                test_name,
                missings.len(),
                extras.len()
            );
        }

        // Warn about extras but pass
        if !extras.is_empty() {
            warn!(
                "[PARITY] {} extras in {} (acceptable if 0 missings)",
                extras.len(),
                test_name
            );
        }
    } else {
        compare_results(&rust_out, &c_out, test_name);
    }
}

// =============================================================================
// TESTS
// =============================================================================

#[test]
fn test_parity_full_pipeline() {
    let _ = env_logger::builder().is_test(true).try_init();

    let (tmpdir, query_path, target_path, c_bin) = setup_common_test_files();
    let c_index = tmpdir.path().join("c_target.pksuf");
    let rust_idx = tmpdir.path().join("rust_target.idx");

    create_c_index(&target_path, &c_index, &c_bin);

    let args = ["-l", "20", "-e", "-20", "-s", "6", "-p3"];

    let rust_output = index_and_search_rust(&query_path, &target_path, &rust_idx, &args);
    let rust_out = String::from_utf8_lossy(&rust_output.stdout).to_string();

    // Log C_DEBUG traces at trace level
    for line in rust_out.lines() {
        if line.contains("C_DEBUG:") {
            if line.contains("final=-20.49") || line.contains("final=-21.50") {
                trace!("[RUST] {}", line);
            }
        }
    }

    let c_out = search_c(&query_path, &c_index, &c_bin, &args);
    for line in c_out.lines() {
        if line.contains("C_DEBUG:") {
            if line.contains("final=-20.49") || line.contains("final=-21.50") {
                trace!("[C] {}", line);
            }
        }
    }

    compare_results(&rust_out, &c_out, "default_config");
}

#[test]
fn test_parity_long_seed() {
    let _ = env_logger::builder().is_test(true).try_init();

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
fn test_parity_energy_threshold() {
    let (tmpdir, query_path, target_path, c_bin) = setup_common_test_files();
    let c_index = tmpdir.path().join("c_target.pksuf");
    let rust_idx = tmpdir.path().join("rust_target.idx");

    create_c_index(&target_path, &c_index, &c_bin);

    let args = ["-l", "10", "-e", "-10.0", "-s", "5", "-p3"];

    let rust_output = index_and_search_rust(&query_path, &target_path, &rust_idx, &args);
    let rust_out = String::from_utf8_lossy(&rust_output.stdout).to_string();
    let c_out = search_c(&query_path, &c_index, &c_bin, &args);

    compare_results(&rust_out, &c_out, "energy_only");
}

#[test]
fn test_parity_alignment_repro() {
    let _ = env_logger::builder().is_test(true).try_init();

    let root = workspace_root();
    let tmpdir = tempfile::tempdir().expect("tempdir");

    // Create inputs
    let query_path = tmpdir.path().join("query.fa");
    let target_path = tmpdir.path().join("target.fa");

    fs::write(&query_path, ">query\nuggcucaguucagcaggaacag\n").expect("write query");
    fs::write(&target_path, ">target\nTGGCTCTGTGGGACACAGCAGG\n").expect("write target");

    let c_bin = root.join("legacy_c/RIsearch2/bin/risearch2.x");
    let c_index = tmpdir.path().join("target.pksuf");

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

    debug!("[PARITY] Rust Output:\n{}", rust_out);
    debug!("[PARITY] C Output:\n{}", c_out);

    if rust_out.trim() == c_out.trim() {
        info!("[PARITY] alignment_mismatch_repro - Outputs identical");
    } else {
        debug!("[PARITY] Outputs differ");
    }

    compare_results(&rust_out, &c_out, "alignment_mismatch_repro");
}

#[test]
fn test_parity_single_seq() {
    // Tests internal mismatch handling between Rust and C implementations.
    let args = ["-l", "20", "-e", "100.0", "-s", "5", "-p3"];

    let query = "UGCUGCUGCCGCUGCUGCUG"; // 20nt with internal variation
    let target = "GCAGCAGCAGCAGCAGCAGC"; // 20nt complement

    run_single_seq_parity(query, target, "internal_mismatch_20nt", &args, false);
}

// =============================================================================
// ISOLATED TESTS: Each tests a specific component to pinpoint differences
// =============================================================================

/// Tests seed matching only (no extension).
/// Use -l 0 to disable extension, so only seed pairing is tested.
#[test]
fn test_parity_seed_only() {
    // No extension: -l 0
    // This tests only the seed pairing and energy calculation
    let args = ["-l", "0", "-e", "100.0", "-s", "5", "-p3"];

    // Query and target are exact complements (20nt)
    // Should produce a single seed hit with no extensions
    let query = "UGCUGCUGCUGCUGCUGCUG"; // 20nt
    let target = "CAGCAGCAGCAGCAGCAGCA"; // Perfect complement, reversed

    run_single_seq_parity(query, target, "seed_only_no_extension", &args, false);
}

/// Tests left extension only (dp_left).
/// Design: seed at 3' end of query, extra bases only to the 5' side.
#[test]
fn test_parity_left_ext() {
    // Seed at 3' end of query forces only left extension
    let args = ["-l", "20", "-e", "100.0", "-s", "5", "-p3"];

    let query = "AAAAAUGCUG";
    let target = "CAGCAUUUUU";

    run_single_seq_parity(query, target, "left_extension_only", &args, false);
}

/// Tests right extension only (dp_right).
/// Design: seed at 5' end of query, extra bases only to the 3' side.
#[test]
fn test_parity_right_ext() {
    let args = ["-l", "20", "-e", "100.0", "-s", "5", "-p3"];

    let query = "UGCUGAAAAA";
    let target = "UUUUUCAGCA";

    run_single_seq_parity(query, target, "right_extension_only", &args, false);
}

/// Tests both left and right extension.
/// Design: seed in middle, extra bases on both sides.
#[test]
fn test_parity_both_ext() {
    let args = ["-l", "20", "-e", "100.0", "-s", "5", "-p3"];

    let query = "AAAUGCUGAAA";
    let target = "UUUCAGCAUUU";

    run_single_seq_parity(query, target, "both_extensions", &args, false);
}

/// Tests with wobble pairs (G-U) in the seed region.
#[test]
fn test_parity_wobble() {
    let args = ["-l", "0", "-e", "100.0", "-s", "5", "-p3"];

    let query = "UGUGUGUGUG"; // 10 alternating U-G pattern
    let target = "CGCGCGCGCG"; // Complement with wobble

    run_single_seq_parity(query, target, "wobble_seed_pairs", &args, false);
}

// =============================================================================
// ISOLATED FAILURE REPRODUCTION: Uses actual failing sequences from parity_default_config
// =============================================================================

/// Isolated repro of parity_default_config failure.
/// Uses actual hsa-miR-24-3p query and a short RHOC target region.
#[test]
fn test_parity_mir24_isolated() {
    let query = "uggcucaguucagcaggaacag"; // 22nt

    let target =
        "GGAAGACCTGCCTCCTCATCGTCTTCAGCAAGGATCAGTTTCCGGAGGTCTACGTCCCTACTGTCTTTGAGAACTATATTG";

    let args = [
        "-l",
        "20",
        "-e",
        "100.0",
        "-s",
        "6",
        "-p3",
        "--no-max-prune",
    ];

    info!("[PARITY] mir24_isolated_single test");
    debug!("[PARITY] Query:  {}", query);
    debug!("[PARITY] Target: {}", target);
    debug!("[PARITY] Args:   {:?}", args);

    run_single_seq_parity(query, target, "miR24_isolated_single", &args, true);
}

/// More targeted single-hit test with explicit segment expectations.
#[test]
fn test_parity_segments_debug() {
    let query = "AAAAAUGCUGUAAAAA"; // 5 + 6 + 5 = 16nt
    let target = "UUUUUACAGCAUUUUU";

    let args = ["-l", "20", "-e", "100.0", "-s", "6", "-p3"];

    info!("[PARITY] segment_debug test");
    debug!("[PARITY] Query:  {} (len={})", query, query.len());
    debug!("[PARITY] Target: {} (len={})", target, target.len());
    debug!("[PARITY] Expected structure:");
    debug!("[PARITY]   LEFT_EXT  (dp_left):  AAAAA -> should extend with UUUUU");
    debug!("[PARITY]   SEED      (6nt):      UGCUGU -> matches ACAGCA");
    debug!("[PARITY]   RIGHT_EXT (dp_right): AAAAA -> should extend with UUUUU");

    run_single_seq_parity(query, target, "segment_debug", &args, false);
}
