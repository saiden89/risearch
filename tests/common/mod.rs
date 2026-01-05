use assert_cmd::cargo::cargo_bin_cmd;
use flate2::read::GzDecoder;
use log::{debug, info, trace, warn};
use std::collections::HashSet;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::Once;
use tabled::{Table, Tabled, builder::Builder, settings::Style};

static INIT: Once = Once::new();

pub fn init_test_logging() {
    INIT.call_once(|| {
        let _ = env_logger::builder()
            .is_test(true)
            .format(|buf, record| {
                use std::io::Write;
                let style = buf.default_level_style(record.level());
                writeln!(
                    buf,
                    "[{}{}{}] {}",
                    style.render(),
                    record.level(),
                    style.render_reset(),
                    record.args()
                )
            })
            .try_init();
    });
}

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
pub struct Rec {
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
                    "{} Skipped line with {} columns (expected >= 10): {}",
                    LogTag::Rec,
                    fields.len(),
                    line
                );
            }
            return None;
        }

        // Helper to strip y/x debug markers from C output
        let strip = |s: &str| -> String {
            let before = char::from(SeedMarker::BeforeSeed);
            let after = char::from(SeedMarker::AfterSeed);
            s.chars().filter(|&c| c != before && c != after).collect()
        };

        Some(Rec {
            q_id: fields[0].to_string(),
            q_start: fields[1].parse().unwrap_or(0),
            q_end: fields[2].parse().unwrap_or(0),
            t_id: fields[3].to_string(),
            t_start: fields[4].parse().unwrap_or(0),
            t_end: fields[5].parse().unwrap_or(0),
            strand: fields[6].to_string(),
            energy: fields[7].to_string(),
            interaction: strip(fields[8]),
            target_seq: strip(fields[9]),
            // Optional fields - use get() for safe access
            flank_5: fields.get(10).map_or(String::new(), |s| strip(s)),
            flank_3: fields.get(11).map_or(String::new(), |s| strip(s)),
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

// =============================================================================
// PARITY DATA MODELS & VIEW
// =============================================================================

#[derive(Debug)]
pub enum ParityKind<'a> {
    Mismatch { rust: &'a Rec, c: &'a Rec },
    RustOnly(&'a Rec),
    COnly(&'a Rec),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ColumnKind {
    Label,
    Ctx5,
    Ext5,
    Seed,
    Ext3,
    Ctx3,
}

impl ColumnKind {
    pub fn header(&self) -> &'static str {
        match self {
            Self::Label => "", // Placeholder
            Self::Ctx5 => "5' Ctx",
            Self::Ext5 => "5' EXT",
            Self::Seed => "SEED",
            Self::Ext3 => "3' EXT",
            Self::Ctx3 => "3' Ctx",
        }
    }
}

pub struct TableConfig {
    pub columns: Vec<ColumnKind>,
}

impl Default for TableConfig {
    fn default() -> Self {
        Self {
            columns: vec![
                ColumnKind::Label,
                ColumnKind::Ctx5,
                ColumnKind::Ext5,
                ColumnKind::Seed,
                ColumnKind::Ext3,
                ColumnKind::Ctx3,
            ],
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum RowLabel {
    SingleFP,
    SingleTarget,
    CompCFP,
    CompCTarget,
    CompDiff,
    CompRFP,
    CompRTarget,
}

impl std::fmt::Display for RowLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SingleFP => write!(f, "FP"),
            Self::SingleTarget => write!(f, "Target"),
            Self::CompCFP => write!(f, "C FP"),
            Self::CompCTarget => write!(f, "C Tgt"),
            Self::CompDiff => write!(f, "DIFF"),
            Self::CompRFP => write!(f, "R FP"),
            Self::CompRTarget => write!(f, "R Tgt"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiffChar {
    Match,    // ' '
    Mismatch, // 'X'
}

impl DiffChar {
    pub fn as_char(&self) -> char {
        match self {
            Self::Match => ' ',
            Self::Mismatch => 'X',
        }
    }
}

impl From<(char, char)> for DiffChar {
    fn from((a, b): (char, char)) -> Self {
        if a == b { Self::Match } else { Self::Mismatch }
    }
}

pub enum SeedMarker {
    BeforeSeed, // 'y'
    AfterSeed,  // 'x'
}

impl From<SeedMarker> for char {
    fn from(marker: SeedMarker) -> char {
        match marker {
            SeedMarker::BeforeSeed => 'y',
            SeedMarker::AfterSeed => 'x',
        }
    }
}

pub struct ParityTable<'a> {
    pub kind: ParityKind<'a>,
    pub config: TableConfig,
}

pub enum LogTag {
    Rec,
    Pair,
    Rust,
    C,
    Parity,
    Summary,
}

impl std::fmt::Display for LogTag {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Rec => write!(f, "[REC]"),
            Self::Pair => write!(f, "[PAIR]"),
            Self::Rust => write!(f, "[RUST]"),
            Self::C => write!(f, "[C]"),
            Self::Parity => write!(f, "[PARITY]"),
            Self::Summary => write!(f, "[SUMMARY]"),
        }
    }
}

/// Status of a hit in parity comparison
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HitStatus {
    Identical,  // Exact match between Rust and C (markers already stripped)
    CoOptimal,  // Equal energy, different trace (co-optimal)
    RustBetter, // Rust has better energy
    RustWorse,  // Rust has worse energy
    Extra,      // Only in Rust
    Missing,    // Only in C
}

impl std::fmt::Display for HitStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Identical => write!(f, "IDENTICAL"),
            Self::CoOptimal => write!(f, "CO-OPTIMAL"),
            Self::RustBetter => write!(f, "RUST BETTER"),
            Self::RustWorse => write!(f, "RUST WORSE"),
            Self::Extra => write!(f, "EXTRA"),
            Self::Missing => write!(f, "MISSING"),
        }
    }
}

impl HitStatus {
    pub fn is_significant(&self) -> bool {
        !matches!(self, Self::Identical)
    }
}

pub struct ParsedInteraction {
    pub ctx_5: String,
    pub ext_5: String,
    pub seed: String,
    pub ext_3: String,
    pub ctx_3: String,
}

fn build_diff(a: &str, b: &str) -> String {
    let a_chars: Vec<char> = a.chars().collect();
    let b_chars: Vec<char> = b.chars().collect();
    let len = a_chars.len().max(b_chars.len());
    let mut diff = String::with_capacity(len);
    for i in 0..len {
        let ac = a_chars.get(i).copied().unwrap_or(' ');
        let bc = b_chars.get(i).copied().unwrap_or(' ');
        diff.push(DiffChar::from((ac, bc)).as_char());
    }
    diff
}

impl ParsedInteraction {
    fn from_rec(rec: &Rec, strip: bool, guide_markers: Option<(usize, usize)>) -> Self {
        let fp_raw = &rec.interaction;

        let (left, seed, right) = if let Some((start, end)) = guide_markers {
            let chars: Vec<char> = fp_raw.chars().collect();
            let len = chars.len();
            if len >= end {
                (
                    chars[..start].iter().collect(),
                    chars[start..end].iter().collect(),
                    chars[end..].iter().collect(),
                )
            } else {
                ("".into(), fp_raw.clone(), "".into())
            }
        } else if strip {
            if let Some((clean, start, end)) = parse_seed_markers(fp_raw) {
                let chars: Vec<char> = clean.chars().collect();
                (
                    chars[..start].iter().collect(),
                    chars[start..end].iter().collect(),
                    chars[end..].iter().collect(),
                )
            } else {
                ("".into(), strip_markers(fp_raw), "".into())
            }
        } else {
            ("".into(), fp_raw.clone(), "".into())
        };

        let ctx_5 = if strip {
            strip_markers(&rec.flank_5)
        } else {
            rec.flank_5.clone()
        };
        let ctx_3 = if strip {
            strip_markers(&rec.flank_3)
        } else {
            rec.flank_3.clone()
        };

        Self {
            ctx_5,
            ext_5: left,
            seed,
            ext_3: right,
            ctx_3,
        }
    }

    fn from_rec_target(rec: &Rec, strip: bool, ref_parts: &ParsedInteraction) -> Self {
        let tgt_raw = if strip {
            strip_markers(&rec.target_seq)
        } else {
            rec.target_seq.clone()
        };
        let chars: Vec<char> = tgt_raw.chars().collect();

        let l_len = ref_parts.ext_5.chars().count();
        let s_len = ref_parts.seed.chars().count();

        let total = chars.len();
        let p1 = l_len.min(total);
        let p2 = (l_len + s_len).min(total);

        Self {
            ctx_5: ref_parts.ctx_5.clone(),
            ext_5: chars[..p1].iter().collect(),
            seed: chars[p1..p2].iter().collect(),
            ext_3: chars[p2..].iter().collect(),
            ctx_3: ref_parts.ctx_3.clone(),
        }
    }
}

impl<'a> std::fmt::Display for ParityTable<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut builder = Builder::default();

        let strand = match &self.kind {
            ParityKind::Mismatch { c, .. } => &c.strand,
            ParityKind::RustOnly(r) => &r.strand,
            ParityKind::COnly(c) => &c.strand,
        };

        let headers: Vec<String> = self
            .config
            .columns
            .iter()
            .map(|c| {
                if matches!(c, ColumnKind::Label) {
                    format!("S={}", strand)
                } else {
                    c.header().to_string()
                }
            })
            .collect();
        builder.push_record(headers);

        let add_row = |b: &mut Builder, label: RowLabel, p: &ParsedInteraction| {
            let label_str = label.to_string();
            let row: Vec<String> = self
                .config
                .columns
                .iter()
                .map(|col| match col {
                    ColumnKind::Label => label_str.clone(),
                    ColumnKind::Ctx5 => p.ctx_5.clone(),
                    ColumnKind::Ext5 => p.ext_5.clone(),
                    ColumnKind::Seed => p.seed.clone(),
                    ColumnKind::Ext3 => p.ext_3.clone(),
                    ColumnKind::Ctx3 => p.ctx_3.clone(),
                })
                .collect();
            b.push_record(row);
        };

        match self.kind {
            ParityKind::RustOnly(r) => {
                let p = ParsedInteraction::from_rec(r, false, None);
                let p_tgt = ParsedInteraction::from_rec_target(r, false, &p);
                add_row(&mut builder, RowLabel::SingleTarget, &p_tgt);
                add_row(&mut builder, RowLabel::SingleFP, &p);
            }
            ParityKind::COnly(c) => {
                let p = ParsedInteraction::from_rec(c, true, None);
                let p_tgt = ParsedInteraction::from_rec_target(c, true, &p);
                add_row(&mut builder, RowLabel::SingleTarget, &p_tgt);
                add_row(&mut builder, RowLabel::SingleFP, &p);
            }
            ParityKind::Mismatch { rust: r, c } => {
                let p_c = ParsedInteraction::from_rec(c, true, None);

                let guide = if let Some((_, start, end)) = parse_seed_markers(&c.interaction) {
                    Some((start, end))
                } else {
                    None
                };
                let p_r = ParsedInteraction::from_rec(r, false, guide);

                let diff_l = build_diff(&p_c.ext_5, &p_r.ext_5);
                let diff_s = build_diff(&p_c.seed, &p_r.seed);
                let diff_r = build_diff(&p_c.ext_3, &p_r.ext_3);

                let p_diff = ParsedInteraction {
                    ctx_5: "".into(),
                    ext_5: diff_l,
                    seed: diff_s,
                    ext_3: diff_r,
                    ctx_3: "".into(),
                };

                let p_c_tgt = ParsedInteraction::from_rec_target(c, true, &p_c);
                let p_r_tgt = ParsedInteraction::from_rec_target(r, false, &p_r);

                add_row(&mut builder, RowLabel::CompCTarget, &p_c_tgt);
                add_row(&mut builder, RowLabel::CompCFP, &p_c);
                add_row(&mut builder, RowLabel::CompDiff, &p_diff);
                add_row(&mut builder, RowLabel::CompRFP, &p_r);
                add_row(&mut builder, RowLabel::CompRTarget, &p_r_tgt);
            }
        }

        write!(f, "{}", builder.build().with(Style::rounded()))
    }
}

// For debug output comparison, we insert y/x markers around the seed portion
// of interaction strings. Format: <left_ext>y<seed>x<right_ext>
// =============================================================================

/// Parse y/x markers from a C debug output string.
/// Returns (interaction_without_markers, seed_start, seed_end) or None if no markers.
pub fn parse_seed_markers(s: &str) -> Option<(String, usize, usize)> {
    let before = char::from(SeedMarker::BeforeSeed);
    let after = char::from(SeedMarker::AfterSeed);

    let y_pos = s.find(before)?;
    let x_pos = s.find(after)?;

    if x_pos <= y_pos {
        return None; // Invalid marker order
    }

    // Remove markers and compute seed indices in the clean string
    let clean: String = s.chars().filter(|&c| c != before && c != after).collect();
    let seed_start = y_pos;
    let seed_end = x_pos - 1; // -1 because y was before this position

    Some((clean, seed_start, seed_end))
}

/// Strip y/x markers from a string (simple cleanup without parsing positions)
/// Strip y/x markers from a string (simple cleanup without parsing positions)
fn strip_markers(s: &str) -> String {
    let before = char::from(SeedMarker::BeforeSeed);
    let after = char::from(SeedMarker::AfterSeed);
    s.chars().filter(|&c| c != before && c != after).collect()
}

/// Format a Rec as a horizontal table for display
/// Columns: 5' Context | 5' EXT | SEED | 3' EXT | 3' Context
/// Rows: FP, Target (and optionally Query)

/// Analyzes differences between Rust and C interaction strings.
///
/// If C output contains y/x seed markers, shows aligned comparison:
/// ```
/// C:    [left_ext] y[seed]x [right_ext]
/// Rust: [left_ext]  [seed]  [right_ext]
/// ```

/// Detailed analysis for focused tests with few hits.
/// Tries to pair EXTRA (Rust-only) and MISSING (C-only) hits that likely represent
/// the same biological interaction but with different coordinates/extensions.
pub fn analyze_hit_pairs(extras: &[&Rec], missings: &[&Rec]) {
    if extras.is_empty() || missings.is_empty() {
        return;
    }

    info!(
        "{} Paired hit analysis ({} extras, {} missings)",
        LogTag::Pair,
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

            debug!("{} Likely pair (distance={})", LogTag::Pair, score);

            if extra.interaction == m.interaction {
                debug!("{}   [NOTE] Interactions Identical!", LogTag::Pair);
                if extra.energy != m.energy {
                    debug!(
                        "{}   Energy Diff: Rust={} vs C={}",
                        LogTag::Pair,
                        extra.energy,
                        m.energy
                    );
                }
                if extra.target_seq != m.target_seq {
                    debug!(
                        "{}   Target Diff: Rust={} vs C={}",
                        LogTag::Pair,
                        extra.target_seq,
                        m.target_seq
                    );
                }
                // Check coords
                if !extra.coords_match(m) {
                    debug!(
                        "{}   Coords Diff: Rust={} vs C={}",
                        LogTag::Pair,
                        extra.fmt_coords(),
                        m.fmt_coords()
                    );
                }
            }

            // Show as a mismatch table to visualize the alignment differences
            let table = ParityTable {
                kind: ParityKind::Mismatch { rust: extra, c: m },
                config: TableConfig::default(),
            };

            for line in table.to_string().lines() {
                debug!("{}   {}", LogTag::Pair, line);
            }
        }
    }
}

pub fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Returns (records after dedup, count before dedup)
pub fn parse_output(output: &str) -> (Vec<Rec>, usize) {
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

pub fn index_and_search_rust(
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

    // Pass through binary stderr directly
    // Note: won't be colored since subprocess output is captured, not a TTY
    if trace_enabled {
        eprint!("{}", String::from_utf8_lossy(&output.stderr));
    }

    if !output.status.success() {
        panic!("Rust search failed: {:?}", output.status);
    }
    output
}

pub fn search_c(query: &Path, c_index: &Path, c_bin: &Path, args: &[&str]) -> String {
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

pub fn setup_common_test_files() -> (tempfile::TempDir, PathBuf, PathBuf, PathBuf) {
    let root = workspace_root();
    let query = root.join("legacy_c/RIsearch2/test_suite/mirnas.fa");
    let target = root.join("legacy_c/RIsearch2/test_suite/RHOC.fa");
    let c_bin = c_binary_path(&root);
    let tmpdir = tempfile::tempdir().expect("tempdir");
    (tmpdir, query, target, c_bin)
}

/// Get path to C risearch2 binary.
/// Prefers debug binary (risearch2.dbg.x) if it exists, otherwise uses release (risearch2.x).
/// Debug binary has x/y seed markers and DP matrix output for parity debugging.
pub fn c_binary_path(root: &Path) -> PathBuf {
    let debug_bin = root.join("legacy_c/RIsearch2/bin/risearch2.dbg.x");
    let release_bin = root.join("legacy_c/RIsearch2/bin/risearch2.x");

    if debug_bin.exists() {
        debug!(
            "{} Using debug C binary: {}",
            LogTag::Parity,
            debug_bin.display()
        );
        debug_bin
    } else {
        debug!(
            "{} Using release C binary: {}",
            LogTag::Parity,
            release_bin.display()
        );
        release_bin
    }
}

pub fn create_c_index(target: &Path, c_index: &Path, c_bin: &Path) {
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

pub fn compare_results(rust_out: &str, c_out: &str, test_name: &str) {
    // Initialize logger if not already initialized
    // Initialize logger if not already initialized
    init_test_logging();

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
        info!("{} {} - PASS (exact match)", LogTag::Parity, test_name);
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
        trace!(
            "{} Wrote raw records to target/debug_parity_report.txt",
            LogTag::Parity
        );
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
        debug!("{} Group [{}:{}]", LogTag::Parity, q, t);

        let mut c_rem_matched = vec![false; c_remaining.len()];

        for r in &r_remaining {
            let mut best_match_idx = None;

            for (i, c) in c_remaining.iter().enumerate() {
                if c_rem_matched[i] {
                    continue;
                }
                if r.coords_match(c) {
                    best_match_idx = Some(i);
                    break;
                }
            }

            if let Some(idx) = best_match_idx {
                c_rem_matched[idx] = true;
                let c = c_remaining[idx];

                // Check if this is an "Improved Energy" case or a "True Mismatch"
                let r_e = r.energy.parse::<f64>().unwrap_or(0.0);
                let c_e = c.energy.parse::<f64>().unwrap_or(0.0);

                let status: HitStatus;

                // If Rust is strictly better, log at warn. If equal, silently accept.
                if r_e < c_e - 0.001 {
                    rust_better_energy += 1;
                    status = HitStatus::RustBetter;
                    warn!(
                        "{} Improved energy: Rust E={} vs C E={}",
                        LogTag::Parity,
                        r.energy,
                        c.energy
                    );
                } else if r_e > c_e + 0.001 {
                    // Rust is worse - this is a true mismatch
                    rust_worse_energy += 1;
                    mismatch_count += 1;
                    status = HitStatus::RustWorse;
                    debug!(
                        "{} {} (coords match, content differs)",
                        LogTag::Parity,
                        status
                    );
                } else {
                    // Equal within tolerance - markers already stripped at parse time
                    energy_equal += 1;
                    if r.interaction == c.interaction {
                        status = HitStatus::Identical;
                    } else {
                        status = HitStatus::CoOptimal;
                    }
                }

                // Use ParityTable for all mismatch display (it handles parsing/diff internally)
                let table = ParityTable {
                    kind: ParityKind::Mismatch { rust: r, c },
                    config: TableConfig::default(),
                };

                // Only show significant differences
                if status.is_significant() {
                    debug!("{}", LogTag::Parity);
                    let energy_msg = if r.energy == c.energy {
                        "".to_string()
                    } else {
                        format!(" | C Energy: {}", c.energy)
                    };
                    debug!(
                        "{}   {} Coords: {}{}",
                        LogTag::Parity,
                        status,
                        r.fmt_coords(),
                        energy_msg
                    );

                    // Print table as single block to preserve alignment
                    debug!("{}   {}\n{}", LogTag::Parity, status, table);
                }

                // Note: Explicit "Differs in: ..." summary removed as it's redundant with the visual diff table.
            } else {
                extra_count += 1;
                extra_len_sum += r.interaction.len();
                extra_energy_sum += r.energy.parse::<f64>().unwrap_or(0.0);
                let table = ParityTable {
                    kind: ParityKind::RustOnly(r),
                    config: TableConfig::default(),
                };
                // Print table as single block
                debug!("{}   {}\n{}", LogTag::Parity, HitStatus::Extra, table);
            }
        }

        for (i, c) in c_remaining.iter().enumerate() {
            if !c_rem_matched[i] {
                missing_count += 1;
                missing_len_sum += c.interaction.len();
                missing_energy_sum += c.energy.parse::<f64>().unwrap_or(0.0);
                let table = ParityTable {
                    kind: ParityKind::COnly(c),
                    config: TableConfig::default(),
                };
                // Print table as single block
                debug!("{}   {}\n{}", LogTag::Parity, HitStatus::Missing, table);
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
        info!("{} {}", LogTag::Summary, line);
    }

    // Handle allowed divergences
    if is_allowed_divergence {
        warn!(
            "{} Allowed divergence for '{}': Missing in Rust expected due to Maximality/Wobble improvements",
            LogTag::Parity,
            test_name
        );
        return;
    }

    if !is_pass {
        panic!(
            "\n{} FAILED: {} ({} rust-worse, {} missing, {} extra)\nRun with RUST_LOG=debug for detailed diff analysis.",
            LogTag::Parity,
            test_name,
            rust_worse_energy,
            missing_count,
            extra_count
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
pub fn run_single_seq_parity(
    query_seq: &str,
    target_seq: &str,
    test_name: &str,
    args: &[&str],
    detailed: bool,
) {
    init_test_logging();

    let root = workspace_root();
    let tmpdir = tempfile::tempdir().expect("tempdir");

    let query_path = tmpdir.path().join("query.fa");
    let target_path = tmpdir.path().join("target.fa");

    // Force uppercase for compatibility
    let q_upper = query_seq.to_uppercase();
    let t_upper = target_seq.to_uppercase();

    fs::write(&query_path, format!(">query\n{}\n", q_upper)).expect("write query");
    fs::write(&target_path, format!(">target\n{}\n", t_upper)).expect("write target");

    let c_bin = c_binary_path(&root);
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
