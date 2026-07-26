//! Parity for the detailed (`-p1`) alignment block: query row, pairing symbols,
//! target row.
//!
//! Every other parity suite goes through `RustRunner::search`, which forces
//! binding-site columns, so the three rendered rows have no coverage anywhere
//! else. Case is normalized before comparing: C echoes the query in the input
//! FASTA's case and the target in lowercase, Rust lowercases both.

mod support;

use std::collections::BTreeMap;

use support::SingleSeqRunner;

/// Fully complementary duplex except for one mismatch at query position 7.
/// Off-center on purpose: a reversed flank cannot coincide with the correct one.
const QUERY: &str = "UGCUGACGUACGU";
const TARGET: &str = "ACGUACAUCAGCA";

const BASE_ARGS: &[&str] = &["-l", "20", "-e", "100", "-p1"];

#[test]
fn seed_at_query_start_extends_right() {
    assert_blocks_match(
        "extend_right",
        &["--seed-start", "1", "--seed-end", "5", "--seed-length", "5"],
    );
}

#[test]
fn seed_at_query_end_extends_left() {
    assert_blocks_match(
        "extend_left",
        &[
            "--seed-start",
            "9",
            "--seed-end",
            "13",
            "--seed-length",
            "5",
        ],
    );
}

#[test]
fn unconstrained_seed() {
    assert_blocks_match("unconstrained", &["--seed-length", "5"]);
}

/// Alignment rows keyed by hit coordinates. Energy is deliberately excluded:
/// the TSV suites own energy parity, this one owns the rendered rows.
type Blocks = BTreeMap<String, Vec<[String; 3]>>;

fn assert_blocks_match(name: &str, seed_args: &[&str]) {
    let args: Vec<&str> = seed_args.iter().chain(BASE_ARGS).copied().collect();
    let (rust_out, c_out) = SingleSeqRunner::new(QUERY, TARGET).rendered_texts(&args);

    let rust = parse_blocks(&rust_out);
    let c = parse_blocks(&c_out);

    assert!(!c.is_empty(), "{name}: C reported no hits to compare");
    assert_eq!(
        rust.keys().collect::<Vec<_>>(),
        c.keys().collect::<Vec<_>>(),
        "{name}: hit coordinates differ"
    );

    let mut report = String::new();
    for (key, c_blocks) in &c {
        let rust_blocks = &rust[key];
        if rust_blocks != c_blocks {
            report.push_str(&format!("\n  {key}\n"));
            for (label, blocks) in [("rust", rust_blocks), ("c", c_blocks)] {
                for rows in blocks {
                    for (n, row) in rows.iter().enumerate() {
                        let tag = if n == 0 { label } else { "" };
                        report.push_str(&format!("    {tag:>4} {row:?}\n"));
                    }
                }
            }
        }
    }
    assert!(
        report.is_empty(),
        "{name}: detailed alignment blocks differ:{report}"
    );
}

/// Detailed output is one 4-line record per hit: the three alignment rows, then
/// the TSV line. Only the TSV line contains tabs.
fn parse_blocks(out: &str) -> Blocks {
    let mut blocks = Blocks::new();
    let mut pending: Vec<String> = Vec::new();

    for line in out.lines() {
        if !line.contains('\t') {
            if !line.is_empty() {
                pending.push(normalize(line));
            }
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(fields.len() >= 7, "short detailed TSV line: {line}");
        assert_eq!(
            pending.len(),
            3,
            "expected 3 alignment rows before {line}, got {pending:?}"
        );
        let rows: [String; 3] = std::mem::take(&mut pending).try_into().unwrap();
        blocks.entry(fields[..7].join(" ")).or_default().push(rows);
    }
    assert!(
        pending.is_empty(),
        "alignment rows with no TSV line: {pending:?}"
    );

    for hits in blocks.values_mut() {
        hits.sort();
    }
    blocks
}

/// Drop the debug C binary's `y`/`x` seed markers and unify case. Whitespace is
/// left alone: a trailing space in the symbol row marks an unpaired last column.
fn normalize(row: &str) -> String {
    row.chars()
        .filter(|&c| c != 'y' && c != 'x')
        .map(|c| c.to_ascii_lowercase())
        .collect()
}
