mod support;

use risearch::types::{Energy, Strand};
use risearch::{Alignment, PairClass, SearchHit};
use rstest::rstest;
use std::path::PathBuf;
use support::comparison::{
    classify_missing, hits_overlap, ranges_overlap, ParityComparator, ParityResult,
};
use support::status::{HitStatus, MissingReason, ParityMode};
use support::table::{build_diff, ColumnKind, DiffChar, RowLabel};

fn make_hit(
    q_start: usize,
    q_end: usize,
    t_start: usize,
    t_end: usize,
    strand: Strand,
    energy: f64,
) -> SearchHit {
    let len = q_end.saturating_sub(q_start).max(1);
    let seed: Vec<PairClass> = (0..len).map(|_| PairClass::Canonical).collect();
    let alignment = Alignment::new(&[], &seed, &[]);

    SearchHit {
        query_idx: 0,
        target_idx: 0,
        q_start,
        q_end,
        t_start,
        t_end,
        strand,
        energy: Energy::from_kcal(energy),
        seed_start: None,
        seed_end: None,
        alignment: Some(alignment),
    }
}

#[test]
fn c_binary_path_structure() {
    let root = PathBuf::from("/fake/root");
    let debug_path = root.join("legacy_c/RIsearch2/bin/risearch2.dbg.x");
    let release_path = root.join("legacy_c/RIsearch2/bin/risearch2.x");

    assert!(debug_path.to_str().unwrap().contains("dbg"));
    assert!(!release_path.to_str().unwrap().contains("dbg"));
}

#[test]
fn ranges_overlap_cases() {
    assert!(ranges_overlap(1, 10, 5, 15));
    assert!(ranges_overlap(1, 10, 1, 10));
    assert!(ranges_overlap(1, 10, 3, 7));
    assert!(!ranges_overlap(1, 10, 11, 20));
    assert!(ranges_overlap(1, 10, 10, 20));
}

#[test]
fn hits_overlap_checks_strand() {
    let a = make_hit(1, 10, 100, 110, Strand::Forward, -10.0);
    let b = make_hit(5, 15, 105, 115, Strand::Forward, -12.0);
    let c = make_hit(1, 10, 100, 110, Strand::Reverse, -10.0);
    assert!(hits_overlap(&a, &b));
    assert!(!hits_overlap(&a, &c));
}

#[test]
fn classify_missing_better_energy() {
    let c_hit = make_hit(1, 10, 100, 110, Strand::Forward, -10.0);
    let rust_hits = [make_hit(1, 10, 100, 110, Strand::Forward, -15.0)];
    let refs: Vec<&SearchHit> = rust_hits.iter().collect();
    let (reason, _) = classify_missing(&c_hit, &refs);
    assert_eq!(reason, MissingReason::BetterEnergy);
}

#[test]
fn comparator_exact_match() {
    let rust = vec![make_hit(1, 10, 100, 110, Strand::Forward, -10.0)];
    let c = vec![make_hit(1, 10, 100, 110, Strand::Forward, -10.0)];
    let result = ParityComparator::new(&rust, &c).compare();
    assert_eq!(result.exact_matches, 1);
    assert!(result.is_pass(ParityMode::Relaxed));
}

#[test]
fn comparator_rust_better() {
    let rust = vec![make_hit(1, 10, 100, 110, Strand::Forward, -15.0)];
    let c = vec![make_hit(1, 10, 100, 110, Strand::Forward, -10.0)];
    let result = ParityComparator::new(&rust, &c).compare();
    assert_eq!(result.rust_better.len(), 1);
    assert!(result.is_pass(ParityMode::Relaxed));
}

#[test]
fn comparator_rust_worse_fails() {
    let rust = vec![make_hit(1, 10, 100, 110, Strand::Forward, -5.0)];
    let c = vec![make_hit(1, 10, 100, 110, Strand::Forward, -10.0)];
    let result = ParityComparator::new(&rust, &c).compare();
    assert_eq!(result.rust_worse.len(), 1);
    assert!(!result.is_pass(ParityMode::Relaxed));
}

#[test]
fn log_details_covers_all_categories() {
    let rust = vec![
        make_hit(1, 10, 100, 110, Strand::Forward, -10.0), // exact match
        make_hit(1, 10, 200, 210, Strand::Forward, -15.0), // rust_better
        make_hit(1, 10, 300, 310, Strand::Forward, -5.0),  // rust_worse
        make_hit(1, 10, 400, 410, Strand::Forward, -8.0),  // extra
    ];
    let c = vec![
        make_hit(1, 10, 100, 110, Strand::Forward, -10.0), // exact match
        make_hit(1, 10, 200, 210, Strand::Forward, -10.0), // covered by rust_better
        make_hit(1, 10, 300, 310, Strand::Forward, -10.0), // covered by rust_worse
        make_hit(1, 10, 500, 510, Strand::Forward, -8.0),  // missing (no overlap)
    ];

    let result = ParityComparator::new(&rust, &c).compare();
    assert_eq!(result.exact_matches, 1);
    assert_eq!(result.rust_better.len(), 1);
    assert_eq!(result.rust_worse.len(), 1);
    assert_eq!(result.extras.len(), 1);
    assert_eq!(result.missings.len(), 1);
    result.log_details("test_all_categories");
}

#[test]
fn log_details_empty_result() {
    let result = ParityResult::default();
    result.log_details("empty_test");
}

#[test]
fn log_details_co_optimal() {
    let rust_hit = make_hit(1, 10, 100, 110, Strand::Forward, -10.0);
    let c_hit = make_hit(1, 10, 100, 110, Strand::Forward, -10.0);
    let result = ParityComparator::new(&[rust_hit], &[c_hit]).compare();
    assert!(result.exact_matches >= 1 || !result.co_optimal.is_empty());
    result.log_details("co_optimal_test");
}

#[rstest]
#[case(ParityMode::Absolute, false, false, false, false)]
#[case(ParityMode::Strict, false, false, false, false)]
#[case(ParityMode::Balanced, false, false, false, false)]
#[case(ParityMode::Relaxed, true, true, false, false)]
fn missing_reason_acceptable(
    #[case] mode: ParityMode,
    #[case] better_ok: bool,
    #[case] equal_ok: bool,
    #[case] worse_ok: bool,
    #[case] no_overlap_ok: bool,
) {
    assert_eq!(MissingReason::BetterEnergy.is_acceptable(mode), better_ok);
    assert_eq!(MissingReason::EqualEnergy.is_acceptable(mode), equal_ok);
    assert_eq!(MissingReason::WorseEnergy.is_acceptable(mode), worse_ok);
    assert_eq!(MissingReason::NoOverlap.is_acceptable(mode), no_overlap_ok);
}

#[test]
fn hit_status_display() {
    assert_eq!(format!("{}", HitStatus::RustBetter), "RUST BETTER");
    assert_eq!(format!("{}", HitStatus::Missing), "MISSING");
    assert_eq!(format!("{}", HitStatus::CoOptimal), "CO-OPTIMAL");
}

#[test]
fn parity_mode_default_is_balanced() {
    assert_eq!(ParityMode::default(), ParityMode::Balanced);
}

#[test]
fn build_diff_identical() {
    assert_eq!(build_diff("PPPP", "PPPP"), "    ");
}

#[test]
fn build_diff_mismatch() {
    assert_eq!(build_diff("PPUP", "PPPP"), "  X ");
}

#[test]
fn build_diff_length_mismatch() {
    let diff = build_diff("PPP", "PPPPP");
    assert_eq!(diff.len(), 5);
}

#[test]
fn diff_char_conversion() {
    assert_eq!(DiffChar::from(('P', 'P')).as_char(), ' ');
    assert_eq!(DiffChar::from(('P', 'U')).as_char(), 'X');
}

#[test]
fn column_kind_headers() {
    assert_eq!(ColumnKind::Seed.header(), "SEED");
    assert_eq!(ColumnKind::Ext5.header(), "5' EXT");
}

#[test]
fn row_label_display() {
    assert_eq!(format!("{}", RowLabel::CompDiff), "DIFF");
    assert_eq!(format!("{}", RowLabel::SingleFP), "FP");
}
