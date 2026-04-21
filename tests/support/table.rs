//! Table rendering for parity debug output.
//!
//! Provides `ParityTable` and related types for visualizing alignment
//! differences between Rust and C implementations.

use crate::support::search_hit::SearchHitExt;
use risearch::index::store::TargetStore;
use risearch::types::{Base, Strand};
use risearch::{QueryRegistry, SearchHit};
use tabled::{builder::Builder, settings::Style, Table, Tabled};

// =============================================================================
// SUMMARY TABLE
// =============================================================================

/// Row for summary tables.
#[derive(Tabled)]
pub(crate) struct SummaryRow {
    #[tabled(rename = "Metric")]
    pub metric: String,
    #[tabled(rename = "Count")]
    pub count: String,
}

impl SummaryRow {
    pub(crate) fn new(metric: &str, count: impl ToString) -> Self {
        Self {
            metric: metric.to_string(),
            count: count.to_string(),
        }
    }
}

/// Render a list of summary rows as a table string.
pub(crate) fn render_summary_table(rows: Vec<SummaryRow>) -> String {
    Table::new(rows).with(Style::rounded()).to_string()
}

// =============================================================================
// PARITY KIND
// =============================================================================

/// The type of parity comparison being displayed.
#[derive(Debug)]
#[allow(dead_code)] // Some variants used only in detailed debugging
pub(crate) enum ParityKind<'a> {
    /// Mismatch between Rust and C results (same coordinates, different content)
    Mismatch {
        rust: &'a SearchHit,
        c: &'a SearchHit,
    },
    /// Hit only in Rust output
    RustOnly(&'a SearchHit),
    /// Hit only in C output (no overlapping Rust hit)
    COnly(&'a SearchHit),
    /// C hit covered by overlapping Rust hit (different coordinates)
    CoveredBy {
        c: &'a SearchHit,
        rust: &'a SearchHit,
    },
}

// =============================================================================
// COLUMN KIND
// =============================================================================

/// Column types for parity tables.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum ColumnKind {
    Label,
    Ctx5,
    Ext5,
    Seed,
    Ext3,
    Ctx3,
}

impl ColumnKind {
    pub(crate) fn header(&self) -> &'static str {
        match self {
            Self::Label => "",
            Self::Ctx5 => "5' Ctx",
            Self::Ext5 => "5' EXT",
            Self::Seed => "SEED",
            Self::Ext3 => "3' EXT",
            Self::Ctx3 => "3' Ctx",
        }
    }
}

// =============================================================================
// TABLE CONFIG
// =============================================================================

/// Configuration for parity table output.
pub(crate) struct TableConfig {
    columns: Vec<ColumnKind>,
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

// =============================================================================
// ROW LABEL
// =============================================================================

/// Labels for table rows.
#[derive(Debug, Clone, Copy)]
pub(crate) enum RowLabel {
    SingleFP,
    SingleTarget,
    SingleQuery,
    CompCFP,
    CompCTarget,
    CompCQuery,
    CompDiff,
    CompRFP,
    CompRTarget,
    CompRQuery,
}

impl std::fmt::Display for RowLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SingleFP => write!(f, "FP"),
            Self::SingleTarget => write!(f, "Target"),
            Self::SingleQuery => write!(f, "Query"),
            Self::CompCFP => write!(f, "C FP"),
            Self::CompCTarget => write!(f, "C Tgt"),
            Self::CompCQuery => write!(f, "C Qry"),
            Self::CompDiff => write!(f, "DIFF"),
            Self::CompRFP => write!(f, "R FP"),
            Self::CompRTarget => write!(f, "R Tgt"),
            Self::CompRQuery => write!(f, "R Qry"),
        }
    }
}

// =============================================================================
// DIFF CHAR
// =============================================================================

/// Character-level diff indicator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DiffChar {
    Match,    // ' '
    Mismatch, // 'X'
}

impl DiffChar {
    pub(crate) fn as_char(&self) -> char {
        match self {
            Self::Match => ' ',
            Self::Mismatch => 'X',
        }
    }
}

impl From<(char, char)> for DiffChar {
    fn from((a, b): (char, char)) -> Self {
        if a == b {
            Self::Match
        } else {
            Self::Mismatch
        }
    }
}

/// Build a diff string comparing two strings character by character.
pub(crate) fn build_diff(a: &str, b: &str) -> String {
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

// =============================================================================
// PARSED INTERACTION
// =============================================================================

/// Parsed interaction string split into components.
struct ParsedInteraction {
    pub ctx_5: String,
    pub ext_5: String,
    pub seed: String,
    pub ext_3: String,
    pub ctx_3: String,
}

impl ParsedInteraction {
    /// Parse from a hit using its seed range.
    fn from_hit(hit: &SearchHit) -> Self {
        Self::from_hit_with_range(hit, hit.seed_start(), hit.seed_end())
    }

    /// Parse with explicit seed range override (for diff alignment).
    fn from_hit_with_range(
        hit: &SearchHit,
        seed_start: Option<usize>,
        seed_end: Option<usize>,
    ) -> Self {
        let fp_raw = hit
            .fingerprint()
            .expect("hit missing alignment for fingerprint");

        let (left, seed, right) = if let (Some(start), Some(end)) = (seed_start, seed_end) {
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
        } else {
            ("".into(), fp_raw.clone(), "".into())
        };

        Self {
            ctx_5: "".into(),
            ext_5: left,
            seed,
            ext_3: right,
            ctx_3: "".into(),
        }
    }

    /// Parse target sequence using reference parts for alignment.
    fn from_hit_target(
        hit: &SearchHit,
        ref_parts: &ParsedInteraction,
        target_store: Option<&TargetStore>,
    ) -> Self {
        let chars: Vec<char> = target_store
            .and_then(|ts| aligned_target_track(hit, ts))
            .unwrap_or_default()
            .chars()
            .collect();

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

    /// Parse query sequence using reference parts for alignment.
    fn from_hit_query(
        hit: &SearchHit,
        ref_parts: &ParsedInteraction,
        query_registry: Option<&QueryRegistry>,
    ) -> Self {
        let chars: Vec<char> = query_registry
            .and_then(|qr| aligned_query_track(hit, qr))
            .unwrap_or_default()
            .chars()
            .collect();

        let l_len = ref_parts.ext_5.chars().count();
        let s_len = ref_parts.seed.chars().count();

        let total = chars.len();
        let p1 = l_len.min(total);
        let p2 = (l_len + s_len).min(total);

        Self {
            ctx_5: "".into(), // Query doesn't have context flanks in same sense
            ext_5: chars[..p1].iter().collect(),
            seed: chars[p1..p2].iter().collect(),
            ext_3: chars[p2..].iter().collect(),
            ctx_3: "".into(),
        }
    }
}

fn hit_query_bases<'a>(hit: &SearchHit, query_registry: &'a QueryRegistry) -> &'a [Base] {
    let q_seq = query_registry.get(hit.query_idx).sequence();
    let seq_len = q_seq.len();
    let consumed = hit
        .alignment
        .as_ref()
        .map(|a| a.steps().iter().filter(|s| s.consumes_query()).count())
        .unwrap_or(0);
    if consumed == 0 || seq_len == 0 {
        return &q_seq[0..0];
    }
    let candidates = [hit.q_start, hit.q_start.saturating_sub(1)];
    let (start, _) = candidates
        .iter()
        .copied()
        .map(|c| c.min(seq_len))
        .map(|start| {
            let expected_end = start.saturating_add(consumed.saturating_sub(1));
            let score = expected_end.abs_diff(hit.q_end)
                + expected_end.saturating_add(1).abs_diff(hit.q_end);
            (start, score)
        })
        .min_by_key(|(_, score)| *score)
        .unwrap_or((0, usize::MAX));
    let end_excl = start.saturating_add(consumed).min(seq_len);
    &q_seq[start..end_excl]
}

fn hit_target_bases<'a>(hit: &SearchHit, target_store: &'a TargetStore) -> &'a [Base] {
    let t_idx = hit.target_idx as usize;
    let (_, t_fwd, t_rc, _) = match target_store.target_seqs(t_idx) {
        Ok(s) => s,
        Err(_) => return &[],
    };
    match hit.strand {
        Strand::Forward => {
            let start = hit.t_start.min(t_fwd.len());
            let end_excl = hit.t_end.saturating_add(1).min(t_fwd.len());
            if end_excl < start {
                &t_fwd[0..0]
            } else {
                &t_fwd[start..end_excl]
            }
        }
        Strand::Reverse => {
            let len = t_fwd.len();
            if len == 0 {
                return &t_rc[0..0];
            }
            let rc_start = len.saturating_sub(hit.t_end.saturating_add(1));
            let rc_end_incl = len.saturating_sub(hit.t_start.saturating_add(1));
            let start = rc_start.min(t_rc.len());
            let end_excl = rc_end_incl.saturating_add(1).min(t_rc.len());
            if end_excl < start {
                &t_rc[0..0]
            } else {
                &t_rc[start..end_excl]
            }
        }
    }
}

fn build_track<F>(hit: &SearchHit, bases: &[Base], consumes: F) -> Option<String>
where
    F: Fn(risearch::PairClass) -> bool,
{
    let alignment = hit.alignment.as_ref()?;
    let mut idx = 0usize;
    let mut out = String::with_capacity(alignment.steps().len());
    for &step in alignment.steps() {
        if consumes(step) {
            let b = bases.get(idx).copied().unwrap_or(Base::Gap);
            out.push(b.to_byte() as char);
            idx += 1;
        } else {
            out.push(Base::Gap.to_byte() as char);
        }
    }
    Some(out)
}

fn aligned_query_track(hit: &SearchHit, query_registry: &QueryRegistry) -> Option<String> {
    let q_bases = hit_query_bases(hit, query_registry);
    build_track(hit, q_bases, |s| s.consumes_query())
}

fn aligned_target_track(hit: &SearchHit, target_store: &TargetStore) -> Option<String> {
    let t_bases = hit_target_bases(hit, target_store);
    build_track(hit, t_bases, |s| s.consumes_target())
}

// =============================================================================
// PARITY TABLE
// =============================================================================

/// Table for displaying parity comparison results.
pub(crate) struct ParityTable<'a> {
    pub(crate) kind: ParityKind<'a>,
    pub(crate) config: TableConfig,
    pub(crate) query_registry: Option<&'a QueryRegistry>,
    pub(crate) target_store: Option<&'a TargetStore>,
}

impl<'a> std::fmt::Display for ParityTable<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut builder = Builder::default();

        let strand = match &self.kind {
            ParityKind::Mismatch { c, .. } => &c.strand,
            ParityKind::RustOnly(r) => &r.strand,
            ParityKind::COnly(c) => &c.strand,
            ParityKind::CoveredBy { c, .. } => &c.strand,
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
                let p = ParsedInteraction::from_hit(r);
                let p_tgt = ParsedInteraction::from_hit_target(r, &p, self.target_store);
                let p_qry = ParsedInteraction::from_hit_query(r, &p, self.query_registry);
                add_row(&mut builder, RowLabel::SingleTarget, &p_tgt);
                add_row(&mut builder, RowLabel::SingleQuery, &p_qry);
                add_row(&mut builder, RowLabel::SingleFP, &p);
            }
            ParityKind::COnly(c) => {
                let p = ParsedInteraction::from_hit(c);
                let p_tgt = ParsedInteraction::from_hit_target(c, &p, self.target_store);
                let p_qry = ParsedInteraction::from_hit_query(c, &p, self.query_registry);
                add_row(&mut builder, RowLabel::SingleTarget, &p_tgt);
                add_row(&mut builder, RowLabel::SingleQuery, &p_qry);
                add_row(&mut builder, RowLabel::SingleFP, &p);
            }
            ParityKind::Mismatch { rust: r, c } => {
                let p_c = ParsedInteraction::from_hit(c);
                let p_r = ParsedInteraction::from_hit(r); // Use Rust's own seed range

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

                let p_c_tgt = ParsedInteraction::from_hit_target(c, &p_c, self.target_store);
                // Use Rust query for C since C output lacks query seq but they're same hit
                let p_c_qry = ParsedInteraction::from_hit_query(r, &p_c, self.query_registry);
                let p_r_tgt = ParsedInteraction::from_hit_target(r, &p_r, self.target_store);
                let p_r_qry = ParsedInteraction::from_hit_query(r, &p_r, self.query_registry);

                add_row(&mut builder, RowLabel::CompCTarget, &p_c_tgt);
                add_row(&mut builder, RowLabel::CompCQuery, &p_c_qry);
                add_row(&mut builder, RowLabel::CompCFP, &p_c);
                add_row(&mut builder, RowLabel::CompDiff, &p_diff);
                add_row(&mut builder, RowLabel::CompRFP, &p_r);
                add_row(&mut builder, RowLabel::CompRQuery, &p_r_qry);
                add_row(&mut builder, RowLabel::CompRTarget, &p_r_tgt);
            }
            ParityKind::CoveredBy { c, rust: r } => {
                // Show C hit (missing), then separator, then overlapping Rust hit
                let p_c = ParsedInteraction::from_hit(c);
                let p_c_tgt = ParsedInteraction::from_hit_target(c, &p_c, self.target_store);

                let p_r = ParsedInteraction::from_hit(r);
                let p_r_tgt = ParsedInteraction::from_hit_target(r, &p_r, self.target_store);
                let p_r_qry = ParsedInteraction::from_hit_query(r, &p_r, self.query_registry);

                // C hit rows
                add_row(&mut builder, RowLabel::CompCTarget, &p_c_tgt);
                add_row(&mut builder, RowLabel::CompCFP, &p_c);
                // Separator row (simple divider)
                let sep = ParsedInteraction {
                    ctx_5: "".into(),
                    ext_5: "───────".into(),
                    seed: "OVERLAP".into(),
                    ext_3: "───────".into(),
                    ctx_3: "".into(),
                };
                add_row(&mut builder, RowLabel::CompDiff, &sep);
                // Rust hit rows
                add_row(&mut builder, RowLabel::CompRFP, &p_r);
                add_row(&mut builder, RowLabel::CompRQuery, &p_r_qry);
                add_row(&mut builder, RowLabel::CompRTarget, &p_r_tgt);
            }
        }

        write!(f, "{}", builder.build().with(Style::rounded()))
    }
}
