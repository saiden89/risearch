//! Table rendering for parity debug output.
//!
//! Provides `ParityTable` and related types for visualizing alignment
//! differences between Rust and C implementations.

use risearch::SearchHit;
use tabled::{Table, Tabled, builder::Builder, settings::Style};

// =============================================================================
// SUMMARY TABLE
// =============================================================================

/// Row for summary tables.
#[derive(Tabled)]
pub struct SummaryRow {
    #[tabled(rename = "Metric")]
    pub metric: String,
    #[tabled(rename = "Count")]
    pub count: String,
}

impl SummaryRow {
    pub fn new(metric: &str, count: impl ToString) -> Self {
        Self {
            metric: metric.to_string(),
            count: count.to_string(),
        }
    }
}

/// Render a list of summary rows as a table string.
pub fn render_summary_table(rows: Vec<SummaryRow>) -> String {
    Table::new(rows).with(Style::rounded()).to_string()
}

// =============================================================================
// PARITY KIND
// =============================================================================

/// The type of parity comparison being displayed.
#[derive(Debug)]
#[allow(dead_code)] // Some variants used only in detailed debugging
pub enum ParityKind<'a> {
    /// Mismatch between Rust and C results
    Mismatch {
        rust: &'a SearchHit,
        c: &'a SearchHit,
    },
    /// Hit only in Rust output
    RustOnly(&'a SearchHit),
    /// Hit only in C output
    COnly(&'a SearchHit),
}

// =============================================================================
// COLUMN KIND
// =============================================================================

/// Column types for parity tables.
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

// =============================================================================
// ROW LABEL
// =============================================================================

/// Labels for table rows.
#[derive(Debug, Clone, Copy)]
pub enum RowLabel {
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

/// Build a diff string comparing two strings character by character.
pub fn build_diff(a: &str, b: &str) -> String {
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
pub struct ParsedInteraction {
    pub ctx_5: String,
    pub ext_5: String,
    pub seed: String,
    pub ext_3: String,
    pub ctx_3: String,
}

impl ParsedInteraction {
    /// Parse from a hit using its seed range.
    pub fn from_hit(hit: &SearchHit) -> Self {
        Self::from_hit_with_range(hit, hit.seed_start(), hit.seed_end())
    }

    /// Parse with explicit seed range override (for diff alignment).
    pub fn from_hit_with_range(
        hit: &SearchHit,
        seed_start: Option<usize>,
        seed_end: Option<usize>,
    ) -> Self {
        let fp_raw = hit.fingerprint();

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
            ctx_5: hit.flank_5.clone(),
            ext_5: left,
            seed,
            ext_3: right,
            ctx_3: hit.flank_3.clone(),
        }
    }

    /// Parse target sequence using reference parts for alignment.
    pub fn from_hit_target(hit: &SearchHit, ref_parts: &ParsedInteraction) -> Self {
        let chars: Vec<char> = hit.target_seq().chars().collect();

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
    pub fn from_hit_query(hit: &SearchHit, ref_parts: &ParsedInteraction) -> Self {
        let chars: Vec<char> = hit.query_seq().chars().collect();

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

// =============================================================================
// PARITY TABLE
// =============================================================================

/// Table for displaying parity comparison results.
pub struct ParityTable<'a> {
    pub kind: ParityKind<'a>,
    pub config: TableConfig,
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
                let p = ParsedInteraction::from_hit(r);
                let p_tgt = ParsedInteraction::from_hit_target(r, &p);
                let p_qry = ParsedInteraction::from_hit_query(r, &p);
                add_row(&mut builder, RowLabel::SingleTarget, &p_tgt);
                add_row(&mut builder, RowLabel::SingleQuery, &p_qry);
                add_row(&mut builder, RowLabel::SingleFP, &p);
            }
            ParityKind::COnly(c) => {
                let p = ParsedInteraction::from_hit(c);
                let p_tgt = ParsedInteraction::from_hit_target(c, &p);
                let p_qry = ParsedInteraction::from_hit_query(c, &p);
                add_row(&mut builder, RowLabel::SingleTarget, &p_tgt);
                add_row(&mut builder, RowLabel::SingleQuery, &p_qry);
                add_row(&mut builder, RowLabel::SingleFP, &p);
            }
            ParityKind::Mismatch { rust: r, c } => {
                let p_c = ParsedInteraction::from_hit(c);
                let p_r = ParsedInteraction::from_hit_with_range(r, c.seed_start(), c.seed_end());

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

                let p_c_tgt = ParsedInteraction::from_hit_target(c, &p_c);
                // Use Rust query for C since C output lacks query seq but they're same hit
                let p_c_qry = ParsedInteraction::from_hit_query(r, &p_c);
                let p_r_tgt = ParsedInteraction::from_hit_target(r, &p_r);
                let p_r_qry = ParsedInteraction::from_hit_query(r, &p_r);

                add_row(&mut builder, RowLabel::CompCTarget, &p_c_tgt);
                add_row(&mut builder, RowLabel::CompCQuery, &p_c_qry);
                add_row(&mut builder, RowLabel::CompCFP, &p_c);
                add_row(&mut builder, RowLabel::CompDiff, &p_diff);
                add_row(&mut builder, RowLabel::CompRFP, &p_r);
                add_row(&mut builder, RowLabel::CompRQuery, &p_r_qry);
                add_row(&mut builder, RowLabel::CompRTarget, &p_r_tgt);
            }
        }

        write!(f, "{}", builder.build().with(Style::rounded()))
    }
}

// =============================================================================
// UNIT TESTS
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_diff_identical() {
        assert_eq!(build_diff("PPPP", "PPPP"), "    ");
    }

    #[test]
    fn test_build_diff_mismatch() {
        assert_eq!(build_diff("PPUP", "PPPP"), "  X ");
    }

    #[test]
    fn test_build_diff_length_mismatch() {
        let diff = build_diff("PPP", "PPPPP");
        assert_eq!(diff.len(), 5);
    }

    #[test]
    fn test_diff_char_conversion() {
        assert_eq!(DiffChar::from(('P', 'P')).as_char(), ' ');
        assert_eq!(DiffChar::from(('P', 'U')).as_char(), 'X');
    }

    #[test]
    fn test_column_kind_headers() {
        assert_eq!(ColumnKind::Seed.header(), "SEED");
        assert_eq!(ColumnKind::Ext5.header(), "5' EXT");
    }

    #[test]
    fn test_row_label_display() {
        assert_eq!(format!("{}", RowLabel::CompDiff), "DIFF");
        assert_eq!(format!("{}", RowLabel::SingleFP), "FP");
    }
}
