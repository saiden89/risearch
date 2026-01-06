//! Table rendering for parity debug output.
//!
//! Provides summary table utilities for parity test results.

use tabled::{Table, Tabled, settings::Style};

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
