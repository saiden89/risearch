//! Record types for parity testing.
//!
//! Contains `Rec` - the unified record type for comparing Rust and C output.

use log::trace;

// =============================================================================
// SEED MARKERS
// =============================================================================

/// Markers used in C debug output to delimit seed regions.
/// Format: `<left_ext>y<seed>x<right_ext>`
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

// =============================================================================
// REC STRUCT
// =============================================================================

/// A unified record representing a hit from either Rust or C output.
///
/// Used for comparing parity between implementations. Fields are kept as
/// strings where exact comparison matters (energy), or parsed where fuzzy
/// comparison is needed.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Rec {
    // Identity fields for grouping
    pub q_id: String,
    pub t_id: String,

    // Coordinate fields
    pub q_start: usize,
    pub q_end: usize,
    pub t_start: usize,
    pub t_end: usize,
    pub strand: String,

    // Data fields
    pub energy: String, // Keep as string for exact comparison, parse for fuzzy
    pub interaction: String,
    pub target_seq: String,
    pub query_seq: String,
    pub flank_5: String,
    pub flank_3: String,

    // Seed range within interaction (parsed from markers if present)
    pub seed_start: Option<usize>,
    pub seed_end: Option<usize>,
}

impl Rec {
    /// Parse a record from a tab-separated line (C or Rust output format).
    ///
    /// Expects at least 10 columns:
    /// `q_id, q_start, q_end, t_id, t_start, t_end, strand, energy, interaction, target_seq`
    ///
    /// Optional columns 11-13: `flank_5, flank_3, query_seq`
    pub fn from_line(line: &str) -> Option<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
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

        // Parse interaction with seed range extraction
        let (interaction, seed_start, seed_end) = Self::parse_and_strip_markers(fields[8]);

        Some(Rec {
            q_id: fields[0].to_string(),
            q_start: fields[1].parse().unwrap_or(0),
            q_end: fields[2].parse().unwrap_or(0),
            t_id: fields[3].to_string(),
            t_start: fields[4].parse().unwrap_or(0),
            t_end: fields[5].parse().unwrap_or(0),
            strand: fields[6].to_string(),
            energy: fields[7].to_string(),
            interaction,
            target_seq: Self::strip_markers(fields[9]),
            // Optional fields - use get() for safe access
            flank_5: fields
                .get(10)
                .map_or(String::new(), |s| Self::strip_markers(s)),
            flank_3: fields
                .get(11)
                .map_or(String::new(), |s| Self::strip_markers(s)),
            query_seq: fields.get(12).map_or(String::new(), |s| s.to_string()),
            seed_start,
            seed_end,
        })
    }

    /// Strip y/x debug markers and extract seed range.
    fn parse_and_strip_markers(s: &str) -> (String, Option<usize>, Option<usize>) {
        let before = char::from(SeedMarker::BeforeSeed);
        let after = char::from(SeedMarker::AfterSeed);

        let y_pos = s.find(before);
        let x_pos = s.find(after);

        let (seed_start, seed_end) = match (y_pos, x_pos) {
            (Some(y), Some(x)) if x > y => {
                // y marks start, x marks end (after y was removed = x-1)
                (Some(y), Some(x - 1))
            }
            _ => (None, None),
        };

        let clean: String = s.chars().filter(|&c| c != before && c != after).collect();
        (clean, seed_start, seed_end)
    }

    /// Strip y/x markers without tracking positions.
    fn strip_markers(s: &str) -> String {
        let before = char::from(SeedMarker::BeforeSeed);
        let after = char::from(SeedMarker::AfterSeed);
        s.chars().filter(|&c| c != before && c != after).collect()
    }

    /// Returns true if coordinates and strand match.
    pub fn coords_match(&self, other: &Self) -> bool {
        self.q_start == other.q_start
            && self.q_end == other.q_end
            && self.t_start == other.t_start
            && self.t_end == other.t_end
            && self.strand == other.strand
    }

    /// Format record for debug output.
    #[allow(dead_code)] // Useful for debugging
    pub fn fmt_coords(&self) -> String {
        format!(
            "q=[{},{}] t=[{},{}] S={} E={}",
            self.q_start, self.q_end, self.t_start, self.t_end, self.strand, self.energy
        )
    }
}

// =============================================================================
// FROM SEARCHHIT
// =============================================================================

impl From<&risearch::SearchHit> for Rec {
    fn from(hit: &risearch::SearchHit) -> Self {
        // Truncate IDs at whitespace to match C output format
        let q_id = hit
            .query_id
            .split_whitespace()
            .next()
            .unwrap_or(&hit.query_id)
            .to_string();
        let t_id = hit
            .target_id
            .split_whitespace()
            .next()
            .unwrap_or(&hit.target_id)
            .to_string();

        Rec {
            q_id,
            q_start: hit.q_start + 1, // Convert 0-based to 1-based
            q_end: hit.q_end + 1,
            t_id,
            t_start: hit.output_t_start,
            t_end: hit.output_t_end,
            strand: hit.strand.to_string(),
            energy: format!("{:.2}", hit.energy),
            interaction: hit
                .alignment
                .fingerprint()
                .replace('T', "U")
                .replace('t', "u"),
            target_seq: hit
                .alignment
                .target_sequence()
                .replace('T', "U")
                .replace('t', "u"),
            flank_5: hit.flank_5.clone(),
            flank_3: hit.flank_3.clone(),
            query_seq: String::new(),
            seed_start: Some(hit.alignment.left_extension().len()),
            seed_end: Some(hit.alignment.left_extension().len() + hit.alignment.seed().len()),
        }
    }
}

// =============================================================================
// UNIT TESTS
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_rec(q_start: usize, q_end: usize, t_start: usize, t_end: usize, strand: &str) -> Rec {
        Rec {
            q_id: "query".into(),
            t_id: "target".into(),
            q_start,
            q_end,
            t_start,
            t_end,
            strand: strand.into(),
            energy: "-10.0".into(),
            interaction: "PPPPP".into(),
            target_seq: "AUGCG".into(),
            query_seq: String::new(),
            flank_5: String::new(),
            flank_3: String::new(),
            seed_start: None,
            seed_end: None,
        }
    }

    #[test]
    fn test_rec_from_line_minimal() {
        let line = "query\t1\t10\ttarget\t100\t110\t+\t-12.5\tPPPPP\tAUGCG";
        let rec = Rec::from_line(line).unwrap();
        assert_eq!(rec.q_start, 1);
        assert_eq!(rec.energy, "-12.5");
        assert_eq!(rec.interaction, "PPPPP");
    }

    #[test]
    fn test_rec_from_line_with_seed_markers() {
        let line = "q\t1\t5\tt\t1\t5\t+\t-5.0\tyPPPx\tAUGCG";
        let rec = Rec::from_line(line).unwrap();
        assert_eq!(rec.interaction, "PPP"); // markers stripped
        assert_eq!(rec.seed_start, Some(0));
        assert_eq!(rec.seed_end, Some(3));
    }

    #[test]
    fn test_rec_from_line_too_few_columns() {
        let line = "query\t1\t10\ttarget";
        assert!(Rec::from_line(line).is_none());
    }

    #[test]
    fn test_rec_coords_match() {
        let r1 = make_rec(1, 10, 100, 110, "+");
        let r2 = make_rec(1, 10, 100, 110, "+");
        let r3 = make_rec(1, 10, 100, 111, "+");
        assert!(r1.coords_match(&r2));
        assert!(!r1.coords_match(&r3));
    }

    #[test]
    fn test_strip_markers() {
        assert_eq!(Rec::strip_markers("yPPPx"), "PPP");
        assert_eq!(Rec::strip_markers("PPPP"), "PPPP");
    }
}
