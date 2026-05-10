#![allow(dead_code)]
//! RNA dinucleotide stacking energy matrices (DSM)
//!
//! Tables encode nearest-neighbor thermodynamic parameters for RNA-RNA interactions.
//! Runtime energy values are in 0.0001 kcal/mol (divide by -10000 for kcal/mol).
//!
//! Indices: `[q1][q2][t1][t2]` = stacking energy for:
//! ```text
//! Query:  5'─ q1 ─ q2 ─ 3'
//!             |    |
//! Target: 3'─ t1 ─ t2 ─ 5'
//! ```

use anyhow::{bail, Context, Result};
use std::collections::HashSet;
use std::path::Path;

use crate::types::{Base, DsmId, Energy, SequenceType, BASE_COUNT};

type DsmTable = [[[[i32; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];

const DSM_FLAT_SIZE: usize = BASE_COUNT * BASE_COUNT * BASE_COUNT * BASE_COUNT;

const MAX_TSV_ENERGY: f64 = 20.0;

pub const DSM_IDS: [&str; DsmId::ALL.len()] = {
    let mut ids = [""; DsmId::ALL.len()];
    let mut i = 0;
    while i < DsmId::ALL.len() {
        ids[i] = DsmId::ALL[i].as_str();
        i += 1;
    }
    ids
};
const CANONICAL_DSM_HEADER: [&str; 5] = ["q1", "q2", "t1", "t2", "delta_g_kcal_per_mol"];

#[derive(Clone, Copy)]
enum Orientation {
    Identity,
    ReverseSwap,
}

struct CanonicalDsms {
    id: DsmId,
    temperature: i32,
    initiation: f64,
    orientation: Orientation,
    tsv: &'static str,
}

static CANONICAL_TABLES: &[CanonicalDsms] = &[
    // t04: Turner 2004 RNA-RNA
    CanonicalDsms { id: DsmId::T04, temperature: 0,  initiation: 6.8982, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/t04/0.tsv")  },
    CanonicalDsms { id: DsmId::T04, temperature: 25, initiation: 6.3786, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/t04/25.tsv") },
    CanonicalDsms { id: DsmId::T04, temperature: 37, initiation: 6.128,  orientation: Orientation::Identity, tsv: include_str!("../data/dsm/t04/37.tsv") },
    CanonicalDsms { id: DsmId::T04, temperature: 42, initiation: 5.9848, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/t04/42.tsv") },
    CanonicalDsms { id: DsmId::T04, temperature: 50, initiation: 5.7678, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/t04/50.tsv") },
    // slh04: SantaLucia-Hicks 2004 DNA-DNA
    CanonicalDsms { id: DsmId::Slh04, temperature: 0,  initiation: 1.858,  orientation: Orientation::Identity, tsv: include_str!("../data/dsm/slh04/0.tsv")  },
    CanonicalDsms { id: DsmId::Slh04, temperature: 25, initiation: 1.9532, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/slh04/25.tsv") },
    CanonicalDsms { id: DsmId::Slh04, temperature: 37, initiation: 2.0438, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/slh04/37.tsv") },
    CanonicalDsms { id: DsmId::Slh04, temperature: 42, initiation: 2.1748, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/slh04/42.tsv") },
    CanonicalDsms { id: DsmId::Slh04, temperature: 50, initiation: 2.4016, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/slh04/50.tsv") },
    // s95-rna-dna: Sugimoto 1995 RNA-DNA (identity)
    CanonicalDsms { id: DsmId::S95RnaDna, temperature: 0,  initiation: 4.3476, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/s95/0.tsv")  },
    CanonicalDsms { id: DsmId::S95RnaDna, temperature: 25, initiation: 4.1588, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/s95/25.tsv") },
    CanonicalDsms { id: DsmId::S95RnaDna, temperature: 37, initiation: 4.035,  orientation: Orientation::Identity, tsv: include_str!("../data/dsm/s95/37.tsv") },
    CanonicalDsms { id: DsmId::S95RnaDna, temperature: 42, initiation: 3.9402, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/s95/42.tsv") },
    CanonicalDsms { id: DsmId::S95RnaDna, temperature: 50, initiation: 3.8672, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/s95/50.tsv") },
    // s95-dna-rna: Sugimoto 1995 DNA-RNA (reverse-swap of RNA-DNA)
    CanonicalDsms { id: DsmId::S95DnaRna, temperature: 0,  initiation: 4.3476, orientation: Orientation::ReverseSwap, tsv: include_str!("../data/dsm/s95/0.tsv")  },
    CanonicalDsms { id: DsmId::S95DnaRna, temperature: 25, initiation: 4.1588, orientation: Orientation::ReverseSwap, tsv: include_str!("../data/dsm/s95/25.tsv") },
    CanonicalDsms { id: DsmId::S95DnaRna, temperature: 37, initiation: 4.035,  orientation: Orientation::ReverseSwap, tsv: include_str!("../data/dsm/s95/37.tsv") },
    CanonicalDsms { id: DsmId::S95DnaRna, temperature: 42, initiation: 3.9402, orientation: Orientation::ReverseSwap, tsv: include_str!("../data/dsm/s95/42.tsv") },
    CanonicalDsms { id: DsmId::S95DnaRna, temperature: 50, initiation: 3.8672, orientation: Orientation::ReverseSwap, tsv: include_str!("../data/dsm/s95/50.tsv") },
    // t99: Turner 1999 RNA-RNA (37°C only)
    CanonicalDsms { id: DsmId::T99, temperature: 37, initiation: 5.59, orientation: Orientation::Identity, tsv: include_str!("../data/dsm/t99/37.tsv") },
];

/// Gap index used for DSM transition queries (linked to Base::Gap).
pub use crate::types::GAP;

// =============================================================================
// SCORING MODEL - Direction-agnostic penalty-adjusted transition table
// =============================================================================

/// Penalty-adjusted flat DSM table for one orientation.
///
/// Each instance is one canonical orientation (right or left).
/// Use `transpose()` to create the other orientation.
#[derive(Clone, Debug)]
pub struct ScoringModel {
    table: [i32; DSM_FLAT_SIZE],
    initiation: i32,
}

impl ScoringModel {
    /// Load from bundled canonical DSM tables, interpolating between bracket temperatures if needed.
    pub fn from_canonical(id: DsmId, temperature: i32, penalty: Energy) -> Result<Self> {
        let penalty = penalty.to_units();
        let entries: Vec<&CanonicalDsms> = CANONICAL_TABLES
            .iter()
            .filter(|e| e.id == id)
            .collect();

        if entries.is_empty() {
            bail!("No canonical DSM for id='{}'", id.as_str());
        }

        if let Some(e) = entries.iter().find(|e| e.temperature == temperature) {
            let (initiation, table) =
                load_canonical_dsm_tsv_text(e.tsv, e.initiation, e.orientation)?;
            return Ok(Self::from_source_table(&table, initiation, penalty));
        }

        let lo = entries.iter().filter(|e| e.temperature < temperature).max_by_key(|e| e.temperature);
        let hi = entries.iter().filter(|e| e.temperature > temperature).min_by_key(|e| e.temperature);

        let (lo, hi) = match (lo, hi) {
            (Some(l), Some(h)) => (l, h),
            _ => bail!(
                "Temperature {} is outside the range of canonical DSM '{}'",
                temperature, id.as_str()
            ),
        };

        let (init1, table1) =
            load_canonical_dsm_tsv_text(lo.tsv, lo.initiation, lo.orientation)?;
        let (init2, table2) =
            load_canonical_dsm_tsv_text(hi.tsv, hi.initiation, hi.orientation)?;

        let (initiation, table) = interpolate_tables(
            init1, &table1, init2, &table2,
            [temperature as f64, lo.temperature as f64, hi.temperature as f64],
        );

        Ok(Self::from_source_table(&table, initiation, penalty))
    }

    fn from_source_table(source_table: &DsmTable, initiation: i32, penalty: i32) -> Self {
        let mut table = [0i32; DSM_FLAT_SIZE];
        for q1 in 0..6 {
            for q2 in 0..6 {
                for t1 in 0..6 {
                    for t2 in 0..6 {
                        let t1_orig = unsafe { Base::from_idx(t1) }.complement().idx();
                        let t2_orig = unsafe { Base::from_idx(t2) }.complement().idx();

                        let idx = q1 * 216 + q2 * 36 + t1 * 6 + t2;

                        let ext_penalty_mult = if q2 == 0 && t2_orig == 0 {
                            0
                        } else if q2 == 0 || t2_orig == 0 {
                            1
                        } else if q1 == 0 && t1_orig == 0 {
                            if unsafe { Base::from_idx(q2) }
                                .pair_type(unsafe { Base::from_idx(t2_orig) }.complement())
                                .is_match(true)
                            {
                                2
                            } else {
                                0
                            }
                        } else {
                            2
                        };

                        table[idx] = source_table[q1][q2][t1_orig][t2_orig]
                            - penalty * ext_penalty_mult;
                    }
                }
            }
        }

        for (i, &val) in table.iter().enumerate() {
            assert!(
                (val as i64).unsigned_abs() <= crate::dp::MAX_ENERGY as u64,
                "table entry {} = {} exceeds MAX_ENERGY bound {}",
                i,
                val,
                crate::dp::MAX_ENERGY,
            );
        }

        Self { table, initiation }
    }

    pub(crate) fn to_energy(&self, units: i64) -> Energy {
        Energy::from((units as f64 - self.initiation as f64) / -Energy::SCALE)
    }

    /// Check if two bases form a valid seed pair in transformed target space.
    #[inline(always)]
    pub fn seed_pair(q: Base, t: Base, allow_wobble: bool) -> bool {
        q.pair_type(t.complement()).is_match(allow_wobble)
    }

    /// Produce a left-canonical (transposed) copy: `[q1][q2][t1][t2] → [q2][q1][t2][t1]`.
    /// Get raw pointer to the internal flat table
    #[inline(always)]
    pub fn table_ptr(&self) -> *const i32 {
        self.table.as_ptr()
    }

    /// Produce a left-canonical (transposed) copy
    pub fn transpose(&self) -> ScoringModel {
        let mut table = [0i32; DSM_FLAT_SIZE];
        for q1 in 0..BASE_COUNT {
            for q2 in 0..BASE_COUNT {
                for t1 in 0..BASE_COUNT {
                    for t2 in 0..BASE_COUNT {
                        let dst = q1 * 216 + q2 * 36 + t1 * 6 + t2;
                        let src = q2 * 216 + q1 * 36 + t2 * 6 + t1;
                        table[dst] = self.table[src];
                    }
                }
            }
        }
        Self {
            table,
            initiation: self.initiation,
        }
    }

    /// Point query for one transition energy in the canonical 4-base tensor.
    #[inline(always)]
    pub fn transition_energy(&self, q1: u8, q2: u8, t1: u8, t2: u8) -> i32 {
        debug_assert!(q1 < 6 && q2 < 6 && t1 < 6 && t2 < 6);
        let idx = (q1 as usize) * 216 + (q2 as usize) * 36 + (t1 as usize) * 6 + (t2 as usize);
        // SAFETY: All args ∈ 0..6. Max idx = 5*216+5*36+5*6+5 = 1295 < 1296.
        unsafe { *self.table.get_unchecked(idx) }
    }

    /// Point query for one transition energy using semantic bases.
    #[inline(always)]
    pub fn transition_energy_bases(&self, q1: Base, q2: Base, t1: Base, t2: Base) -> i32 {
        self.transition_energy(q1 as u8, q2 as u8, t1 as u8, t2 as u8)
    }

    /// Check if two bases form a valid pair.
    #[inline(always)]
    pub fn is_pair(&self, q: Base, t: Base) -> bool {
        q.pair_type(t.complement()).is_match(true)
    }

    /// Seed energy calculation with antiparallel indexing.
    pub fn energy(
        &self,
        query: &[Base],
        target: &[Base],
        q_pos: usize,
        t_pos: usize,
        len: usize,
    ) -> i32 {
        if len <= 1 {
            return 0;
        }
        let mut score = 0;
        let t_match_end = t_pos + len - 1;
        for i in 0..(len - 1) {
            score += self.transition_energy_bases(
                query[q_pos + i],
                query[q_pos + i + 1],
                target[t_match_end - i],
                target[t_match_end - i - 1],
            );
        }
        score
    }
}


fn load_canonical_dsm_tsv_text(
    text: &str,
    initiation: f64,
    orientation: Orientation,
) -> Result<(i32, DsmTable)> {

    let mut lines = text.lines();
    let header = lines.next().context("DSM TSV is empty")?;
    let header_cols = header.split_whitespace().collect::<Vec<_>>();
    if header_cols != CANONICAL_DSM_HEADER {
        bail!(
            "Invalid DSM TSV header {:?}, expected {:?}",
            header_cols,
            CANONICAL_DSM_HEADER
        );
    }

    let reverse_swap = matches!(orientation, Orientation::ReverseSwap);

    let mut table = [[[[0i32; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];
    let mut seen = [false; DSM_FLAT_SIZE];
    let mut row_count = 0usize;
    for (line_no, line) in lines.enumerate() {
        let cols = line.split_whitespace().collect::<Vec<_>>();
        if cols.len() != 5 {
            bail!(
                "Invalid DSM TSV row {}: expected 5 columns, got {}",
                line_no + 2,
                cols.len()
            );
        }

        let mut q1 = parse_canonical_base(cols[0])
            .with_context(|| format!("Invalid q1 at DSM TSV row {}", line_no + 2))?;
        let mut q2 = parse_canonical_base(cols[1])
            .with_context(|| format!("Invalid q2 at DSM TSV row {}", line_no + 2))?;
        let mut t1 = parse_canonical_base(cols[2])
            .with_context(|| format!("Invalid t1 at DSM TSV row {}", line_no + 2))?;
        let mut t2 = parse_canonical_base(cols[3])
            .with_context(|| format!("Invalid t2 at DSM TSV row {}", line_no + 2))?;
        if reverse_swap {
            (q1, q2, t1, t2) = (t2, t1, q2, q1);
        }

        let delta_g = cols[4]
            .parse::<f64>()
            .with_context(|| format!("Invalid delta_g at DSM TSV row {}", line_no + 2))?;
        if !delta_g.is_finite() {
            bail!("Non-finite delta_g at DSM TSV row {}", line_no + 2);
        }
        if delta_g.abs() > MAX_TSV_ENERGY {
            bail!(
                "DSM value {} at row {} exceeds max energy {}",
                delta_g,
                line_no + 2,
                MAX_TSV_ENERGY
            );
        }

        let idx = q1 * 216 + q2 * 36 + t1 * 6 + t2;
        if std::mem::replace(&mut seen[idx], true) {
            bail!("Duplicate DSM coordinate at row {}", line_no + 2);
        }
        table[q1][q2][t1][t2] = kcal_to_units(-delta_g);
        row_count += 1;
    }

    if row_count != DSM_FLAT_SIZE {
        bail!(
            "Read {} DSM rows from canonical TSV, expected {}",
            row_count,
            DSM_FLAT_SIZE
        );
    }

    Ok((kcal_to_units(initiation), table))
}

fn parse_canonical_base(raw: &str) -> Result<usize> {
    match raw {
        "A" => Ok(Base::A.idx()),
        "C" => Ok(Base::C.idx()),
        "G" => Ok(Base::G.idx()),
        "U" => Ok(Base::U.idx()),
        "N" => Ok(Base::N.idx()),
        "-" => Ok(Base::Gap.idx()),
        other => bail!("Invalid canonical DSM base '{}'", other),
    }
}

fn validate_canonical_manifest_text(text: &str, data_root: &Path) -> Result<()> {
    let doc = text
        .parse::<toml_edit::DocumentMut>()
        .context("Failed to parse DSM manifest")?;
    let entries = doc["dsm"]
        .as_array_of_tables()
        .context("DSM manifest must contain [[dsm]] entries")?;
    if entries.is_empty() {
        bail!("DSM manifest has no entries");
    }

    let mut keys = HashSet::new();
    for entry in entries {
        let id = required_table_str(entry, "id")?;
        let query = required_table_str(entry, "query")?;
        let target = required_table_str(entry, "target")?;
        validate_sequence_type(query).with_context(|| format!("Invalid query for DSM '{}'", id))?;
        validate_sequence_type(target)
            .with_context(|| format!("Invalid target for DSM '{}'", id))?;
        if !keys.insert((query.to_owned(), target.to_owned(), id.to_owned())) {
            bail!("Duplicate DSM manifest entry for ({query}, {target}, {id})");
        }

        required_table_str(entry, "family")?;
        required_table_str(entry, "publication")?;
        required_table_str(entry, "doi")?;
        let orientation = parse_orientation(required_table_str(entry, "orientation")?)
            .with_context(|| format!("Invalid orientation for DSM '{}'", id))?;
        let default_temperature = required_table_int(entry, "default_temperature")?;
        let temperatures = entry
            .get("temperatures")
            .and_then(toml_edit::Item::as_array)
            .context("DSM manifest entry must contain temperatures array")?;
        if temperatures.is_empty() {
            bail!("DSM '{}' has no temperatures", id);
        }

        let mut seen = HashSet::new();
        let mut previous = None;
        let mut has_default = false;
        for value in temperatures.iter() {
            let temp = value
                .as_inline_table()
                .context("DSM temperature entry must be an inline table")?;
            let temperature = required_inline_int(temp, "temperature")?;
            if !seen.insert(temperature) {
                bail!("DSM '{}' has duplicate temperature {}", id, temperature);
            }
            if let Some(prev) = previous {
                if temperature <= prev {
                    bail!("DSM '{}' temperatures must be strictly ascending", id);
                }
            }
            previous = Some(temperature);
            has_default |= temperature == default_temperature;

            let file = required_inline_str(temp, "file")?;
            let initiation = required_inline_number(temp, "initiation_kcal")?;
            if !initiation.is_finite() {
                bail!("DSM '{}' has non-finite initiation", id);
            }
            let path = data_root.join(file);
            let table_text = std::fs::read_to_string(&path)
                .with_context(|| format!("Failed to read canonical DSM table {}", path.display()))?;
            load_canonical_dsm_tsv_text(&table_text, initiation, orientation)
                .with_context(|| format!("Invalid canonical DSM table {}", path.display()))?;
        }
        if !has_default {
            bail!(
                "DSM '{}' default temperature {} is not present",
                id, default_temperature
            );
        }
    }

    Ok(())
}

fn parse_orientation(raw: &str) -> Result<Orientation> {
    match raw {
        "identity" => Ok(Orientation::Identity),
        "reverse-swap" => Ok(Orientation::ReverseSwap),
        other => bail!("unknown orientation '{}'", other),
    }
}

fn validate_sequence_type(raw: &str) -> Result<SequenceType> {
    SequenceType::try_from(raw).map_err(|e| anyhow::anyhow!(e))
}

fn required_table_str<'a>(table: &'a toml_edit::Table, key: &str) -> Result<&'a str> {
    table
        .get(key)
        .and_then(toml_edit::Item::as_str)
        .with_context(|| format!("DSM manifest entry requires string field '{}'", key))
}

fn required_table_int(table: &toml_edit::Table, key: &str) -> Result<i64> {
    table
        .get(key)
        .and_then(toml_edit::Item::as_integer)
        .with_context(|| format!("DSM manifest entry requires integer field '{}'", key))
}

fn required_inline_str<'a>(table: &'a toml_edit::InlineTable, key: &str) -> Result<&'a str> {
    table
        .get(key)
        .and_then(toml_edit::Value::as_str)
        .with_context(|| format!("DSM temperature entry requires string field '{}'", key))
}

fn required_inline_int(table: &toml_edit::InlineTable, key: &str) -> Result<i64> {
    table
        .get(key)
        .and_then(toml_edit::Value::as_integer)
        .with_context(|| format!("DSM temperature entry requires integer field '{}'", key))
}

fn required_inline_number(table: &toml_edit::InlineTable, key: &str) -> Result<f64> {
    let value = table
        .get(key)
        .with_context(|| format!("DSM temperature entry requires numeric field '{}'", key))?;
    value
        .as_float()
        .or_else(|| value.as_integer().map(|integer| integer as f64))
        .with_context(|| format!("DSM temperature entry field '{}' must be numeric", key))
}

fn interpolate_tables(
    offset_1: i32,
    table_1: &DsmTable,
    offset_2: i32,
    table_2: &DsmTable,
    temps: [f64; 3],
) -> (i32, DsmTable) {
    let mut table = [[[[0i32; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];
    for q1 in 0..BASE_COUNT {
        for q2 in 0..BASE_COUNT {
            for t1 in 0..BASE_COUNT {
                for t2 in 0..BASE_COUNT {
                    table[q1][q2][t1][t2] =
                        lerp(table_1[q1][q2][t1][t2], table_2[q1][q2][t1][t2], temps);
                }
            }
        }
    }
    (lerp(offset_1, offset_2, temps), table)
}

fn lerp(v1: i32, v2: i32, temps: [f64; 3]) -> i32 {
    let [t0, t1, t2] = temps;
    let diff = i64::from(v1) - i64::from(v2);
    ((t0 - t2) / (t1 - t2) * diff as f64 + f64::from(v2)).round() as i32
}

fn kcal_to_units(value: f64) -> i32 {
    (value * Energy::SCALE).round() as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pairing_distinguishes_seed_and_extension_modes() {
        let model = ScoringModel::from_canonical(DsmId::T04, 37, Energy::from(0.0)).unwrap();

        // DP/display pairing follows the scoring matrix semantics.
        assert!(model.is_pair(Base::G, Base::G));
        assert!(model.is_pair(Base::G, Base::A));

        // Seed pairing can be stricter than the extension model.
        assert!(!ScoringModel::seed_pair(Base::G, Base::A, false));
        assert!(ScoringModel::seed_pair(Base::G, Base::A, true));
    }

    #[test]
    fn transition_energy_returns_nonzero_for_valid_pairs() {
        let model = ScoringModel::from_canonical(DsmId::T04, 37, Energy::from(0.0)).unwrap();
        // Original target was U-A, index space is A-U
        let energy = model.transition_energy_bases(Base::A, Base::U, Base::A, Base::U);
        let gap_energy = model.transition_energy_bases(
            Base::Gap,
            Base::A,
            Base::Gap,
            Base::A, // original was U, index is A
        );
        assert!(
            gap_energy != 0 || energy != 0,
            "At least one transition query should be non-zero"
        );
    }

    #[test]
    fn gg_cc_stack_is_strongest() {
        let model = ScoringModel::from_canonical(DsmId::T04, 37, Energy::from(0.0)).unwrap();
        let actual = model.transition_energy_bases(Base::G, Base::G, Base::G, Base::G);
        assert_eq!(actual, 33016);
    }

    #[test]
    fn energy_conversion_roundtrips() {
        let model = ScoringModel::from_canonical(DsmId::T04, 37, Energy::from(0.0)).unwrap();
        let energy = model.to_energy(33016_i64);
        let kcal = f64::from(energy);
        assert!((kcal - 2.8264).abs() < 0.001);
    }

    #[test]
    fn transpose_swaps_both_pairs() {
        let right = ScoringModel::from_canonical(DsmId::T04, 37, Energy::from(0.005)).unwrap();
        let left = right.transpose();
        for q1 in 0u8..6 {
            for q2 in 0u8..6 {
                for t1 in 0u8..6 {
                    for t2 in 0u8..6 {
                        assert_eq!(
                            left.transition_energy(q1, q2, t1, t2),
                            right.transition_energy(q2, q1, t2, t1),
                            "transpose mismatch at ({},{},{},{})",
                            q1,
                            q2,
                            t1,
                            t2
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn canonical_dsm_data_package_is_valid() {
        let data_root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("data/dsm");
        validate_canonical_manifest_text(
            include_str!("../data/dsm/manifest.toml"),
            data_root.as_path(),
        )
        .unwrap();
    }

    #[test]
    fn canonical_dsm_parser_rejects_bad_header() {
        let tsv = "q1\tq2\tt1\tt2\tenergy\n";
        assert!(load_canonical_dsm_tsv_text(tsv, 1.0, Orientation::Identity).is_err());
    }

    #[test]
    fn canonical_dsm_parser_rejects_duplicate_coordinate() {
        let tsv = concat!(
            "q1\tq2\tt1\tt2\tdelta_g_kcal_per_mol\n",
            "A\tA\tA\tA\t1.0\n",
            "A\tA\tA\tA\t2.0\n",
        );
        assert!(load_canonical_dsm_tsv_text(tsv, 1.0, Orientation::Identity).is_err());
    }

}
