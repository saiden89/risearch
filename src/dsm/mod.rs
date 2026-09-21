#![allow(dead_code)]
//! RNA dinucleotide stacking energy matrices (DSM)
//!
//! Tables encode nearest-neighbor thermodynamic parameters for RNA-RNA interactions.
//! Canonical table values are stored in kcal/mol and converted to score units on load.
//!
//! Indices: `[q1][q2][t1][t2]` = stacking energy for:
//! ```text
//! Query:  5'─ q1 ─ q2 ─ 3'
//!             |    |
//! Target: 3'─ t1 ─ t2 ─ 5'
//! ```

use anyhow::{bail, Context};
use std::collections::HashSet;
use std::path::Path;

use crate::error::{Error, Result};
use crate::types::{DsmId, Energy, SequenceType, BASE_COUNT};

mod model;
mod parse;

use parse::load_dsm_tsv_text;

pub use model::ScoringModel;

pub(crate) type DsmTable = [[[[i32; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];

pub(crate) const DSM_FLAT_SIZE: usize = BASE_COUNT * BASE_COUNT * BASE_COUNT * BASE_COUNT;

pub(crate) const DSM_HEADER: [&str; 5] = ["q1", "q2", "t1", "t2", "delta_g_kcal_per_mol"];

#[inline(always)]
pub(crate) const fn flat_idx(q1: u8, q2: u8, t1: u8, t2: u8) -> usize {
    (q1 as usize) * BASE_COUNT.pow(3)
        + (q2 as usize) * BASE_COUNT.pow(2)
        + (t1 as usize) * BASE_COUNT
        + t2 as usize
}

#[derive(Clone, Copy)]
pub(crate) enum Orientation {
    Identity,
    ReverseSwap,
}

struct Dsm {
    id: &'static str,
    temperature: i32,
    initiation: f64,
    invalid_transition: f64,
    orientation: Orientation,
    tsv_content: &'static str,
}

include!(concat!(env!("OUT_DIR"), "/generated_canonical_tables.rs"));

/// Service for managing and loading dinucleotide stacking models (DSM).
pub struct DsmRegistry;

impl DsmRegistry {
    /// List all available DSM string identifiers (for CLI/UI).
    pub const fn all_names() -> &'static [&'static str] {
        BUILTIN_NAMES
    }

    /// Map a string identifier to a DsmId.
    pub(crate) fn parse_id(s: &str) -> Result<DsmId> {
        if BUILTIN_TABLES.iter().any(|e| e.id == s) {
            Ok(DsmId(s.to_string()))
        } else {
            Err(Error::Dsm(format!("unknown DSM id '{s}'")))
        }
    }

    /// Load a bundled DSM table, interpolating between bracket temperatures if needed.
    pub fn load(id: &DsmId, temperature: i32) -> Result<(Energy, DsmTable)> {
        let entries: Vec<&Dsm> = BUILTIN_TABLES.iter().filter(|e| e.id == id.0).collect();

        if entries.is_empty() {
            return Err(Error::Dsm(format!("No bundled DSM for id='{}'", id.0)));
        }

        if let Some(e) = entries.iter().find(|e| e.temperature == temperature) {
            return load_dsm_tsv_text(
                e.tsv_content,
                e.initiation,
                e.invalid_transition,
                e.orientation,
            );
        }

        let lo = entries
            .iter()
            .filter(|e| e.temperature < temperature)
            .max_by_key(|e| e.temperature);
        let hi = entries
            .iter()
            .filter(|e| e.temperature > temperature)
            .min_by_key(|e| e.temperature);

        let (lo, hi) = match (lo, hi) {
            (Some(l), Some(h)) => (l, h),
            _ => {
                return Err(Error::Dsm(format!(
                    "Temperature {} is outside range for bundled DSM '{}'",
                    temperature, id.0
                )))
            }
        };

        let (init1, table1) = load_dsm_tsv_text(
            lo.tsv_content,
            lo.initiation,
            lo.invalid_transition,
            lo.orientation,
        )?;
        let (init2, table2) = load_dsm_tsv_text(
            hi.tsv_content,
            hi.initiation,
            hi.invalid_transition,
            hi.orientation,
        )?;

        Ok(interpolate_tables(
            init1,
            &table1,
            init2,
            &table2,
            [
                temperature as f64,
                lo.temperature as f64,
                hi.temperature as f64,
            ],
        ))
    }
}

fn validate_canonical_manifest_text(text: &str, data_root: &Path) -> anyhow::Result<()> {
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
        let invalid_transition = required_table_number(entry, "invalid_transition_kcal")?;
        if !invalid_transition.is_finite() {
            bail!("DSM '{}' has non-finite invalid_transition_kcal", id);
        }
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
            let table_text = fs_err::read_to_string(&path).with_context(|| {
                format!("Failed to read canonical DSM table {}", path.display())
            })?;
            load_dsm_tsv_text(&table_text, initiation, invalid_transition, orientation)
                .with_context(|| format!("Invalid DSM table {}", path.display()))?;
        }
        if !has_default {
            bail!(
                "DSM '{}' default temperature {} is not present",
                id,
                default_temperature
            );
        }
    }

    Ok(())
}

fn parse_orientation(raw: &str) -> anyhow::Result<Orientation> {
    match raw {
        "identity" => Ok(Orientation::Identity),
        "reverse-swap" => Ok(Orientation::ReverseSwap),
        other => bail!("unknown orientation '{}'", other),
    }
}

fn validate_sequence_type(raw: &str) -> anyhow::Result<SequenceType> {
    SequenceType::try_from(raw).map_err(|e| anyhow::anyhow!(e))
}

fn required_table_str<'a>(table: &'a toml_edit::Table, key: &str) -> anyhow::Result<&'a str> {
    table
        .get(key)
        .and_then(toml_edit::Item::as_str)
        .with_context(|| format!("DSM manifest entry requires string field '{}'", key))
}

fn required_table_int(table: &toml_edit::Table, key: &str) -> anyhow::Result<i64> {
    table
        .get(key)
        .and_then(toml_edit::Item::as_integer)
        .with_context(|| format!("DSM manifest entry requires integer field '{}'", key))
}

fn required_table_number(table: &toml_edit::Table, key: &str) -> anyhow::Result<f64> {
    let value = table
        .get(key)
        .with_context(|| format!("DSM manifest entry requires numeric field '{}'", key))?;
    value
        .as_float()
        .or_else(|| value.as_integer().map(|integer| integer as f64))
        .with_context(|| format!("DSM manifest entry field '{}' must be numeric", key))
}

fn required_inline_str<'a>(
    table: &'a toml_edit::InlineTable,
    key: &str,
) -> anyhow::Result<&'a str> {
    table
        .get(key)
        .and_then(toml_edit::Value::as_str)
        .with_context(|| format!("DSM temperature entry requires string field '{}'", key))
}

fn required_inline_int(table: &toml_edit::InlineTable, key: &str) -> anyhow::Result<i64> {
    table
        .get(key)
        .and_then(toml_edit::Value::as_integer)
        .with_context(|| format!("DSM temperature entry requires integer field '{}'", key))
}

fn required_inline_number(table: &toml_edit::InlineTable, key: &str) -> anyhow::Result<f64> {
    let value = table
        .get(key)
        .with_context(|| format!("DSM temperature entry requires numeric field '{}'", key))?;
    value
        .as_float()
        .or_else(|| value.as_integer().map(|integer| integer as f64))
        .with_context(|| format!("DSM temperature entry field '{}' must be numeric", key))
}

fn interpolate_tables(
    offset_1: Energy,
    table_1: &DsmTable,
    offset_2: Energy,
    table_2: &DsmTable,
    temps: [f64; 3],
) -> (Energy, DsmTable) {
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
    (Energy(lerp(offset_1.0, offset_2.0, temps)), table)
}

fn lerp(v1: i32, v2: i32, temps: [f64; 3]) -> i32 {
    let [t0, t1, t2] = temps;
    let diff = i64::from(v1) - i64::from(v2);
    ((t0 - t2) / (t1 - t2) * diff as f64 + f64::from(v2)).round() as i32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Base;

    type ManifestEntry = (String, Orientation, f64, Vec<(i32, String, f64)>);

    fn canonical_manifest() -> Vec<ManifestEntry> {
        let doc = include_str!("../../data/dsm/manifest.toml")
            .parse::<toml_edit::DocumentMut>()
            .expect("canonical DSM manifest parses");
        doc["dsm"]
            .as_array_of_tables()
            .expect("canonical DSM manifest has dsm entries")
            .iter()
            .map(|entry| {
                let id = entry["id"].as_str().expect("DSM id").to_owned();
                let orientation = match entry["orientation"].as_str().expect("orientation") {
                    "identity" => Orientation::Identity,
                    "reverse-swap" => Orientation::ReverseSwap,
                    other => panic!("unknown canonical orientation {other}"),
                };
                let invalid = entry["invalid_transition_kcal"]
                    .as_float()
                    .or_else(|| {
                        entry["invalid_transition_kcal"]
                            .as_integer()
                            .map(|v| v as f64)
                    })
                    .expect("invalid transition");
                let temperatures = entry["temperatures"]
                    .as_array()
                    .expect("DSM temperatures")
                    .iter()
                    .map(|value| {
                        let table = value.as_inline_table().expect("inline temperature table");
                        let temperature =
                            table["temperature"].as_integer().expect("temperature") as i32;
                        let file = table["file"].as_str().expect("TSV file").to_owned();
                        let initiation = table["initiation_kcal"]
                            .as_float()
                            .or_else(|| table["initiation_kcal"].as_integer().map(|v| v as f64))
                            .expect("initiation");
                        (temperature, file, initiation)
                    })
                    .collect();
                (id, orientation, invalid, temperatures)
            })
            .collect()
    }

    #[test]
    fn pairing_policy_uses_physical_duplex_bases() {
        assert!(Base::A.pair_type(Base::U).is_match(false));
        assert!(Base::G.pair_type(Base::U).is_match(true));
        assert!(!Base::G.pair_type(Base::U).is_match(false));
        assert!(!Base::G.pair_type(Base::A).is_match(true));
    }

    #[test]
    fn score_returns_nonzero_for_valid_pairs() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
        let score = model.score_bases(Base::A, Base::C, Base::U, Base::G);
        let gap_score = model.score_bases(Base::Gap, Base::A, Base::Gap, Base::U);
        assert!(
            gap_score != 0 || score != 0,
            "At least one transition query should be non-zero"
        );
    }

    #[test]
    fn gg_cc_stack_is_strongest() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
        let actual = model.score_bases(Base::G, Base::G, Base::C, Base::C);
        assert_eq!(actual, 33016);
    }

    #[test]
    fn bundled_models_have_explicit_stack_and_initiation_energies() {
        // Values read from the canonical TSV/manifest, in integer score units.
        // An AC/UG stack is asymmetric in the mixed-strand model.
        for (name, init, ac_ug, gg_cc) in [
            ("t04", 61280, 21805, 33016),
            ("t99", 55900, 22000, 33000),
            ("slh04", 20438, 13852, 18007),
            ("s95-rna-dna", 40350, 10714, 21010),
        ] {
            let (actual_init, table) = DsmRegistry::load(&DsmId::from(name), 37).unwrap();
            assert_eq!(actual_init, Energy(init), "{name}");
            assert_eq!(table[1][2][5][3], ac_ug, "{name}");
            assert_eq!(table[3][3][2][2], gg_cc, "{name}");
            assert_eq!(table[0][0][0][0], -200000, "omitted transition in {name}");
        }
        let (init, reversed) = DsmRegistry::load(&DsmId::from("s95-dna-rna"), 37).unwrap();
        assert_eq!(init, Energy(40350));
        // AC/UG (RNA/DNA) becomes GU/CA (DNA/RNA), not a plain transpose.
        assert_eq!(reversed[3][5][2][1], 10714);
    }

    #[test]
    fn temperature_interpolation_uses_bracketing_integer_scores() {
        let (init, table) = DsmRegistry::load(&DsmId::from("t04"), 31).unwrap();
        // Midpoint of 25C and 37C: initiation (63786 + 61280)/2;
        // AC/UG stack (25373 + 21805)/2; GG/CC (36920 + 33016)/2.
        assert_eq!(init, Energy(62533));
        assert_eq!(table[1][2][5][3], 23589);
        assert_eq!(table[3][3][2][2], 34968);
        assert_eq!(table[1][0][5][0], 5000);
        assert!(DsmRegistry::load(&DsmId::from("t99"), 31).is_err());
        assert!(DsmRegistry::load(&DsmId::from("t04"), 51).is_err());
    }

    #[test]
    fn interpolation_rounds_exact_half_away_from_zero() {
        // temps: [target, v1's knot, v2's knot]; 31C sits exactly midway in 25..37.
        assert_eq!(lerp(1, 2, [31.0, 25.0, 37.0]), 2);
        assert_eq!(lerp(-1, -2, [31.0, 25.0, 37.0]), -2);
        assert_eq!(lerp(1_200_000, 1_200_001, [31.0, 25.0, 37.0]), 1_200_001);
    }

    #[test]
    fn energy_conversion_roundtrips() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let model = ScoringModel::new(&source, init, Energy::from_kcal(0.0));
        let energy = model.binding_energy(Energy(33016));
        assert!((energy.to_kcal() - 2.8264).abs() < 0.001);
    }

    #[test]
    fn transpose_swaps_both_pairs() {
        let (init, source) = DsmRegistry::load(&DsmId::from("t04"), 37).unwrap();
        let right = ScoringModel::new(&source, init, Energy::from_kcal(0.005));
        let left = right.transpose();
        for q1 in 0u8..6 {
            for q2 in 0u8..6 {
                for t1 in 0u8..6 {
                    for t2 in 0u8..6 {
                        assert_eq!(
                            left.score(q1, q2, t1, t2),
                            right.score(q2, q1, t2, t1),
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
            include_str!("../../data/dsm/manifest.toml"),
            data_root.as_path(),
        )
        .unwrap();
    }

    #[test]
    fn every_bundled_model_loads_at_its_source_temperatures() {
        let entries = canonical_manifest();
        let mut file_paths = HashSet::new();
        let mut ids = HashSet::new();
        let mut source_count = 0;

        for (id, _, _, temperatures) in &entries {
            ids.insert(id.clone());
            for (temperature, file, initiation) in temperatures {
                source_count += 1;
                file_paths.insert(file.clone());
                let (actual_init, _) = DsmRegistry::load(&DsmId::from(id.as_str()), *temperature)
                    .unwrap_or_else(|err| panic!("loading {id} at {temperature}C: {err}"));
                assert_eq!(
                    actual_init,
                    Energy::from_kcal(*initiation),
                    "initiation {id} at {temperature}C"
                );
            }
        }

        assert_eq!(file_paths.len(), 16, "all canonical TSV files are covered");
        assert_eq!(
            source_count, 21,
            "all model-temperature source entries are covered"
        );
        assert_eq!(ids.len(), 5, "all bundled DSM identifiers are covered");
        assert_eq!(DsmRegistry::all_names().len(), ids.len());
    }

    #[test]
    fn canonical_dsm_parser_rejects_bad_header() {
        let tsv = "q1\tq2\tt1\tt2\tenergy\n";
        assert!(load_dsm_tsv_text(tsv, 1.0, 20.0, Orientation::Identity).is_err());
    }

    #[test]
    fn canonical_dsm_parser_rejects_duplicate_coordinate() {
        let tsv = concat!(
            "q1\tq2\tt1\tt2\tdelta_g_kcal_per_mol\n",
            "A\tA\tA\tA\t1.0\n",
            "A\tA\tA\tA\t2.0\n",
        );
        assert!(load_dsm_tsv_text(tsv, 1.0, 20.0, Orientation::Identity).is_err());
    }

    #[test]
    fn manifest_numeric_fields_must_parse() {
        let data_root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("data/dsm");
        let manifest = |transition: &str, initiation: &str| {
            format!(
                r#"[[dsm]]
id = "t04"
family = "turner-2004"
query = "rna"
target = "rna"
publication = "p"
doi = "d"
invalid_transition_kcal = {transition}
orientation = "identity"
default_temperature = 37
temperatures = [{{ temperature = 37, file = "t04/37.tsv", initiation_kcal = {initiation} }}]
"#
            )
        };

        validate_canonical_manifest_text(&manifest("20", "6.128"), data_root.as_path()).unwrap();
        assert!(validate_canonical_manifest_text(
            &manifest("\"twenty\"", "6.128"),
            data_root.as_path()
        )
        .is_err());
        assert!(
            validate_canonical_manifest_text(&manifest("20", "\"six\""), data_root.as_path())
                .is_err()
        );
    }

    fn table_hash(table: &DsmTable) -> u64 {
        table
            .iter()
            .flatten()
            .flatten()
            .flatten()
            .fold(0u64, |h, &v| {
                h.wrapping_mul(31).wrapping_add(v as u32 as u64)
            })
    }

    /// Pin every cell of every bundled model's 6^4 tensor at its default temperature.
    ///
    /// A single flipped cell changes the hash. Regenerate after an intentional
    /// table update with `cargo test bundled_model_tables -- --ignored --nocapture`.
    #[test]
    fn bundled_model_tables_are_bitwise_stable() {
        for (name, expected) in [
            ("t04", 0x3F7C_035F_D32B_A9BC_u64),
            ("t99", 0x405F_24E9_8C52_81C8),
            ("slh04", 0x1D57_D0DA_DC39_6F2E),
            ("s95-rna-dna", 0xA32E_410C_EEC0_4931),
            ("s95-dna-rna", 0xC15F_83D7_9D7C_AB3F),
        ] {
            let (_, table) = DsmRegistry::load(&DsmId::from(name), 37).unwrap();
            assert_eq!(
                table_hash(&table),
                expected,
                "{name} table changed — if intentional, update the hash",
            );
        }
    }

    #[test]
    fn unsupported_model_temperatures_are_rejected() {
        for (id, _, _, temperatures) in canonical_manifest() {
            let first = temperatures.first().unwrap().0;
            let last = temperatures.last().unwrap().0;
            for temperature in [first - 1, last + 1, -1, 101] {
                assert!(
                    DsmRegistry::load(&DsmId::from(id.as_str()), temperature).is_err(),
                    "expected {id} at {temperature}C to be rejected"
                );
            }
        }
        assert!(DsmRegistry::load(&DsmId::from("missing"), 37).is_err());
    }
}
