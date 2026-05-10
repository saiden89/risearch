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

use anyhow::{bail, Context, Result};
use std::collections::HashSet;
use std::path::Path;

use crate::types::{DsmId, Energy, SequenceType, BASE_COUNT};

mod model;
mod parse;

use parse::load_canonical_dsm_tsv_text;

pub use model::ScoringModel;

pub(crate) type DsmTable = [[[[i32; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];

pub(crate) const DSM_FLAT_SIZE: usize = BASE_COUNT * BASE_COUNT * BASE_COUNT * BASE_COUNT;

pub(crate) const MAX_TSV_ENERGY: f64 = 20.0;

pub const DSM_IDS: [&str; DsmId::ALL.len()] = {
    let mut ids = [""; DsmId::ALL.len()];
    let mut i = 0;
    while i < DsmId::ALL.len() {
        ids[i] = DsmId::ALL[i].as_str();
        i += 1;
    }
    ids
};
pub(crate) const CANONICAL_DSM_HEADER: [&str; 5] = ["q1", "q2", "t1", "t2", "delta_g_kcal_per_mol"];

#[derive(Clone, Copy)]
pub(crate) enum Orientation {
    Identity,
    ReverseSwap,
}

struct CanonicalDsms {
    id: DsmId,
    temperature: i32,
    initiation: f64,
    invalid_transition: f64,
    orientation: Orientation,
    tsv: &'static str,
}

static CANONICAL_TABLES: &[CanonicalDsms] = &[
    // t04: Turner 2004 RNA-RNA
    CanonicalDsms {
        id: DsmId::T04,
        temperature: 0,
        initiation: 6.8982,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/t04/0.tsv"),
    },
    CanonicalDsms {
        id: DsmId::T04,
        temperature: 25,
        initiation: 6.3786,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/t04/25.tsv"),
    },
    CanonicalDsms {
        id: DsmId::T04,
        temperature: 37,
        initiation: 6.128,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/t04/37.tsv"),
    },
    CanonicalDsms {
        id: DsmId::T04,
        temperature: 42,
        initiation: 5.9848,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/t04/42.tsv"),
    },
    CanonicalDsms {
        id: DsmId::T04,
        temperature: 50,
        initiation: 5.7678,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/t04/50.tsv"),
    },
    // slh04: SantaLucia-Hicks 2004 DNA-DNA
    CanonicalDsms {
        id: DsmId::Slh04,
        temperature: 0,
        initiation: 1.858,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/slh04/0.tsv"),
    },
    CanonicalDsms {
        id: DsmId::Slh04,
        temperature: 25,
        initiation: 1.9532,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/slh04/25.tsv"),
    },
    CanonicalDsms {
        id: DsmId::Slh04,
        temperature: 37,
        initiation: 2.0438,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/slh04/37.tsv"),
    },
    CanonicalDsms {
        id: DsmId::Slh04,
        temperature: 42,
        initiation: 2.1748,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/slh04/42.tsv"),
    },
    CanonicalDsms {
        id: DsmId::Slh04,
        temperature: 50,
        initiation: 2.4016,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/slh04/50.tsv"),
    },
    // s95-rna-dna: Sugimoto 1995 RNA-DNA (identity)
    CanonicalDsms {
        id: DsmId::S95RnaDna,
        temperature: 0,
        initiation: 4.3476,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/s95/0.tsv"),
    },
    CanonicalDsms {
        id: DsmId::S95RnaDna,
        temperature: 25,
        initiation: 4.1588,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/s95/25.tsv"),
    },
    CanonicalDsms {
        id: DsmId::S95RnaDna,
        temperature: 37,
        initiation: 4.035,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/s95/37.tsv"),
    },
    CanonicalDsms {
        id: DsmId::S95RnaDna,
        temperature: 42,
        initiation: 3.9402,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/s95/42.tsv"),
    },
    CanonicalDsms {
        id: DsmId::S95RnaDna,
        temperature: 50,
        initiation: 3.8672,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/s95/50.tsv"),
    },
    // s95-dna-rna: Sugimoto 1995 DNA-RNA (reverse-swap of RNA-DNA)
    CanonicalDsms {
        id: DsmId::S95DnaRna,
        temperature: 0,
        initiation: 4.3476,
        invalid_transition: 20.0,
        orientation: Orientation::ReverseSwap,
        tsv: include_str!("../../data/dsm/s95/0.tsv"),
    },
    CanonicalDsms {
        id: DsmId::S95DnaRna,
        temperature: 25,
        initiation: 4.1588,
        invalid_transition: 20.0,
        orientation: Orientation::ReverseSwap,
        tsv: include_str!("../../data/dsm/s95/25.tsv"),
    },
    CanonicalDsms {
        id: DsmId::S95DnaRna,
        temperature: 37,
        initiation: 4.035,
        invalid_transition: 20.0,
        orientation: Orientation::ReverseSwap,
        tsv: include_str!("../../data/dsm/s95/37.tsv"),
    },
    CanonicalDsms {
        id: DsmId::S95DnaRna,
        temperature: 42,
        initiation: 3.9402,
        invalid_transition: 20.0,
        orientation: Orientation::ReverseSwap,
        tsv: include_str!("../../data/dsm/s95/42.tsv"),
    },
    CanonicalDsms {
        id: DsmId::S95DnaRna,
        temperature: 50,
        initiation: 3.8672,
        invalid_transition: 20.0,
        orientation: Orientation::ReverseSwap,
        tsv: include_str!("../../data/dsm/s95/50.tsv"),
    },
    // t99: Turner 1999 RNA-RNA (37°C only)
    CanonicalDsms {
        id: DsmId::T99,
        temperature: 37,
        initiation: 5.59,
        invalid_transition: 20.0,
        orientation: Orientation::Identity,
        tsv: include_str!("../../data/dsm/t99/37.tsv"),
    },
];

/// Gap index used for DSM transition queries (linked to Base::Gap).
pub use crate::types::GAP;

pub(crate) fn load_canonical_table(id: DsmId, temperature: i32) -> Result<(Energy, DsmTable)> {
    let entries: Vec<&CanonicalDsms> = CANONICAL_TABLES.iter().filter(|e| e.id == id).collect();

    if entries.is_empty() {
        bail!("No canonical DSM for id='{}'", id.as_str());
    }

    if let Some(e) = entries.iter().find(|e| e.temperature == temperature) {
        return load_canonical_dsm_tsv_text(
            e.tsv,
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
        _ => bail!(
            "Temperature {} is outside the range of canonical DSM '{}'",
            temperature,
            id.as_str()
        ),
    };

    let (init1, table1) =
        load_canonical_dsm_tsv_text(lo.tsv, lo.initiation, lo.invalid_transition, lo.orientation)?;
    let (init2, table2) =
        load_canonical_dsm_tsv_text(hi.tsv, hi.initiation, hi.invalid_transition, hi.orientation)?;

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
            let table_text = std::fs::read_to_string(&path).with_context(|| {
                format!("Failed to read canonical DSM table {}", path.display())
            })?;
            load_canonical_dsm_tsv_text(&table_text, initiation, invalid_transition, orientation)
                .with_context(|| format!("Invalid canonical DSM table {}", path.display()))?;
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

fn required_table_number(table: &toml_edit::Table, key: &str) -> Result<f64> {
    let value = table
        .get(key)
        .with_context(|| format!("DSM manifest entry requires numeric field '{}'", key))?;
    value
        .as_float()
        .or_else(|| value.as_integer().map(|integer| integer as f64))
        .with_context(|| format!("DSM manifest entry field '{}' must be numeric", key))
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

    #[test]
    fn pairing_distinguishes_seed_and_extension_modes() {
        let model = ScoringModel::from_canonical(DsmId::T04, 37, Energy::from(0.0)).unwrap();

        assert!(model.is_pair(Base::G, Base::G));
        assert!(model.is_pair(Base::G, Base::A));
        assert!(!ScoringModel::seed_pair(Base::G, Base::A, false));
        assert!(ScoringModel::seed_pair(Base::G, Base::A, true));
    }

    #[test]
    fn transition_score_returns_nonzero_for_valid_pairs() {
        let model = ScoringModel::from_canonical(DsmId::T04, 37, Energy::from(0.0)).unwrap();
        let score = model.transition_score_bases(Base::A, Base::U, Base::A, Base::U);
        let gap_score = model.transition_score_bases(Base::Gap, Base::A, Base::Gap, Base::A);
        assert!(
            gap_score != 0 || score != 0,
            "At least one transition query should be non-zero"
        );
    }

    #[test]
    fn gg_cc_stack_is_strongest() {
        let model = ScoringModel::from_canonical(DsmId::T04, 37, Energy::from(0.0)).unwrap();
        let actual = model.transition_score_bases(Base::G, Base::G, Base::G, Base::G);
        assert_eq!(actual, 33016);
    }

    #[test]
    fn energy_conversion_roundtrips() {
        let model = ScoringModel::from_canonical(DsmId::T04, 37, Energy::from(0.0)).unwrap();
        let energy = model.hit_energy(33016, 0, 0, 0);
        assert!((energy.to_kcal() - 2.8264).abs() < 0.001);
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
                            left.transition_score(q1, q2, t1, t2),
                            right.transition_score(q2, q1, t2, t1),
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
    fn canonical_dsm_parser_rejects_bad_header() {
        let tsv = "q1\tq2\tt1\tt2\tenergy\n";
        assert!(load_canonical_dsm_tsv_text(tsv, 1.0, 20.0, Orientation::Identity).is_err());
    }

    #[test]
    fn canonical_dsm_parser_rejects_duplicate_coordinate() {
        let tsv = concat!(
            "q1\tq2\tt1\tt2\tdelta_g_kcal_per_mol\n",
            "A\tA\tA\tA\t1.0\n",
            "A\tA\tA\tA\t2.0\n",
        );
        assert!(load_canonical_dsm_tsv_text(tsv, 1.0, 20.0, Orientation::Identity).is_err());
    }
}
