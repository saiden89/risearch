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
use std::path::{Path, PathBuf};

use crate::config::{Matrix, ScoreConfig};
use crate::types::{Base, Energy, BASE_COUNT};

/// DSM table type: 4D array [q1][q2][t1][t2]
pub(crate) type DsmTable = [[[[i32; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];

/// Number of entries in a flattened DSM table (6^4 = 1296).
pub(crate) const DSM_FLAT_SIZE: usize = BASE_COUNT * BASE_COUNT * BASE_COUNT * BASE_COUNT;

const DSM_TSV_VALUE_COUNT: usize = DSM_FLAT_SIZE + 1;
const BUILTIN_SCALE: i32 = 100;
/// Duplex initiation free energy for built-in Turner 2004/1999 parameters.
/// Applied once per interaction: ΔG = (raw_score - init) / -SCALE.
/// Value: 5.59 kcal/mol × 10000 = 55,900 raw units.
const INITIATION_ENERGY_RAW: i32 = 55_900;
const MAX_TSV_ENERGY: f64 = 20.0;
const DEFAULT_TEMPERATURES: [f64; 3] = [310.15, 310.15, 315.15];
const TSV_BASE_TO_RUST: [usize; BASE_COUNT] = [
    Base::A as usize,
    Base::C as usize,
    Base::G as usize,
    Base::U as usize,
    Base::N as usize,
    Base::Gap as usize,
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
    initiation_raw: i32,
}

impl ScoringModel {
    /// Create a new scoring table from a base matrix, baking in the extension penalty.
    pub fn new(matrix: Matrix, penalty: i32) -> Self {
        let source_table = match matrix {
            Matrix::T04 => &T04,
            Matrix::T99 => &T99,
        };

        Self::from_source_table(source_table, INITIATION_ENERGY_RAW, penalty, BUILTIN_SCALE)
    }

    pub(crate) fn from_score_config(score: &ScoreConfig) -> Result<Self> {
        let penalty = score.penalty.to_raw();
        let Some(matpath) = &score.matpath else {
            return Ok(Self::new(score.matrix, penalty));
        };
        let matrix = matrix_file_stem(score.matrix);
        let temps = parse_temperatures(score.temperature.as_deref())?;
        let (initiation_raw, source_table) =
            load_temperature_table(Path::new(matpath), matrix, score.matrix2.as_deref(), temps)?;
        Ok(Self::from_source_table(
            &source_table,
            initiation_raw,
            penalty,
            1,
        ))
    }

    fn from_source_table(
        source_table: &DsmTable,
        initiation_raw: i32,
        penalty: i32,
        source_scale: i32,
    ) -> Self {
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

                        table[idx] = source_table[q1][q2][t1_orig][t2_orig] * source_scale
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

        Self { table, initiation_raw }
    }

    pub(crate) fn energy_from_raw(&self, raw: i64) -> Energy {
        Energy::from((raw as f64 - self.initiation_raw as f64) / -Energy::RAW_SCALE)
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
        for q1 in 0..6 {
            for q2 in 0..6 {
                for t1 in 0..6 {
                    for t2 in 0..6 {
                        let dst = q1 * 216 + q2 * 36 + t1 * 6 + t2;
                        let src = q2 * 216 + q1 * 36 + t2 * 6 + t1;
                        table[dst] = self.table[src];
                    }
                }
            }
        }
        Self {
            table,
            initiation_raw: self.initiation_raw,
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

fn matrix_file_stem(matrix: Matrix) -> &'static str {
    match matrix {
        Matrix::T04 => "t04.v4",
        Matrix::T99 => "t99.v2",
    }
}

fn parse_temperatures(raw: Option<&str>) -> Result<[f64; 3]> {
    let mut temps = DEFAULT_TEMPERATURES;
    let Some(raw) = raw else {
        return Ok(temps);
    };
    for (i, part) in raw.split(',').enumerate() {
        if i >= temps.len() {
            bail!("Expected at most 3 temperatures, got '{}'", raw);
        }
        if part.is_empty() {
            bail!("Empty temperature in '{}'", raw);
        }
        temps[i] = part
            .parse::<f64>()
            .with_context(|| format!("Invalid temperature '{}'", part))?;
    }
    for (i, temp) in temps.iter().enumerate() {
        if !temp.is_finite() || !(273.15..=373.15).contains(temp) {
            bail!(
                "Temperatures must be in the range 273.15 - 373.15K; T{} is {}",
                i,
                temp
            );
        }
    }
    if temps[1] >= temps[2] {
        bail!("Temperature interpolation requires T1 < T2");
    }
    Ok(temps)
}

fn load_temperature_table(
    matpath: &Path,
    matrix: &str,
    matrix2: Option<&str>,
    temps: [f64; 3],
) -> Result<(i32, DsmTable)> {
    let exact_path = matrix_path(matpath, temps[0], matrix);
    if exact_path.is_file() {
        return load_dsm_tsv(&exact_path);
    }

    let t1_path = matrix_path(matpath, temps[1], matrix);
    let matrix2 = matrix2.unwrap_or(matrix);
    let t2_path = matrix_path(matpath, temps[2], matrix2);
    let (offset_1, table_1) = load_dsm_tsv(&t1_path)
        .with_context(|| format!("Failed to load T1 matrix {}", t1_path.display()))?;
    let (offset_2, table_2) = load_dsm_tsv(&t2_path)
        .with_context(|| format!("Failed to load T2 matrix {}", t2_path.display()))?;
    Ok(interpolate_tables(
        offset_1, &table_1, offset_2, &table_2, temps,
    ))
}

fn matrix_path(matpath: &Path, temp: f64, matrix: &str) -> PathBuf {
    matpath
        .join("RNA/RNA")
        .join(format!("{temp:.2}"))
        .join(format!("{matrix}.tsv"))
}

fn load_dsm_tsv(path: &Path) -> Result<(i32, DsmTable)> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read DSM TSV {}", path.display()))?;
    let values = text
        .split_whitespace()
        .map(|token| {
            token
                .parse::<f64>()
                .with_context(|| format!("Invalid DSM value '{}' in {}", token, path.display()))
        })
        .collect::<Result<Vec<_>>>()?;

    if values.len() != DSM_TSV_VALUE_COUNT {
        bail!(
            "Read {} DSM values from {}, expected {}",
            values.len(),
            path.display(),
            DSM_TSV_VALUE_COUNT
        );
    }

    let offset = values[0];
    if !offset.is_finite() {
        bail!("Non-finite DSM offset in {}", path.display());
    }
    let initiation_raw = scale_raw(offset);
    let mut table = [[[[0i32; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];
    for (flat, value) in values.iter().skip(1).enumerate() {
        if !value.is_finite() {
            bail!(
                "Non-finite DSM value at position {} in {}",
                flat + 1,
                path.display()
            );
        }
        if value.abs() > MAX_TSV_ENERGY {
            bail!(
                "DSM value {} at position {} exceeds max energy {} in {}",
                value,
                flat + 1,
                MAX_TSV_ENERGY,
                path.display()
            );
        }
        let t2 = TSV_BASE_TO_RUST[flat % BASE_COUNT];
        let t1 = TSV_BASE_TO_RUST[(flat / BASE_COUNT) % BASE_COUNT];
        let q2 = TSV_BASE_TO_RUST[(flat / (BASE_COUNT * BASE_COUNT)) % BASE_COUNT];
        let q1 = TSV_BASE_TO_RUST[(flat / (BASE_COUNT * BASE_COUNT * BASE_COUNT)) % BASE_COUNT];
        table[q1][q2][t1][t2] = scale_raw(-*value);
    }
    Ok((initiation_raw, table))
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
                        interpolate_raw(table_1[q1][q2][t1][t2], table_2[q1][q2][t1][t2], temps);
                }
            }
        }
    }
    (interpolate_raw(offset_1, offset_2, temps), table)
}

fn interpolate_raw(raw_1: i32, raw_2: i32, temps: [f64; 3]) -> i32 {
    let [t0, t1, t2] = temps;
    let diff = i64::from(raw_1) - i64::from(raw_2);
    ((t0 - t2) / (t1 - t2) * diff as f64 + f64::from(raw_2)).round() as i32
}

fn scale_raw(value: f64) -> i32 {
    (value * Energy::RAW_SCALE).round() as i32
}

const T04: DsmTable = [
    [
        [
            [-2000, -123, -123, -123, -123, -123],
            [-123, -40, -40, -40, -40, -40],
            [-123, -40, -40, -40, -40, -40],
            [-123, -40, -40, -40, -40, -40],
            [-123, -40, -40, -40, -40, -40],
            [-123, -40, -40, -40, -40, -40],
        ],
        [
            [-123, 0, 0, 0, 0, 105],
            [-123, -22, -22, -22, -22, -45],
            [-123, -22, -22, -22, -22, -45],
            [-123, -22, -22, -22, -22, -45],
            [-123, -22, -22, -22, -22, -45],
            [-123, -22, -22, -22, -22, -45],
        ],
        [
            [-123, 0, 0, 150, 0, 0],
            [-123, -22, -22, 0, -22, -22],
            [-123, -22, -22, 0, -22, -22],
            [-123, -22, -22, 0, -22, -22],
            [-123, -22, -22, 0, -22, -22],
            [-123, -22, -22, 0, -22, -22],
        ],
        [
            [-123, 0, 150, 0, 0, 105],
            [-123, -22, 0, -22, -22, -45],
            [-123, -22, 0, -22, -22, -45],
            [-123, -22, 0, -22, -22, -45],
            [-123, -22, 0, -22, -22, -45],
            [-123, -22, 0, -22, -22, -45],
        ],
        [
            [-123, 0, 0, 0, 0, 0],
            [-123, -22, -22, -22, -22, -22],
            [-123, -22, -22, -22, -22, -22],
            [-123, -22, -22, -22, -22, -22],
            [-123, -22, -22, -22, -22, -22],
            [-123, -22, -22, -22, -22, -22],
        ],
        [
            [-123, 105, 0, 105, 0, 0],
            [-123, -45, -22, -45, -22, -22],
            [-123, -45, -22, -45, -22, -22],
            [-123, -45, -22, -45, -22, -22],
            [-123, -45, -22, -45, -22, -22],
            [-123, -45, -22, -45, -22, -22],
        ],
    ],
    [
        [
            [-123, -123, -123, -123, -123, -123],
            [-2000, -71, -71, -71, -71, -71],
            [-2000, -71, -71, -71, -71, -71],
            [-2000, -71, -71, -71, -71, -71],
            [-2000, -71, -71, -71, -71, -71],
            [-45, -285, -285, -285, -285, -285],
        ],
        [
            [-40, -22, -22, -22, -22, -45],
            [-71, -22, -22, -22, -22, -70],
            [-71, -22, -22, -22, -22, -70],
            [-71, -22, -22, -22, -22, 30],
            [-71, -22, -22, -22, -22, -70],
            [-285, -230, -230, -150, -230, 90],
        ],
        [
            [-40, -22, -22, 0, -22, -22],
            [-71, -22, -22, 0, -22, -22],
            [-71, -22, -22, 0, -22, -22],
            [-71, -22, -22, 100, -22, -22],
            [-71, -22, -22, 0, -22, -22],
            [-285, -230, -230, 220, -230, -230],
        ],
        [
            [-40, -22, 0, -22, -22, -45],
            [-71, -22, 0, -22, -22, -70],
            [-71, -22, 0, -22, -22, -70],
            [-71, -22, 100, -22, -22, 30],
            [-71, -22, 0, -22, -22, -70],
            [-285, -130, 210, -110, -230, 60],
        ],
        [
            [-40, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
            [-285, -230, -230, -230, -230, -230],
        ],
        [
            [-40, -45, -22, -45, -22, -22],
            [-71, -70, -22, -70, -22, -22],
            [-71, -70, -22, -70, -22, -22],
            [-71, 30, -22, 30, -22, -22],
            [-71, -70, -22, -70, -22, -22],
            [-285, 110, -230, 140, -230, -160],
        ],
    ],
    [
        [
            [-123, -123, -123, -123, -123, -123],
            [-2000, -71, -71, -71, -71, -71],
            [-2000, -71, -71, -71, -71, -71],
            [0, -240, -240, -240, -240, -240],
            [-2000, -71, -71, -71, -71, -71],
            [-2000, -71, -71, -71, -71, -71],
        ],
        [
            [-40, -22, -22, -22, -22, -45],
            [-71, -22, -22, -22, -22, -70],
            [-71, -22, -22, -22, -22, -70],
            [-240, -160, -160, -80, -160, 210],
            [-71, -22, -22, -22, -22, -70],
            [-71, -22, -22, -22, -22, -70],
        ],
        [
            [-40, -22, -22, 0, -22, -22],
            [-71, -22, -22, 0, -22, -22],
            [-71, -22, -22, 0, -22, -22],
            [-240, -160, -160, 330, -160, -160],
            [-71, -22, -22, 0, -22, -22],
            [-71, -22, -22, 0, -22, -22],
        ],
        [
            [-40, -22, 0, -22, -22, -45],
            [-71, -22, 0, -22, -22, -70],
            [-71, -22, 0, -22, -22, -70],
            [-240, -60, 240, -40, -160, 140],
            [-71, -22, 0, -22, -22, -70],
            [-71, -22, 0, -22, -22, -70],
        ],
        [
            [-40, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
            [-240, -160, -160, -160, -160, -160],
            [-71, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
        ],
        [
            [-40, -45, -22, -45, -22, -22],
            [-71, -70, -22, -70, -22, -22],
            [-71, -70, -22, -70, -22, -22],
            [-240, 210, -160, 210, -160, -90],
            [-71, -70, -22, -70, -22, -22],
            [-71, -70, -22, -70, -22, -22],
        ],
    ],
    [
        [
            [-123, -123, -123, -123, -123, -123],
            [-2000, -71, -71, -71, -71, -71],
            [0, -240, -240, -240, -240, -240],
            [-2000, -71, -71, -71, -71, -71],
            [-2000, -71, -71, -71, -71, -71],
            [-45, -285, -285, -285, -285, -285],
        ],
        [
            [-40, -22, -22, -22, -22, -45],
            [-71, -22, -22, -22, -22, 10],
            [-240, -160, -160, -80, -160, 240],
            [-71, -22, -22, -22, -22, 50],
            [-71, -22, -22, -22, -22, -70],
            [-285, -230, -230, -150, -230, 130],
        ],
        [
            [-40, -22, -22, 0, -22, -22],
            [-71, -22, -22, 80, -22, -22],
            [-240, -160, -160, 340, -160, -160],
            [-71, -22, -22, 120, -22, -22],
            [-71, -22, -22, 0, -22, -22],
            [-285, -230, -230, 250, -230, -230],
        ],
        [
            [-40, -22, 0, -22, -22, -45],
            [-71, -22, 80, -22, -22, 10],
            [-240, -60, 330, -40, -160, 150],
            [-71, -22, 120, -22, -22, 50],
            [-71, -22, 0, -22, -22, -70],
            [-285, -130, 210, -110, -230, 50],
        ],
        [
            [-40, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
            [-240, -160, -160, -160, -160, -160],
            [-71, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
            [-285, -230, -230, -230, -230, -230],
        ],
        [
            [-40, -45, -22, -45, -22, -22],
            [-71, 10, -22, 10, -22, -22],
            [-240, 220, -160, 250, -160, -90],
            [-71, 50, -22, 50, -22, -22],
            [-71, -70, -22, -70, -22, -22],
            [-285, 140, -230, -130, -230, -160],
        ],
    ],
    [
        [
            [-123, -123, -123, -123, -123, -123],
            [-2000, -71, -71, -71, -71, -71],
            [-2000, -71, -71, -71, -71, -71],
            [-2000, -71, -71, -71, -71, -71],
            [-2000, -71, -71, -71, -71, -71],
            [-2000, -71, -71, -71, -71, -71],
        ],
        [
            [-40, -22, -22, -22, -22, -45],
            [-71, -22, -22, -22, -22, -70],
            [-71, -22, -22, -22, -22, -70],
            [-71, -22, -22, -22, -22, -70],
            [-71, -22, -22, -22, -22, -70],
            [-71, -22, -22, -22, -22, -70],
        ],
        [
            [-40, -22, -22, 0, -22, -22],
            [-71, -22, -22, 0, -22, -22],
            [-71, -22, -22, 0, -22, -22],
            [-71, -22, -22, 0, -22, -22],
            [-71, -22, -22, 0, -22, -22],
            [-71, -22, -22, 0, -22, -22],
        ],
        [
            [-40, -22, 0, -22, -22, -45],
            [-71, -22, 0, -22, -22, -70],
            [-71, -22, 0, -22, -22, -70],
            [-71, -22, 0, -22, -22, -70],
            [-71, -22, 0, -22, -22, -70],
            [-71, -22, 0, -22, -22, -70],
        ],
        [
            [-40, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
        ],
        [
            [-40, -45, -22, -45, -22, -22],
            [-71, -70, -22, -70, -22, -22],
            [-71, -70, -22, -70, -22, -22],
            [-71, -70, -22, -70, -22, -22],
            [-71, -70, -22, -70, -22, -22],
            [-71, -70, -22, -70, -22, -22],
        ],
    ],
    [
        [
            [-123, -123, -123, -123, -123, -123],
            [-45, -285, -285, -285, -285, -285],
            [-2000, -71, -71, -71, -71, -71],
            [-45, -285, -285, -285, -285, -285],
            [-2000, -71, -71, -71, -71, -71],
            [-2000, -71, -71, -71, -71, -71],
        ],
        [
            [-40, -22, -22, -22, -22, -45],
            [-285, -230, -230, -150, -230, 130],
            [-71, -22, -22, -22, -22, -70],
            [-285, -230, -230, -150, -230, 100],
            [-71, -22, -22, -22, -22, -70],
            [-71, -22, -22, -22, -22, 0],
        ],
        [
            [-40, -22, -22, 0, -22, -22],
            [-285, -230, -230, 240, -230, -230],
            [-71, -22, -22, 0, -22, -22],
            [-285, -230, -230, 150, -230, -230],
            [-71, -22, -22, 0, -22, -22],
            [-71, -22, -22, 70, -22, -22],
        ],
        [
            [-40, -22, 0, -22, -22, -45],
            [-285, -130, 210, -110, -230, 100],
            [-71, -22, 0, -22, -22, -70],
            [-285, -130, 140, -110, -230, -30],
            [-71, -22, 0, -22, -22, -70],
            [-71, -22, 70, -22, -22, 0],
        ],
        [
            [-40, -22, -22, -22, -22, -22],
            [-285, -230, -230, -230, -230, -230],
            [-71, -22, -22, -22, -22, -22],
            [-285, -230, -230, -230, -230, -230],
            [-71, -22, -22, -22, -22, -22],
            [-71, -22, -22, -22, -22, -22],
        ],
        [
            [-40, -45, -22, -45, -22, -22],
            [-285, 90, -230, 130, -230, -160],
            [-71, -70, -22, -70, -22, -22],
            [-285, 60, -230, 50, -230, -160],
            [-71, -70, -22, -70, -22, -22],
            [-71, 0, -22, 0, -22, -22],
        ],
    ],
];

const T99: DsmTable = [
    [
        [
            [-2000, -123, -123, -123, -123, -123],
            [-123, -40, -40, -40, -40, -40],
            [-123, -40, -40, -40, -40, -40],
            [-123, -40, -40, -40, -40, -40],
            [-123, -40, -40, -40, -40, -40],
            [-123, -40, -40, -40, -40, -40],
        ],
        [
            [-123, 0, 0, 0, 0, 105],
            [-123, -24, -24, -24, -24, -45],
            [-123, -24, -24, -24, -24, -45],
            [-123, -24, -24, -24, -24, -45],
            [-123, -24, -24, -24, -24, -45],
            [-123, -24, -24, -24, -24, -45],
        ],
        [
            [-123, 0, 0, 150, 0, 0],
            [-123, -24, -24, 0, -24, -24],
            [-123, -24, -24, 0, -24, -24],
            [-123, -24, -24, 0, -24, -24],
            [-123, -24, -24, 0, -24, -24],
            [-123, -24, -24, 0, -24, -24],
        ],
        [
            [-123, 0, 150, 0, 0, 105],
            [-123, -24, 0, -24, -24, -45],
            [-123, -24, 0, -24, -24, -45],
            [-123, -24, 0, -24, -24, -45],
            [-123, -24, 0, -24, -24, -45],
            [-123, -24, 0, -24, -24, -45],
        ],
        [
            [-123, 0, 0, 0, 0, 0],
            [-123, -24, -24, -24, -24, -24],
            [-123, -24, -24, -24, -24, -24],
            [-123, -24, -24, -24, -24, -24],
            [-123, -24, -24, -24, -24, -24],
            [-123, -24, -24, -24, -24, -24],
        ],
        [
            [-123, 105, 0, 105, 0, 0],
            [-123, -45, -24, -45, -24, -24],
            [-123, -45, -24, -45, -24, -24],
            [-123, -45, -24, -45, -24, -24],
            [-123, -45, -24, -45, -24, -24],
            [-123, -45, -24, -45, -24, -24],
        ],
    ],
    [
        [
            [-123, -123, -123, -123, -123, -123],
            [-2000, -60, -60, -60, -60, -60],
            [-2000, -60, -60, -60, -60, -60],
            [-2000, -60, -60, -60, -60, -60],
            [-2000, -60, -60, -60, -60, -60],
            [-45, -285, -285, -285, -285, -285],
        ],
        [
            [-40, -24, -24, -24, -24, -45],
            [-60, -24, -24, -24, -24, -65],
            [-60, -24, -24, -24, -24, -65],
            [-60, -24, -24, -24, -24, 45],
            [-60, -24, -24, -24, -24, -65],
            [-285, -217, -217, -107, -217, 90],
        ],
        [
            [-40, -24, -24, 0, -24, -24],
            [-60, -24, -24, 0, -24, -24],
            [-60, -24, -24, 0, -24, -24],
            [-60, -24, -24, 110, -24, -24],
            [-60, -24, -24, 0, -24, -24],
            [-285, -217, -217, 220, -217, -217],
        ],
        [
            [-40, -24, 0, -24, -24, -45],
            [-60, -24, 0, -24, -24, -65],
            [-60, -24, 0, -24, -24, -65],
            [-60, -24, 110, -24, -24, 45],
            [-60, -24, 0, -24, -24, -65],
            [-285, -107, 210, -217, -217, 60],
        ],
        [
            [-40, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
            [-285, -217, -217, -217, -217, -217],
        ],
        [
            [-40, -45, -24, -45, -24, -24],
            [-60, -65, -24, -65, -24, -24],
            [-60, -65, -24, -65, -24, -24],
            [-60, 45, -24, 45, -24, -24],
            [-60, -65, -24, -65, -24, -24],
            [-285, 110, -217, 140, -217, -147],
        ],
    ],
    [
        [
            [-123, -123, -123, -123, -123, -123],
            [-2000, -60, -60, -60, -60, -60],
            [-2000, -60, -60, -60, -60, -60],
            [0, -240, -240, -240, -240, -240],
            [-2000, -60, -60, -60, -60, -60],
            [-2000, -60, -60, -60, -60, -60],
        ],
        [
            [-40, -24, -24, -24, -24, -45],
            [-60, -24, -24, -24, -24, -65],
            [-60, -24, -24, -24, -24, -65],
            [-240, -152, -152, -42, -152, 210],
            [-60, -24, -24, -24, -24, -65],
            [-60, -24, -24, -24, -24, -65],
        ],
        [
            [-40, -24, -24, 0, -24, -24],
            [-60, -24, -24, 0, -24, -24],
            [-60, -24, -24, 0, -24, -24],
            [-240, -152, -152, 330, -152, -152],
            [-60, -24, -24, 0, -24, -24],
            [-60, -24, -24, 0, -24, -24],
        ],
        [
            [-40, -24, 0, -24, -24, -45],
            [-60, -24, 0, -24, -24, -65],
            [-60, -24, 0, -24, -24, -65],
            [-240, -42, 240, -152, -152, 140],
            [-60, -24, 0, -24, -24, -65],
            [-60, -24, 0, -24, -24, -65],
        ],
        [
            [-40, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
            [-240, -152, -152, -152, -152, -152],
            [-60, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
        ],
        [
            [-40, -45, -24, -45, -24, -24],
            [-60, -65, -24, -65, -24, -24],
            [-60, -65, -24, -65, -24, -24],
            [-240, 210, -152, 210, -152, -82],
            [-60, -65, -24, -65, -24, -24],
            [-60, -65, -24, -65, -24, -24],
        ],
    ],
    [
        [
            [-123, -123, -123, -123, -123, -123],
            [-2000, -60, -60, -60, -60, -60],
            [0, -240, -240, -240, -240, -240],
            [-2000, -60, -60, -60, -60, -60],
            [-2000, -60, -60, -60, -60, -60],
            [-45, -285, -285, -285, -285, -285],
        ],
        [
            [-40, -24, -24, -24, -24, -45],
            [-60, -24, -24, -24, -24, 45],
            [-240, -152, -152, -42, -152, 240],
            [-60, -24, -24, -24, -24, -65],
            [-60, -24, -24, -24, -24, -65],
            [-285, -217, -217, -107, -217, 130],
        ],
        [
            [-40, -24, -24, 0, -24, -24],
            [-60, -24, -24, 110, -24, -24],
            [-240, -152, -152, 340, -152, -152],
            [-60, -24, -24, 0, -24, -24],
            [-60, -24, -24, 0, -24, -24],
            [-285, -217, -217, 250, -217, -217],
        ],
        [
            [-40, -24, 0, -24, -24, -45],
            [-60, -24, 110, -24, -24, 45],
            [-240, -42, 330, -152, -152, 150],
            [-60, -24, 0, -24, -24, -65],
            [-60, -24, 0, -24, -24, -65],
            [-285, -107, 210, -217, -217, 50],
        ],
        [
            [-40, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
            [-240, -152, -152, -152, -152, -152],
            [-60, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
            [-285, -217, -217, -217, -217, -217],
        ],
        [
            [-40, -45, -24, -45, -24, -24],
            [-60, 45, -24, 45, -24, -24],
            [-240, 220, -152, 250, -152, -82],
            [-60, -65, -24, -65, -24, -24],
            [-60, -65, -24, -65, -24, -24],
            [-285, 140, -217, -130, -217, -147],
        ],
    ],
    [
        [
            [-123, -123, -123, -123, -123, -123],
            [-2000, -60, -60, -60, -60, -60],
            [-2000, -60, -60, -60, -60, -60],
            [-2000, -60, -60, -60, -60, -60],
            [-2000, -60, -60, -60, -60, -60],
            [-2000, -60, -60, -60, -60, -60],
        ],
        [
            [-40, -24, -24, -24, -24, -45],
            [-60, -24, -24, -24, -24, -65],
            [-60, -24, -24, -24, -24, -65],
            [-60, -24, -24, -24, -24, -65],
            [-60, -24, -24, -24, -24, -65],
            [-60, -24, -24, -24, -24, -65],
        ],
        [
            [-40, -24, -24, 0, -24, -24],
            [-60, -24, -24, 0, -24, -24],
            [-60, -24, -24, 0, -24, -24],
            [-60, -24, -24, 0, -24, -24],
            [-60, -24, -24, 0, -24, -24],
            [-60, -24, -24, 0, -24, -24],
        ],
        [
            [-40, -24, 0, -24, -24, -45],
            [-60, -24, 0, -24, -24, -65],
            [-60, -24, 0, -24, -24, -65],
            [-60, -24, 0, -24, -24, -65],
            [-60, -24, 0, -24, -24, -65],
            [-60, -24, 0, -24, -24, -65],
        ],
        [
            [-40, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
        ],
        [
            [-40, -45, -24, -45, -24, -24],
            [-60, -65, -24, -65, -24, -24],
            [-60, -65, -24, -65, -24, -24],
            [-60, -65, -24, -65, -24, -24],
            [-60, -65, -24, -65, -24, -24],
            [-60, -65, -24, -65, -24, -24],
        ],
    ],
    [
        [
            [-123, -123, -123, -123, -123, -123],
            [-45, -285, -285, -285, -285, -285],
            [-2000, -60, -60, -60, -60, -60],
            [-45, -285, -285, -285, -285, -285],
            [-2000, -60, -60, -60, -60, -60],
            [-2000, -60, -60, -60, -60, -60],
        ],
        [
            [-40, -24, -24, -24, -24, -45],
            [-285, -217, -217, -107, -217, 130],
            [-60, -24, -24, -24, -24, -65],
            [-285, -217, -217, -107, -217, 100],
            [-60, -24, -24, -24, -24, -65],
            [-60, -24, -24, -24, -24, 5],
        ],
        [
            [-40, -24, -24, 0, -24, -24],
            [-285, -217, -217, 240, -217, -217],
            [-60, -24, -24, 0, -24, -24],
            [-285, -217, -217, 150, -217, -217],
            [-60, -24, -24, 0, -24, -24],
            [-60, -24, -24, 70, -24, -24],
        ],
        [
            [-40, -24, 0, -24, -24, -45],
            [-285, -107, 210, -217, -217, 100],
            [-60, -24, 0, -24, -24, -65],
            [-285, -107, 140, -217, -217, -30],
            [-60, -24, 0, -24, -24, -65],
            [-60, -24, 70, -24, -24, 5],
        ],
        [
            [-40, -24, -24, -24, -24, -24],
            [-285, -217, -217, -217, -217, -217],
            [-60, -24, -24, -24, -24, -24],
            [-285, -217, -217, -217, -217, -217],
            [-60, -24, -24, -24, -24, -24],
            [-60, -24, -24, -24, -24, -24],
        ],
        [
            [-40, -45, -24, -45, -24, -24],
            [-285, 90, -217, 130, -217, -147],
            [-60, -65, -24, -65, -24, -24],
            [-285, 60, -217, 50, -217, -147],
            [-60, -65, -24, -65, -24, -24],
            [-60, 5, -24, 5, -24, -24],
        ],
    ],
];

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn pairing_distinguishes_seed_and_extension_modes() {
        let model = ScoringModel::new(Matrix::T04, 0);

        // DP/display pairing follows the scoring matrix semantics.
        assert!(model.is_pair(Base::G, Base::G));
        assert!(model.is_pair(Base::G, Base::A));

        // Seed pairing can be stricter than the extension model.
        assert!(!ScoringModel::seed_pair(Base::G, Base::A, false));
        assert!(ScoringModel::seed_pair(Base::G, Base::A, true));
    }

    #[test]
    fn transition_energy_returns_nonzero_for_valid_pairs() {
        let model = ScoringModel::new(Matrix::T04, 0);
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
    fn gg_cc_stack_is_strongest_at_330() {
        let model = ScoringModel::new(Matrix::T04, 0);
        // GG/CC stack is the strongest at 330 (3.30 kcal/mol)
        // Original target was CC, index space is GG
        let actual = model.transition_energy_bases(Base::G, Base::G, Base::G, Base::G);
        assert_eq!(actual, 33_000);
    }

    #[test]
    fn built_in_energy_conversion_preserves_previous_kcal_output() {
        let model = ScoringModel::new(Matrix::T04, 0);
        let energy = model.energy_from_raw(33_000_i64);
        assert_eq!(f64::from(energy), 2.29);
    }

    #[test]
    fn transpose_swaps_both_pairs() {
        let right = ScoringModel::new(Matrix::T04, 50);
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
    fn load_dsm_tsv_parses_offset_and_flat_order() {
        let file = write_temp_tsv(&tsv_with_overrides(
            1.2345,
            0.0,
            &[(0, 0.1), (1, 0.2), (DSM_FLAT_SIZE - 1, 3.0)],
        ));
        let (offset, table) = load_dsm_tsv(file.path()).unwrap();
        assert_eq!(offset, 12_345);
        assert_eq!(
            table[Base::A.idx()][Base::A.idx()][Base::A.idx()][Base::A.idx()],
            -1_000
        );
        assert_eq!(
            table[Base::A.idx()][Base::A.idx()][Base::A.idx()][Base::C.idx()],
            -2_000
        );
        assert_eq!(
            table[Base::Gap.idx()][Base::Gap.idx()][Base::Gap.idx()][Base::Gap.idx()],
            -30_000
        );
    }

    #[test]
    fn load_dsm_tsv_rejects_wrong_count_and_invalid_values() {
        let too_short = write_temp_tsv("0.0\n1.0\n");
        assert!(load_dsm_tsv(too_short.path()).is_err());

        let nan = write_temp_tsv(&tsv_with_overrides(0.0, 0.0, &[(0, f64::NAN)]));
        assert!(load_dsm_tsv(nan.path()).is_err());

        let too_large = write_temp_tsv(&tsv_with_overrides(0.0, 0.0, &[(0, 20.0001)]));
        assert!(load_dsm_tsv(too_large.path()).is_err());
    }

    #[test]
    fn score_config_loads_exact_temperature_table_when_present() {
        let dir = tempfile::tempdir().unwrap();
        write_matrix(dir.path(), 310.15, "t04.v4", &uniform_tsv(1.0, 1.0));
        write_matrix(dir.path(), 315.15, "t04.v4", &uniform_tsv(3.0, 3.0));

        let score = test_score_config(dir.path(), Some("310.15"));
        let model = ScoringModel::from_score_config(&score).unwrap();
        assert_eq!(model.initiation_raw, 10_000);
        assert_eq!(
            model.transition_energy_bases(Base::A, Base::A, Base::U, Base::U),
            -10_000
        );
    }

    #[test]
    fn score_config_interpolates_missing_temperature_table() {
        let dir = tempfile::tempdir().unwrap();
        write_matrix(dir.path(), 310.15, "t04.v4", &uniform_tsv(1.0, 1.0));
        write_matrix(dir.path(), 315.15, "t04.v4", &uniform_tsv(3.0, 3.0));

        let score = test_score_config(dir.path(), Some("312.65"));
        let model = ScoringModel::from_score_config(&score).unwrap();
        assert_eq!(model.initiation_raw, 20_000);
        assert_eq!(
            model.transition_energy_bases(Base::A, Base::A, Base::U, Base::U),
            -20_000
        );
    }

    #[test]
    fn parse_temperatures_rejects_bad_input() {
        assert!(parse_temperatures(Some("310.15,315.15,320.15,325.15")).is_err());
        assert!(parse_temperatures(Some("310.15,,315.15")).is_err());
        assert!(parse_temperatures(Some("272.0")).is_err());
        assert!(parse_temperatures(Some("315.15,315.15,310.15")).is_err());
    }

    fn test_score_config(path: &std::path::Path, temperature: Option<&str>) -> ScoreConfig {
        ScoreConfig {
            matrix: Matrix::T04,
            penalty: Energy::from(0.0),
            matrix2: None,
            matpath: Some(path.display().to_string()),
            temperature: temperature.map(str::to_owned),
            weights: None,
        }
    }

    fn write_matrix(root: &std::path::Path, temp: f64, name: &str, content: &str) {
        let dir = root.join("RNA/RNA").join(format!("{temp:.2}"));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join(format!("{name}.tsv")), content).unwrap();
    }

    fn write_temp_tsv(content: &str) -> tempfile::NamedTempFile {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file
    }

    fn uniform_tsv(offset: f64, value: f64) -> String {
        tsv_with_overrides(offset, value, &[])
    }

    fn tsv_with_overrides(offset: f64, value: f64, overrides: &[(usize, f64)]) -> String {
        let mut values = vec![value; DSM_FLAT_SIZE];
        for &(idx, override_value) in overrides {
            values[idx] = override_value;
        }
        let mut out = format!("{offset:.4}");
        for value in values {
            out.push('\t');
            if value.is_nan() {
                out.push_str("NaN");
            } else {
                out.push_str(&format!("{value:.4}"));
            }
        }
        out
    }
}
