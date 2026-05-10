use anyhow::{bail, Context, Result};

use crate::types::{Base, Energy};

use super::{
    DsmTable, Orientation, BASE_COUNT, CANONICAL_DSM_HEADER, DSM_FLAT_SIZE, MAX_TSV_ENERGY,
};

pub(super) fn load_canonical_dsm_tsv_text(
    text: &str,
    initiation: f64,
    invalid_transition: f64,
    orientation: Orientation,
) -> Result<(Energy, DsmTable)> {
    if !invalid_transition.is_finite() {
        bail!("Non-finite invalid transition energy");
    }

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

    let invalid_score = Energy::from_kcal(-invalid_transition).0;
    let mut table = [[[[invalid_score; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];
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
        table[q1][q2][t1][t2] = Energy::from_kcal(-delta_g).0;
        row_count += 1;
    }

    if row_count == 0 {
        bail!("DSM TSV has no transition rows");
    }

    Ok((Energy::from_kcal(initiation), table))
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
