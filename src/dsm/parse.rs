use anyhow::{bail, Context, Result};

use crate::types::{Base, Energy};

use super::{
    DsmTable, Orientation, BASE_COUNT, DSM_FLAT_SIZE, DSM_HEADER,
};

pub(crate) fn load_dsm_tsv_text(
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
    if header_cols != DSM_HEADER {
        bail!(
            "Invalid DSM TSV header {:?}, expected {:?}",
            header_cols,
            DSM_HEADER
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

        let mut q1 = Base::try_from(cols[0].chars().next().unwrap_or('?'))
            .map_err(|e| anyhow::anyhow!(e))
            .with_context(|| format!("Invalid q1 at DSM TSV row {}", line_no + 2))?
            .idx();
        let mut q2 = Base::try_from(cols[1].chars().next().unwrap_or('?'))
            .map_err(|e| anyhow::anyhow!(e))
            .with_context(|| format!("Invalid q2 at DSM TSV row {}", line_no + 2))?
            .idx();
        let mut t1 = Base::try_from(cols[2].chars().next().unwrap_or('?'))
            .map_err(|e| anyhow::anyhow!(e))
            .with_context(|| format!("Invalid t1 at DSM TSV row {}", line_no + 2))?
            .idx();
        let mut t2 = Base::try_from(cols[3].chars().next().unwrap_or('?'))
            .map_err(|e| anyhow::anyhow!(e))
            .with_context(|| format!("Invalid t2 at DSM TSV row {}", line_no + 2))?
            .idx();
        if reverse_swap {
            (q1, q2, t1, t2) = (t2, t1, q2, q1);
        }

        let idx = q1 * 216 + q2 * 36 + t1 * 6 + t2;
        if std::mem::replace(&mut seen[idx], true) {
            bail!("Duplicate DSM coordinate at row {}", line_no + 2);
        }

        table[q1][q2][t1][t2] = Energy::try_from_kcal(-cols[4].parse::<f64>().with_context(|| {
            format!("Invalid delta_g '{}' at row {}", cols[4], line_no + 2)
        })?)
        .map_err(|e| anyhow::anyhow!(e))
        .with_context(|| format!("Invalid energy at row {}", line_no + 2))?
        .0;
        row_count += 1;
    }

    if row_count == 0 {
        bail!("DSM TSV has no transition rows");
    }

    Ok((
        Energy::try_from_kcal(initiation)
            .map_err(|e| anyhow::anyhow!(e))
            .context("Invalid initiation energy")?,
        table,
    ))
}
