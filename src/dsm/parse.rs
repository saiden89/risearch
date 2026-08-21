use crate::error::{Error, Result};
use crate::types::{Base, Energy};

use super::{flat_idx, DsmTable, Orientation, BASE_COUNT, DSM_FLAT_SIZE, DSM_HEADER};

pub(crate) fn load_dsm_tsv_text(
    text: &str,
    initiation: f64,
    invalid_transition: f64,
    orientation: Orientation,
) -> Result<(Energy, DsmTable)> {
    if !invalid_transition.is_finite() {
        return Err(Error::Dsm("Non-finite invalid transition energy".into()));
    }

    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| Error::Dsm("DSM TSV is empty".into()))?;
    let header_cols = header.split_whitespace().collect::<Vec<_>>();
    if header_cols != DSM_HEADER {
        return Err(Error::Dsm(format!(
            "Invalid DSM TSV header {header_cols:?}, expected {DSM_HEADER:?}"
        )));
    }

    let reverse_swap = matches!(orientation, Orientation::ReverseSwap);

    let invalid_score = Energy::try_from(-invalid_transition)
        .map_err(|e| Error::Dsm(format!("Invalid transition energy: {e}")))?
        .0;
    let mut table = [[[[invalid_score; BASE_COUNT]; BASE_COUNT]; BASE_COUNT]; BASE_COUNT];
    let mut seen = [false; DSM_FLAT_SIZE];
    let mut row_count = 0usize;
    for (line_no, line) in lines.enumerate() {
        let row = line_no + 2;
        let cols = line.split_whitespace().collect::<Vec<_>>();
        if cols.len() != 5 {
            return Err(Error::Dsm(format!(
                "Invalid DSM TSV row {row}: expected 5 columns, got {}",
                cols.len()
            )));
        }

        let base = |col: usize, field: &str| -> Result<usize> {
            Base::try_from(cols[col].chars().next().unwrap_or('?'))
                .map(Base::as_usize)
                .map_err(|e| Error::Dsm(format!("Invalid {field} at DSM TSV row {row}: {e}")))
        };
        let mut q1 = base(0, "q1")?;
        let mut q2 = base(1, "q2")?;
        let mut t1 = base(2, "t1")?;
        let mut t2 = base(3, "t2")?;
        if reverse_swap {
            (q1, q2, t1, t2) = (t2, t1, q2, q1);
        }

        let idx = flat_idx(q1 as u8, q2 as u8, t1 as u8, t2 as u8);
        if std::mem::replace(&mut seen[idx], true) {
            return Err(Error::Dsm(format!("Duplicate DSM coordinate at row {row}")));
        }

        let delta_g = cols[4]
            .parse::<f64>()
            .map_err(|e| Error::Dsm(format!("Invalid delta_g '{}' at row {row}: {e}", cols[4])))?;
        table[q1][q2][t1][t2] = Energy::try_from(-delta_g)
            .map_err(|e| Error::Dsm(format!("Invalid energy at row {row}: {e}")))?
            .0;
        row_count += 1;
    }

    if row_count == 0 {
        return Err(Error::Dsm("DSM TSV has no transition rows".into()));
    }

    Ok((
        Energy::try_from(initiation)
            .map_err(|e| Error::Dsm(format!("Invalid initiation energy: {e}")))?,
        table,
    ))
}
