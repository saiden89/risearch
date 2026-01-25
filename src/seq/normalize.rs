use anyhow::{Result, bail};

#[derive(Default, Debug, Clone, Copy)]
pub(crate) struct NormalizationStats {
    pub(crate) removed_gaps: usize,
    pub(crate) converted_to_n: usize,
}

pub(crate) fn normalize_rna_sequence(
    id: &str,
    seq: &[u8],
) -> Result<(Vec<u8>, NormalizationStats)> {
    let mut out = Vec::with_capacity(seq.len());
    let mut stats = NormalizationStats::default();

    for &b in seq {
        let c = b.to_ascii_lowercase();
        match c {
            b'a' | b'c' | b'g' => out.push(c),
            b'u' | b't' => out.push(b't'),
            b'n' => out.push(b'n'),
            b'-' | b'.' => {
                stats.removed_gaps += 1;
            }
            _ if c.is_ascii_alphabetic() => {
                // Preserve behavior similar to legacy RIsearch2: map ambiguous bases to N.
                out.push(b'n');
                stats.converted_to_n += 1;
            }
            _ => {
                bail!(
                    "Invalid character in sequence '{}': byte=0x{:02X} ('{}')",
                    id,
                    b,
                    b as char
                );
            }
        }
    }

    Ok((out, stats))
}
