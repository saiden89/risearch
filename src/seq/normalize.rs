use anyhow::{Result, bail};
use std::sync::OnceLock;

#[derive(Default, Debug, Clone, Copy)]
pub(crate) struct NormalizationStats {
    pub(crate) removed_gaps: usize,
    pub(crate) converted_to_n: usize,
}

#[derive(Copy, Clone)]
enum NormalizeKind {
    /// Keep the normalized base as is (lowercase DNA or `n`).
    Copy,
    /// Map `u`/`t` to `t`.
    MapToT,
    /// Map ambiguous alphabetic bases to `n`.
    MapToN,
    /// Skip gap characters.
    SkipGap,
    /// Invalid/unsupported byte.
    Error,
}

#[derive(Copy, Clone)]
struct NormalizeEntry {
    norm: u8,
    kind: NormalizeKind,
}

fn classify_byte(c: u8) -> NormalizeEntry {
    let lower = c.to_ascii_lowercase();
    match lower {
        b'a' | b'c' | b'g' | b'n' => NormalizeEntry {
            norm: lower,
            kind: NormalizeKind::Copy,
        },
        b'u' | b't' => NormalizeEntry {
            norm: b't',
            kind: NormalizeKind::MapToT,
        },
        b'-' | b'.' => NormalizeEntry {
            norm: 0,
            kind: NormalizeKind::SkipGap,
        },
        _ if lower.is_ascii_alphabetic() => NormalizeEntry {
            norm: b'n',
            kind: NormalizeKind::MapToN,
        },
        _ => NormalizeEntry {
            norm: 0,
            kind: NormalizeKind::Error,
        },
    }
}
fn lookup_table() -> &'static [NormalizeEntry; 256] {
    static LOOKUP: OnceLock<[NormalizeEntry; 256]> = OnceLock::new();
    LOOKUP.get_or_init(|| {
        let mut table = [NormalizeEntry {
            norm: 0,
            kind: NormalizeKind::Error,
        }; 256];
        let mut i = 0;
        while i < 256 {
            table[i] = classify_byte(i as u8);
            i += 1;
        }
        table
    })
}

pub(crate) fn normalize_rna_sequence(
    id: &str,
    seq: &[u8],
) -> Result<(Vec<u8>, NormalizationStats)> {
    let mut out = Vec::with_capacity(seq.len());
    let mut stats = NormalizationStats::default();
    let table = lookup_table();

    for &b in seq {
        let entry = table[b as usize];
        match entry.kind {
            NormalizeKind::Copy | NormalizeKind::MapToT => out.push(entry.norm),
            NormalizeKind::MapToN => {
                out.push(entry.norm);
                stats.converted_to_n += 1;
            }
            NormalizeKind::SkipGap => {
                stats.removed_gaps += 1;
            }
            NormalizeKind::Error => {
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
