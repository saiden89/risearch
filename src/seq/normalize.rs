use crate::error::{Error, Result};
use std::sync::OnceLock;

use crate::types::Base;

/// What [`Sequence::normalize`](crate::Sequence::normalize) changed while
/// reading a raw record.
#[derive(Default, Debug, Clone, Copy)]
pub struct NormalizationStats {
    /// Gap characters dropped from the input.
    pub removed_gaps: usize,
    /// Unrecognised characters rewritten as `N`.
    pub converted_to_n: usize,
}

#[derive(Copy, Clone)]
enum NormalizeKind {
    /// Keep the normalized base as is (lowercase DNA or `n`).
    Copy,
    /// Map `u`/`t` to U.
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
    norm: Base,
    kind: NormalizeKind,
}

fn classify_byte(c: u8) -> NormalizeEntry {
    let lower = c.to_ascii_lowercase();
    let upper = lower.to_ascii_uppercase();
    match upper {
        b'A' => NormalizeEntry {
            norm: Base::A,
            kind: NormalizeKind::Copy,
        },
        b'C' => NormalizeEntry {
            norm: Base::C,
            kind: NormalizeKind::Copy,
        },
        b'G' => NormalizeEntry {
            norm: Base::G,
            kind: NormalizeKind::Copy,
        },
        b'N' => NormalizeEntry {
            norm: Base::N,
            kind: NormalizeKind::Copy,
        },
        b'U' | b'T' => NormalizeEntry {
            norm: Base::U,
            kind: NormalizeKind::MapToT,
        },
        b'-' | b'.' => NormalizeEntry {
            norm: Base::Gap,
            kind: NormalizeKind::SkipGap,
        },
        _ if upper.is_ascii_alphabetic() => NormalizeEntry {
            norm: Base::N,
            kind: NormalizeKind::MapToN,
        },
        _ => NormalizeEntry {
            norm: Base::Gap,
            kind: NormalizeKind::Error,
        },
    }
}
fn lookup_table() -> &'static [NormalizeEntry; 256] {
    static LOOKUP: OnceLock<[NormalizeEntry; 256]> = OnceLock::new();
    LOOKUP.get_or_init(|| {
        let mut table = [NormalizeEntry {
            norm: Base::Gap,
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
) -> Result<(Vec<Base>, NormalizationStats)> {
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
                return Err(Error::Input(format!(
                    "Invalid character in sequence '{}': byte=0x{:02X} ('{}')",
                    id, b, b as char
                )));
            }
        }
    }

    Ok((out, stats))
}
