use crate::seq::SeqView;
use crate::types::Base;

/// Normalize a single base for RNA output (T→U, uppercase).
#[inline]
fn normalize_base_for_rna(b: u8) -> u8 {
    match b {
        b'T' | b't' => b'U',
        other => other.to_ascii_uppercase(),
    }
}

/// Normalize a Base for RNA output (uppercase, U not T).
#[inline]
fn normalize_base_enum_for_rna(b: Base) -> u8 {
    b.to_u8_upper()
}

/// Convert bytes to an RNA string (T→U, uppercase) and optionally reverse.
pub fn bytes_to_rna_string(s: &[u8], reverse: bool) -> String {
    let mut result = String::with_capacity(s.len());
    if reverse {
        for &b in s.iter().rev() {
            result.push(normalize_base_for_rna(b) as char);
        }
    } else {
        for &b in s {
            result.push(normalize_base_for_rna(b) as char);
        }
    }
    result
}

/// Convert Base slice to an RNA string (uppercase) and optionally reverse.
pub fn bases_to_rna_string(s: SeqView<'_>, reverse: bool) -> String {
    let s = s.as_slice();
    let mut result = String::with_capacity(s.len());
    if reverse {
        for &b in s.iter().rev() {
            result.push(normalize_base_enum_for_rna(b) as char);
        }
    } else {
        for &b in s {
            result.push(normalize_base_enum_for_rna(b) as char);
        }
    }
    result
}

/// Append RNA-formatted bytes (uppercased, T→U) to an existing buffer.
pub fn push_bytes_as_rna(buf: &mut Vec<u8>, s: &[u8], reverse: bool) {
    if reverse {
        for &b in s.iter().rev() {
            buf.push(normalize_base_for_rna(b));
        }
    } else {
        for &b in s {
            buf.push(normalize_base_for_rna(b));
        }
    }
}

/// Append RNA-formatted bases (uppercased) to an existing buffer.
pub fn push_bases_as_rna(buf: &mut Vec<u8>, s: SeqView<'_>, reverse: bool) {
    let s = s.as_slice();
    if reverse {
        for &b in s.iter().rev() {
            buf.push(normalize_base_enum_for_rna(b));
        }
    } else {
        for &b in s {
            buf.push(normalize_base_enum_for_rna(b));
        }
    }
}
