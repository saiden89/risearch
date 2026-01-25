/// Normalize a single base for RNA output (T→U, uppercase).
#[inline]
fn normalize_base_for_rna(b: u8) -> u8 {
    match b {
        b'T' | b't' => b'U',
        other => other.to_ascii_uppercase(),
    }
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
