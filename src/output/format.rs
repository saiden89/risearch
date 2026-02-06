use std::io::Write;

use crate::alignment::Alignment;
use crate::config::OutputFormat;
use crate::registry::{QueryRegistry, TargetRegistry};
use crate::search::SearchHit;
use crate::seq::utils::push_bases_as_rna;
use crate::seq::Sequence;

/// Reusable buffers for hit formatting (avoids per-hit allocation).
pub struct OutputBuffers {
    line: Vec<u8>,
    itoa: itoa::Buffer,
    zmij: zmij::Buffer,
}

impl OutputBuffers {
    pub fn new() -> Self {
        Self {
            line: Vec::with_capacity(256), // Preallocate typical line size
            itoa: itoa::Buffer::new(),
            zmij: zmij::Buffer::new(),
        }
    }

    /// Clear buffers for reuse (keeps capacity allocated)
    fn clear(&mut self) {
        self.line.clear();
        // itoa and zmij buffers are stack-allocated, no need to clear
    }
}

impl Default for OutputBuffers {
    fn default() -> Self {
        Self::new()
    }
}

#[inline]
fn push_usize(buf: &mut Vec<u8>, itoa_buf: &mut itoa::Buffer, val: usize) {
    buf.extend_from_slice(itoa_buf.format(val).as_bytes());
}

#[inline]
fn push_score_fixed_2(
    buf: &mut Vec<u8>,
    itoa_buf: &mut itoa::Buffer,
    zmij_buf: &mut zmij::Buffer,
    score: f64,
) {
    // Round to 2 decimals with ties-to-even to match std formatting.
    let rounded = (score * 100.0).round_ties_even() / 100.0;
    let start = buf.len();
    let s = zmij_buf.format_finite(rounded);
    let bytes = s.as_bytes();

    // Avoid scientific notation in output; fall back to integer-based formatting.
    if bytes.iter().any(|&b| b == b'e' || b == b'E') {
        buf.truncate(start);
        push_score_fixed_2_int(buf, itoa_buf, score);
        return;
    }

    buf.extend_from_slice(bytes);

    if let Some(dot) = bytes.iter().position(|&b| b == b'.') {
        let decimals = bytes.len() - dot - 1;
        match decimals {
            0 => buf.extend_from_slice(b"00"),
            1 => buf.push(b'0'),
            2 => {}
            _ => {
                // Rare: if shortest representation has more decimals, use fixed formatter.
                buf.truncate(start);
                push_score_fixed_2_int(buf, itoa_buf, score);
            }
        }
    } else {
        buf.extend_from_slice(b".00");
    }
}

#[inline]
fn push_score_fixed_2_int(buf: &mut Vec<u8>, itoa_buf: &mut itoa::Buffer, score: f64) {
    let scaled = (score * 100.0).round_ties_even() as i64;
    let mut v = scaled;
    if v < 0 {
        buf.push(b'-');
        v = -v;
    }
    let int_part = (v / 100) as u64;
    let frac = (v % 100) as u8;
    buf.extend_from_slice(itoa_buf.format(int_part).as_bytes());
    buf.push(b'.');
    buf.push(b'0' + (frac / 10));
    buf.push(b'0' + (frac % 10));
}

#[allow(clippy::too_many_arguments)]
#[inline]
pub(super) fn push_result_fields(
    buf: &mut Vec<u8>,
    itoa_buf: &mut itoa::Buffer,
    zmij_buf: &mut zmij::Buffer,
    q_id: &str,
    q_start: usize,
    q_end: usize,
    t_id: &str,
    t_start: usize,
    t_end: usize,
    strand_char: char,
    score: f64,
) {
    buf.extend_from_slice(q_id.as_bytes());
    buf.push(b'\t');
    push_usize(buf, itoa_buf, q_start);
    buf.push(b'\t');
    push_usize(buf, itoa_buf, q_end);
    buf.push(b'\t');
    buf.extend_from_slice(t_id.as_bytes());
    buf.push(b'\t');
    push_usize(buf, itoa_buf, t_start);
    buf.push(b'\t');
    push_usize(buf, itoa_buf, t_end);
    buf.push(b'\t');
    buf.push(strand_char as u8);
    buf.push(b'\t');
    push_score_fixed_2(buf, itoa_buf, zmij_buf, score);
}

#[inline]
fn truncate_id(id: &str, max_len: Option<usize>) -> &str {
    let base = id.split_whitespace().next().unwrap_or(id);
    if let Some(max) = max_len {
        if base.len() > max {
            &base[..max]
        } else {
            base
        }
    } else {
        base
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn fill_line_buf(
    line_buf: &mut Vec<u8>,
    itoa_buf: &mut itoa::Buffer,
    zmij_buf: &mut zmij::Buffer,
    format: OutputFormat,
    q_id: &str,
    q_start: usize,
    q_end: usize,
    t_id: &str,
    t_start: usize,
    t_end: usize,
    strand_char: char,
    score: f64,
    alignment: Option<&Alignment>,
    flank_5: (&Sequence, std::ops::Range<usize>, bool),
    flank_3: (&Sequence, std::ops::Range<usize>, bool),
    id_max_len: Option<usize>,
) {
    let q_id_trunc = truncate_id(q_id, id_max_len);
    let t_id_trunc = truncate_id(t_id, id_max_len);
    let steps_len = alignment.map(|a| a.steps().len()).unwrap_or(0);

    let approx = q_id_trunc.len()
        + t_id_trunc.len()
        + (steps_len * 2)
        + flank_5.1.end.saturating_sub(flank_5.1.start)
        + flank_3.1.end.saturating_sub(flank_3.1.start)
        + 96;
    line_buf.clear();
    if line_buf.capacity() < approx {
        line_buf.reserve(approx - line_buf.capacity());
    }

    if format == OutputFormat::Detailed {
        if let Some(align) = alignment {
            align.write_query_seq(line_buf);
            line_buf.push(b'\n');
            align.write_alignment_line(line_buf);
            line_buf.push(b'\n');
            align.write_target_seq(line_buf);
            line_buf.push(b'\n');
        }
    }

    push_result_fields(
        line_buf,
        itoa_buf,
        zmij_buf,
        q_id_trunc,
        q_start,
        q_end,
        t_id_trunc,
        t_start,
        t_end,
        strand_char,
        score,
    );

    match format {
        OutputFormat::Minimal | OutputFormat::Detailed => {
            line_buf.push(b'\n');
        }
        OutputFormat::Cigar => {
            line_buf.push(b'\t');
            if let Some(align) = alignment {
                align.write_pairing_string(line_buf);
            }
            line_buf.push(b'\n');
        }
        OutputFormat::BindingSite => {
            line_buf.push(b'\t');
            if let Some(align) = alignment {
                align.write_pairing_string(line_buf);
                line_buf.push(b'\t');
                align.write_target_seq(line_buf);
            } else {
                line_buf.push(b'\t'); // empty pairing
                                      // empty target seq
            }
            line_buf.push(b'\t');
            push_bases_as_rna(
                line_buf,
                &flank_5.0[flank_5.1.start..flank_5.1.end],
                flank_5.2,
            );
            line_buf.push(b'\t');
            push_bases_as_rna(
                line_buf,
                &flank_3.0[flank_3.1.start..flank_3.1.end],
                flank_3.2,
            );
            line_buf.push(b'\n');
        }
    }
}

/// Write a hit with format using provided buffers (for efficient reuse).
///
/// Buffers should be cleared between calls using `bufs.clear()`.
pub fn write_hit_with_format<W: Write + ?Sized>(
    bufs: &mut OutputBuffers,
    hit: &SearchHit,
    format: OutputFormat,
    writer: &mut W,
    query_registry: &QueryRegistry,
    target_registry: &TargetRegistry,
) -> std::io::Result<()> {
    bufs.clear(); // Clear from previous use
    let q_name = query_registry.get_name(hit.query_idx);
    let t_name = target_registry.get_name(hit.target_idx);
    fill_line_buf(
        &mut bufs.line,
        &mut bufs.itoa,
        &mut bufs.zmij,
        format,
        q_name,
        hit.output_q_start,
        hit.output_q_end,
        t_name,
        hit.output_t_start,
        hit.output_t_end,
        char::from(hit.strand),
        hit.energy.as_f64(),
        hit.alignment.as_ref(),
        (&hit.flank_5, 0..hit.flank_5.len(), false),
        (&hit.flank_3, 0..hit.flank_3.len(), false),
        None,
    );
    writer.write_all(&bufs.line)
}

impl SearchHit {
    pub fn write(
        &self,
        w: &mut dyn Write,
        query_registry: &QueryRegistry,
        target_registry: &TargetRegistry,
    ) -> std::io::Result<()> {
        self.write_with_format(
            w,
            OutputFormat::BindingSite,
            query_registry,
            target_registry,
        )
    }

    pub fn write_with_format(
        &self,
        w: &mut dyn Write,
        format: OutputFormat,
        query_registry: &QueryRegistry,
        target_registry: &TargetRegistry,
    ) -> std::io::Result<()> {
        let mut bufs = OutputBuffers::new();
        write_hit_with_format(&mut bufs, self, format, w, query_registry, target_registry)
    }
}
