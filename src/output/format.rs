use std::io::Write;

use crate::config::OutputFormat;
use crate::registry::{QueryRegistry, TargetRegistry};
use crate::search::SearchHit;
use crate::seq::utils::push_bases_as_rna;

/// Reusable buffers for hit formatting (avoids per-hit allocation).
pub struct OutputBuffers {
    line: Vec<u8>,
    itoa: itoa::Buffer,
}

impl OutputBuffers {
    pub fn new() -> Self {
        Self {
            line: Vec::with_capacity(256), // Preallocate typical line size
            itoa: itoa::Buffer::new(),
        }
    }
}

impl Default for OutputBuffers {
    fn default() -> Self {
        Self::new()
    }
}

#[inline]
fn append_score_2dp(buf: &mut Vec<u8>, itoa_buf: &mut itoa::Buffer, score: f64) {
    let scaled = (score * 100.0).round_ties_even() as i64;
    if scaled < 0 {
        buf.push(b'-');
    }
    let v = scaled.unsigned_abs();
    let int_part = v / 100;
    let frac = (v % 100) as u8;
    buf.extend_from_slice(itoa_buf.format(int_part).as_bytes());
    buf.push(b'.');
    buf.push(b'0' + (frac / 10));
    buf.push(b'0' + (frac % 10));
}

fn build_line(
    line_buf: &mut Vec<u8>,
    itoa_buf: &mut itoa::Buffer,
    hit: &SearchHit,
    q_id: &str,
    t_id: &str,
    format: OutputFormat,
) {
    let alignment = hit.alignment.as_ref();
    let steps_len = alignment.map(|a| a.steps().len()).unwrap_or(0);
    let flank_5 = (&hit.flank_5, 0..hit.flank_5.len(), false);
    let flank_3 = (&hit.flank_3, 0..hit.flank_3.len(), false);

    let approx = q_id.len()
        + t_id.len()
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

    line_buf.extend_from_slice(q_id.as_bytes());
    line_buf.push(b'\t');
    line_buf.extend_from_slice(itoa_buf.format(hit.output_q_start).as_bytes());
    line_buf.push(b'\t');
    line_buf.extend_from_slice(itoa_buf.format(hit.output_q_end).as_bytes());
    line_buf.push(b'\t');
    line_buf.extend_from_slice(t_id.as_bytes());
    line_buf.push(b'\t');
    line_buf.extend_from_slice(itoa_buf.format(hit.output_t_start).as_bytes());
    line_buf.push(b'\t');
    line_buf.extend_from_slice(itoa_buf.format(hit.output_t_end).as_bytes());
    line_buf.push(b'\t');
    line_buf.push(char::from(hit.strand) as u8);
    line_buf.push(b'\t');
    append_score_2dp(line_buf, itoa_buf, hit.energy.as_f64());

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

/// Write a hit with format using provided reusable buffers.
pub fn write_hit<W: Write + ?Sized>(
    bufs: &mut OutputBuffers,
    hit: &SearchHit,
    format: OutputFormat,
    writer: &mut W,
    query_registry: &QueryRegistry,
    target_registry: &TargetRegistry,
) -> std::io::Result<()> {
    bufs.line.clear();
    let q_name = query_registry.get_name(hit.query_idx);
    let t_name = target_registry.get_name(hit.target_idx);
    build_line(&mut bufs.line, &mut bufs.itoa, hit, q_name, t_name, format);
    writer.write_all(&bufs.line)
}

impl SearchHit {
    pub fn write(
        &self,
        w: &mut dyn Write,
        format: OutputFormat,
        query_registry: &QueryRegistry,
        target_registry: &TargetRegistry,
    ) -> std::io::Result<()> {
        let mut bufs = OutputBuffers::new();
        write_hit(&mut bufs, self, format, w, query_registry, target_registry)
    }
}
