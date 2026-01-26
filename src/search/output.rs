use std::io::Write;

use anyhow::Result;
use log::debug;

use crate::config::OutputFormat;
use crate::seq::utils::push_bytes_as_rna;
use crate::types::{Alignment, Pairing};

use super::{SearchHit, SearchStage};

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
pub(super) fn push_alignment_line(buf: &mut Vec<u8>, steps: &[Pairing]) {
    for p in steps {
        let c = match p {
            Pairing::Match(_, _) => b'|',
            Pairing::Wobble(_, _) => b':',
            _ => b' ',
        };
        buf.push(c);
    }
}

#[inline]
pub(super) fn push_query_seq(buf: &mut Vec<u8>, steps: &[Pairing]) {
    for p in steps {
        let c = p.query_char();
        let normalized = match c {
            'T' => 'U',
            't' => 'u',
            other => other,
        };
        buf.push(normalized as u8);
    }
}

#[inline]
pub(super) fn push_target_seq(buf: &mut Vec<u8>, steps: &[Pairing]) {
    for p in steps {
        let c = p.target_char();
        let normalized = match c {
            'T' => 'U',
            't' => 'u',
            other => other,
        };
        buf.push(normalized as u8);
    }
}

#[inline]
fn truncate_id(id: &str, max_len: Option<usize>) -> &str {
    let base = id.split_whitespace().next().unwrap_or(id);
    if let Some(max) = max_len {
        if base.len() > max { &base[..max] } else { base }
    } else {
        base
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn fill_line_buf(
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
    alignment: &Alignment,
    flank_5: (&[u8], bool),
    flank_3: (&[u8], bool),
    id_max_len: Option<usize>,
) {
    let q_id_trunc = truncate_id(q_id, id_max_len);
    let t_id_trunc = truncate_id(t_id, id_max_len);
    let steps = alignment.steps();

    let approx = q_id_trunc.len()
        + t_id_trunc.len()
        + (steps.len() * 2)
        + flank_5.0.len()
        + flank_3.0.len()
        + 96;
    line_buf.clear();
    if line_buf.capacity() < approx {
        line_buf.reserve(approx - line_buf.capacity());
    }

    match format {
        OutputFormat::Detailed => {
            push_query_seq(line_buf, steps);
            line_buf.push(b'\n');
            push_alignment_line(line_buf, steps);
            line_buf.push(b'\n');
            push_target_seq(line_buf, steps);
            line_buf.push(b'\n');
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
            line_buf.push(b'\n');
        }
        OutputFormat::Cigar => {
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
            line_buf.push(b'\t');
            for p in steps {
                line_buf.push(p.to_char() as u8);
            }
            line_buf.push(b'\n');
        }
        OutputFormat::BindingSite => {
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
            line_buf.push(b'\t');
            for p in steps {
                line_buf.push(p.to_char() as u8);
            }
            line_buf.push(b'\t');
            push_target_seq(line_buf, steps);
            line_buf.push(b'\t');
            push_bytes_as_rna(line_buf, flank_5.0, flank_5.1);
            line_buf.push(b'\t');
            push_bytes_as_rna(line_buf, flank_3.0, flank_3.1);
            line_buf.push(b'\n');
        }
        OutputFormat::Minimal => {
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
            line_buf.push(b'\n');
        }
    }
}

impl SearchHit {
    pub fn write(&self, w: &mut dyn Write) -> std::io::Result<()> {
        self.write_with_format(w, OutputFormat::BindingSite)
    }

    pub fn write_with_format(
        &self,
        w: &mut dyn Write,
        format: OutputFormat,
    ) -> std::io::Result<()> {
        let mut line_buf = Vec::new();
        let mut itoa_buf = itoa::Buffer::new();
        let mut zmij_buf = zmij::Buffer::new();
        fill_line_buf(
            &mut line_buf,
            &mut itoa_buf,
            &mut zmij_buf,
            format,
            self.query_id.as_str(),
            self.output_q_start,
            self.output_q_end,
            self.target_id.as_str(),
            self.output_t_start,
            self.output_t_end,
            char::from(self.strand),
            self.energy.as_f64(),
            &self.alignment,
            (self.flank_5.as_bytes(), false),
            (self.flank_3.as_bytes(), false),
            None,
        );
        w.write_all(&line_buf)
    }
}

pub fn write_results_to<W: Write>(hits: &[SearchHit], writer: &mut W) -> Result<()> {
    debug!("{} output=<writer>", SearchStage::Output);
    for hit in hits {
        hit.write(writer)?;
    }
    Ok(())
}

pub fn write_results_with_format_to<W: Write>(
    hits: &[SearchHit],
    writer: &mut W,
    format: OutputFormat,
) -> Result<()> {
    debug!(
        "{} output=<writer> format={:?}",
        SearchStage::Output,
        format
    );
    for hit in hits {
        hit.write_with_format(writer, format)?;
    }
    Ok(())
}
