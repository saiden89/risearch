use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use log::debug;

use crate::args::OutputFormat;
use crate::types::Pairing;

use super::{SearchHit, SearchStage};

/// Convert bytes to RNA string (T->U), optionally reversed. Single-pass, one allocation.
#[inline]
pub(super) fn bytes_to_rna_string(s: &[u8], reverse: bool) -> String {
    let mut result = String::with_capacity(s.len());
    if reverse {
        for &b in s.iter().rev() {
            result.push(match b {
                b'T' => 'U',
                b't' => 'u',
                _ => b as char,
            });
        }
    } else {
        for &b in s {
            result.push(match b {
                b'T' => 'U',
                b't' => 'u',
                _ => b as char,
            });
        }
    }
    result
}

#[inline]
pub(super) fn push_bytes_as_rna(buf: &mut Vec<u8>, s: &[u8], reverse: bool) {
    if reverse {
        for &b in s.iter().rev() {
            let c = match b {
                b'T' => b'U',
                b't' => b'u',
                _ => b,
            };
            buf.push(c);
        }
    } else {
        for &b in s {
            let c = match b {
                b'T' => b'U',
                b't' => b'u',
                _ => b,
            };
            buf.push(c);
        }
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

impl SearchHit {
    pub fn write(&self, w: &mut dyn Write) -> std::io::Result<()> {
        self.write_with_format(w, OutputFormat::BindingSite)
    }

    pub fn write_with_format(
        &self,
        w: &mut dyn Write,
        format: OutputFormat,
    ) -> std::io::Result<()> {
        let q_id = self.query_id.truncated();
        let t_id = self.target_id.truncated();

        match format {
            OutputFormat::Detailed => {
                let query = self.alignment.query_sequence();
                let aln = self.alignment.alignment_string();
                let target = self.alignment.target_sequence();

                writeln!(w, "{}", query)?;
                writeln!(w, "{}", aln)?;
                writeln!(w, "{}", target)?;
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}",
                    q_id,
                    self.output_q_start,
                    self.output_q_end,
                    t_id,
                    self.output_t_start,
                    self.output_t_end,
                    self.strand,
                    self.energy.as_f64()
                )
            }
            OutputFormat::Cigar => {
                write!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t",
                    q_id,
                    self.output_q_start,
                    self.output_q_end,
                    t_id,
                    self.output_t_start,
                    self.output_t_end,
                    self.strand,
                    self.energy.as_f64()
                )?;
                for p in self.alignment.steps() {
                    write!(w, "{}", p.to_char())?;
                }
                writeln!(w)
            }
            OutputFormat::BindingSite => {
                write!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t",
                    q_id,
                    self.output_q_start,
                    self.output_q_end,
                    t_id,
                    self.output_t_start,
                    self.output_t_end,
                    self.strand,
                    self.energy.as_f64()
                )?;
                for p in self.alignment.steps() {
                    write!(w, "{}", p.to_char())?;
                }
                write!(w, "\t")?;
                for p in self.alignment.steps() {
                    let c = p.target_char();
                    let normalized = match c {
                        'T' => 'U',
                        't' => 'u',
                        other => other,
                    };
                    write!(w, "{}", normalized)?;
                }
                writeln!(w, "\t{}\t{}", self.flank_5, self.flank_3)
            }
            OutputFormat::Minimal => writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}",
                q_id,
                self.output_q_start,
                self.output_q_end,
                t_id,
                self.output_t_start,
                self.output_t_end,
                self.strand,
                self.energy.as_f64()
            ),
        }
    }
}

pub fn write_results(hits: &[SearchHit], output: impl AsRef<Path>) -> Result<()> {
    let inner: Box<dyn Write> = if output.as_ref() == Path::new("-") {
        Box::new(std::io::stdout())
    } else {
        Box::new(std::fs::File::create(output.as_ref()).context("Failed to create output file")?)
    };

    // Larger buffer reduces syscall overhead for high-volume output.
    let mut writer = BufWriter::with_capacity(256 * 1024, inner);

    debug!("{} output={:?}", SearchStage::Output, output.as_ref());

    for hit in hits {
        hit.write(&mut writer)?;
    }

    // BufWriter flushes on drop, but explicit flush ensures errors are caught
    writer.flush().context("Failed to flush output")?;
    Ok(())
}

pub fn write_results_with_format(
    hits: &[SearchHit],
    output: impl AsRef<Path>,
    format: OutputFormat,
) -> Result<()> {
    let inner: Box<dyn Write> = if output.as_ref() == Path::new("-") {
        Box::new(std::io::stdout())
    } else {
        Box::new(std::fs::File::create(output.as_ref()).context("Failed to create output file")?)
    };

    // Larger buffer reduces syscall overhead for high-volume output.
    let mut writer = BufWriter::with_capacity(256 * 1024, inner);

    debug!(
        "{} output={:?} format={:?}",
        SearchStage::Output,
        output.as_ref(),
        format
    );

    for hit in hits {
        hit.write_with_format(&mut writer, format)?;
    }

    writer.flush().context("Failed to flush output")?;
    Ok(())
}
