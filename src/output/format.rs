use std::io::Write;

use crate::alignment::{Alignment, Pairing};
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

#[derive(Clone, Copy)]
enum FieldKind {
    QueryId,
    QStart,
    QEnd,
    TargetId,
    TStart,
    TEnd,
    Strand,
    Energy2dp,
    Pairing,
    TargetSeq,
    Flank5,
    Flank3,
}

#[derive(Clone, Copy)]
enum PreludeKind {
    None,
    DetailedAlignment,
}

struct FormatSpec {
    prelude: PreludeKind,
    fields: &'static [FieldKind],
}

const FIELDS_BASE: &[FieldKind] = &[
    FieldKind::QueryId,
    FieldKind::QStart,
    FieldKind::QEnd,
    FieldKind::TargetId,
    FieldKind::TStart,
    FieldKind::TEnd,
    FieldKind::Strand,
    FieldKind::Energy2dp,
];

const FIELDS_CIGAR: &[FieldKind] = &[
    FieldKind::QueryId,
    FieldKind::QStart,
    FieldKind::QEnd,
    FieldKind::TargetId,
    FieldKind::TStart,
    FieldKind::TEnd,
    FieldKind::Strand,
    FieldKind::Energy2dp,
    FieldKind::Pairing,
];

const FIELDS_BINDING_SITE: &[FieldKind] = &[
    FieldKind::QueryId,
    FieldKind::QStart,
    FieldKind::QEnd,
    FieldKind::TargetId,
    FieldKind::TStart,
    FieldKind::TEnd,
    FieldKind::Strand,
    FieldKind::Energy2dp,
    FieldKind::Pairing,
    FieldKind::TargetSeq,
    FieldKind::Flank5,
    FieldKind::Flank3,
];

#[inline]
fn format_spec(format: OutputFormat) -> FormatSpec {
    match format {
        OutputFormat::Minimal => FormatSpec {
            prelude: PreludeKind::None,
            fields: FIELDS_BASE,
        },
        OutputFormat::Detailed => FormatSpec {
            prelude: PreludeKind::DetailedAlignment,
            fields: FIELDS_BASE,
        },
        OutputFormat::Cigar => FormatSpec {
            prelude: PreludeKind::None,
            fields: FIELDS_CIGAR,
        },
        OutputFormat::BindingSite => FormatSpec {
            prelude: PreludeKind::None,
            fields: FIELDS_BINDING_SITE,
        },
    }
}

#[inline]
fn push_alignment_mapped(buf: &mut Vec<u8>, alignment: &Alignment, map: fn(Pairing) -> u8) {
    for &p in alignment.steps() {
        buf.push(map(p));
    }
}

#[inline]
fn push_alignment_query_seq(buf: &mut Vec<u8>, alignment: &Alignment) {
    push_alignment_mapped(buf, alignment, |p| p.query_char() as u8);
}

#[inline]
fn push_alignment_target_seq(buf: &mut Vec<u8>, alignment: &Alignment) {
    push_alignment_mapped(buf, alignment, |p| p.target_char() as u8);
}

#[inline]
fn push_alignment_line(buf: &mut Vec<u8>, alignment: &Alignment) {
    push_alignment_mapped(buf, alignment, |p| p.class().alignment_symbol() as u8);
}

#[inline]
fn push_pairing_string(buf: &mut Vec<u8>, alignment: &Alignment) {
    push_alignment_mapped(buf, alignment, |p| p.class().symbol() as u8);
}

fn build_line(
    line_buf: &mut Vec<u8>,
    itoa_buf: &mut itoa::Buffer,
    hit: &SearchHit,
    q_id: &str,
    t_id: &str,
    format: OutputFormat,
) {
    let spec = format_spec(format);
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

    if matches!(spec.prelude, PreludeKind::DetailedAlignment) {
        if let Some(align) = alignment {
            push_alignment_query_seq(line_buf, align);
            line_buf.push(b'\n');
            push_alignment_line(line_buf, align);
            line_buf.push(b'\n');
            push_alignment_target_seq(line_buf, align);
            line_buf.push(b'\n');
        }
    }

    for (idx, field) in spec.fields.iter().enumerate() {
        if idx > 0 {
            line_buf.push(b'\t');
        }
        match field {
            FieldKind::QueryId => line_buf.extend_from_slice(q_id.as_bytes()),
            FieldKind::QStart => {
                line_buf.extend_from_slice(itoa_buf.format(hit.q_start + 1).as_bytes())
            }
            FieldKind::QEnd => {
                line_buf.extend_from_slice(itoa_buf.format(hit.q_end + 1).as_bytes())
            }
            FieldKind::TargetId => line_buf.extend_from_slice(t_id.as_bytes()),
            FieldKind::TStart => {
                line_buf.extend_from_slice(itoa_buf.format(hit.t_start + 1).as_bytes())
            }
            FieldKind::TEnd => {
                line_buf.extend_from_slice(itoa_buf.format(hit.t_end + 1).as_bytes())
            }
            FieldKind::Strand => line_buf.push(char::from(hit.strand) as u8),
            FieldKind::Energy2dp => append_score_2dp(line_buf, itoa_buf, hit.energy.as_f64()),
            FieldKind::Pairing => {
                if let Some(align) = alignment {
                    push_pairing_string(line_buf, align);
                }
            }
            FieldKind::TargetSeq => {
                if let Some(align) = alignment {
                    push_alignment_target_seq(line_buf, align);
                }
            }
            FieldKind::Flank5 => push_bases_as_rna(
                line_buf,
                &flank_5.0[flank_5.1.start..flank_5.1.end],
                flank_5.2,
            ),
            FieldKind::Flank3 => push_bases_as_rna(
                line_buf,
                &flank_3.0[flank_3.1.start..flank_3.1.end],
                flank_3.2,
            ),
        }
    }
    line_buf.push(b'\n');
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
    let q_name = query_registry.get_name(hit.query_idx);
    let t_name = target_registry.get_name(hit.target_idx);
    write_hit_names(bufs, hit, format, writer, q_name, t_name)
}

/// Write a hit with explicit query/target names.
pub fn write_hit_names<W: Write + ?Sized>(
    bufs: &mut OutputBuffers,
    hit: &SearchHit,
    format: OutputFormat,
    writer: &mut W,
    query_name: &str,
    target_name: &str,
) -> std::io::Result<()> {
    bufs.line.clear();
    build_line(
        &mut bufs.line,
        &mut bufs.itoa,
        hit,
        query_name,
        target_name,
        format,
    );
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
