use std::io::Write;

use crate::alignment::{Alignment, PairClass};
use crate::config::OutputFormat;
use crate::index::store::TargetStore;
use crate::registry::QueryRegistry;
use crate::search::SearchHit;
use crate::seq::SeqView;
use crate::types::Base;

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

/// Append a formatted hit line to an output Vec, choosing between minimal and full formats.
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn append_hit_names_vec(
    bufs: &mut OutputBuffers,
    hit: &SearchHit,
    format: OutputFormat,
    out: &mut Vec<u8>,
    query_name: &str,
    target_name: &str,
    q_seq: SeqView<'_>,
    t_fwd: SeqView<'_>,
    t_rc: SeqView<'_>,
) {
    if format == OutputFormat::Minimal {
        append_hit_minimal_names_vec(bufs, hit, out, query_name, target_name);
    } else {
        build_line(
            &mut bufs.line,
            &mut bufs.itoa,
            hit,
            query_name,
            target_name,
            q_seq,
            t_fwd,
            t_rc,
            format,
        );
        out.extend_from_slice(&bufs.line);
    }
}

/// Append one minimal-format hit line directly into an output Vec.
///
/// This avoids the generic field loop and intermediate line buffer copy used by
/// `write_hit_names`, and is intended for high-volume minimal output paths.
#[inline]
pub fn append_hit_minimal_names_vec(
    bufs: &mut OutputBuffers,
    hit: &SearchHit,
    out: &mut Vec<u8>,
    query_name: &str,
    target_name: &str,
) {
    // q_id, q_start, q_end, t_id, t_start, t_end, strand, energy (8 fields + 7 tabs + '\n')
    out.reserve(query_name.len() + target_name.len() + 72);
    out.extend_from_slice(query_name.as_bytes());
    out.push(b'\t');
    out.extend_from_slice(bufs.itoa.format(hit.q_start + 1).as_bytes());
    out.push(b'\t');
    out.extend_from_slice(bufs.itoa.format(hit.q_end + 1).as_bytes());
    out.push(b'\t');
    out.extend_from_slice(target_name.as_bytes());
    out.push(b'\t');
    out.extend_from_slice(bufs.itoa.format(hit.t_start + 1).as_bytes());
    out.push(b'\t');
    out.extend_from_slice(bufs.itoa.format(hit.t_end + 1).as_bytes());
    out.push(b'\t');
    out.push(char::from(hit.strand) as u8);
    out.push(b'\t');
    append_score_2dp(out, &mut bufs.itoa, hit.energy.as_f64());
    out.push(b'\n');
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
fn push_alignment_mapped(buf: &mut Vec<u8>, alignment: &Alignment, map: fn(PairClass) -> u8) {
    for &p in alignment.steps() {
        buf.push(map(p));
    }
}

#[inline]
fn push_alignment_query_seq(buf: &mut Vec<u8>, alignment: &Alignment, q_bases: &[Base]) {
    let mut q_idx = 0usize;
    for &step in alignment.steps() {
        if step.consumes_query() {
            let b = q_bases.get(q_idx).copied().unwrap_or(Base::Gap);
            buf.push(b.to_byte());
            q_idx += 1;
        } else {
            buf.push(Base::Gap.to_byte());
        }
    }
}

#[inline]
fn push_alignment_target_seq(buf: &mut Vec<u8>, alignment: &Alignment, t_bases: &[Base]) {
    let mut t_idx = 0usize;
    for &step in alignment.steps() {
        if step.consumes_target() {
            let b = t_bases.get(t_idx).copied().unwrap_or(Base::Gap);
            buf.push(b.to_byte());
            t_idx += 1;
        } else {
            buf.push(Base::Gap.to_byte());
        }
    }
}

#[inline]
fn push_base_as_rna_lower_complement(buf: &mut Vec<u8>, b: Base) {
    let out = match b {
        Base::A => b'u',
        Base::G => b'c',
        Base::C => b'g',
        Base::U => b'a',
        Base::N => b'n',
        Base::Gap => b'-',
    };
    buf.push(out);
}

#[inline]
fn push_bases_as_rna_lower_complement(buf: &mut Vec<u8>, s: SeqView<'_>, reverse: bool) {
    let s = s.as_slice();
    if reverse {
        for &b in s.iter().rev() {
            push_base_as_rna_lower_complement(buf, b);
        }
    } else {
        for &b in s {
            push_base_as_rna_lower_complement(buf, b);
        }
    }
}

/// Emit binding-site target track in C `-p3` orientation.
///
/// The current internal target track is opposite-orientation transformed-space.
/// C expects the reverse-complemented target track (RNA lowercase).
#[inline]
fn push_alignment_target_seq_bindingsite(
    buf: &mut Vec<u8>,
    alignment: &Alignment,
    t_bases: &[Base],
) {
    let mut t_idx = alignment
        .steps()
        .iter()
        .filter(|step| step.consumes_target())
        .count();

    for &step in alignment.steps().iter().rev() {
        if step.consumes_target() {
            t_idx = t_idx.saturating_sub(1);
            let b = t_bases.get(t_idx).copied().unwrap_or(Base::Gap);
            push_base_as_rna_lower_complement(buf, b);
        } else {
            buf.push(b'-');
        }
    }
}

#[inline]
fn push_alignment_line(buf: &mut Vec<u8>, alignment: &Alignment) {
    push_alignment_mapped(buf, alignment, |p| p.alignment_symbol() as u8);
}

#[inline]
fn push_pairing_string(buf: &mut Vec<u8>, alignment: &Alignment) {
    push_alignment_mapped(buf, alignment, |p| p.symbol() as u8);
}

fn hit_query_bases<'a>(hit: &SearchHit, q_seq: &'a [Base]) -> &'a [Base] {
    let start = hit.q_start.min(q_seq.len());
    let end_excl = hit.q_end.saturating_add(1).min(q_seq.len());
    if end_excl < start {
        &q_seq[0..0]
    } else {
        &q_seq[start..end_excl]
    }
}

fn hit_target_bases<'a>(hit: &SearchHit, t_fwd: &'a [Base], t_rc: &'a [Base]) -> &'a [Base] {
    match hit.strand {
        crate::types::Strand::Forward => {
            let start = hit.t_start.min(t_fwd.len());
            let end_excl = hit.t_end.saturating_add(1).min(t_fwd.len());
            if end_excl < start {
                &t_fwd[0..0]
            } else {
                &t_fwd[start..end_excl]
            }
        }
        crate::types::Strand::Reverse => {
            let len = t_fwd.len();
            if len == 0 {
                return &t_rc[0..0];
            }
            let rc_start = len.saturating_sub(hit.t_end.saturating_add(1));
            let rc_end_incl = len.saturating_sub(hit.t_start.saturating_add(1));
            let start = rc_start.min(t_rc.len());
            let end_excl = rc_end_incl.saturating_add(1).min(t_rc.len());
            if end_excl < start {
                &t_rc[0..0]
            } else {
                &t_rc[start..end_excl]
            }
        }
    }
}

const BINDING_SITE_FLANK_LEN: usize = 20;

#[inline]
fn hit_target_flanks<'a>(
    hit: &SearchHit,
    t_fwd: &'a [Base],
    t_rc: &'a [Base],
) -> (&'a [Base], bool, &'a [Base], bool) {
    let len = t_fwd.len();
    if len == 0 {
        return (&t_fwd[0..0], false, &t_fwd[0..0], false);
    }

    let (oriented, start, end) = match hit.strand {
        crate::types::Strand::Forward => {
            let start = hit.t_start.min(len);
            let end = hit.t_end.min(len.saturating_sub(1));
            (t_fwd, start, end)
        }
        crate::types::Strand::Reverse => {
            if t_rc.is_empty() {
                return (&t_fwd[0..0], false, &t_fwd[0..0], false);
            }
            let start = len
                .saturating_sub(hit.t_end.saturating_add(1))
                .min(t_rc.len());
            let end = len
                .saturating_sub(hit.t_start.saturating_add(1))
                .min(t_rc.len().saturating_sub(1));
            (t_rc, start, end)
        }
    };

    if start >= oriented.len() || end >= oriented.len() || start > end {
        return (&oriented[0..0], false, &oriented[0..0], false);
    }

    // In binding-site output, flank_5 is emitted "outward" from the interaction
    // end, and flank_3 from the interaction start.
    let right_start = end.saturating_add(1).min(oriented.len());
    let right_end = right_start
        .saturating_add(BINDING_SITE_FLANK_LEN)
        .min(oriented.len());
    let left_end = start;
    let left_start = left_end.saturating_sub(BINDING_SITE_FLANK_LEN);

    let flank_5 = &oriented[right_start..right_end];
    let flank_3 = &oriented[left_start..left_end];
    (flank_5, false, flank_3, true)
}

#[allow(clippy::too_many_arguments)]
fn build_line(
    line_buf: &mut Vec<u8>,
    itoa_buf: &mut itoa::Buffer,
    hit: &SearchHit,
    q_id: &str,
    t_id: &str,
    q_seq: SeqView<'_>,
    t_fwd: SeqView<'_>,
    t_rc: SeqView<'_>,
    format: OutputFormat,
) {
    let q_seq = q_seq.as_slice();
    let t_fwd = t_fwd.as_slice();
    let t_rc = t_rc.as_slice();
    let spec = format_spec(format);
    let alignment = hit.alignment.as_ref();
    let steps_len = alignment.map(|a| a.steps().len()).unwrap_or(0);
    let (flank_5, flank_5_rev, flank_3, flank_3_rev) = hit_target_flanks(hit, t_fwd, t_rc);

    let approx = q_id.len() + t_id.len() + (steps_len * 2) + flank_5.len() + flank_3.len() + 96;
    line_buf.clear();
    if line_buf.capacity() < approx {
        line_buf.reserve(approx - line_buf.capacity());
    }

    if matches!(spec.prelude, PreludeKind::DetailedAlignment) {
        if let Some(align) = alignment {
            let q_bases = hit_query_bases(hit, q_seq);
            let t_bases = hit_target_bases(hit, t_fwd, t_rc);
            push_alignment_query_seq(line_buf, align, q_bases);
            line_buf.push(b'\n');
            push_alignment_line(line_buf, align);
            line_buf.push(b'\n');
            push_alignment_target_seq(line_buf, align, t_bases);
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
                    let t_bases = hit_target_bases(hit, t_fwd, t_rc);
                    push_alignment_target_seq_bindingsite(line_buf, align, t_bases);
                }
            }
            // C `-p3` semantics:
            // - flank_5 is complemented output from our current flank_3 side.
            // - flank_3 is complemented output from our current flank_5 side.
            FieldKind::Flank5 => {
                push_bases_as_rna_lower_complement(line_buf, SeqView::from(flank_3), flank_3_rev)
            }
            FieldKind::Flank3 => {
                push_bases_as_rna_lower_complement(line_buf, SeqView::from(flank_5), flank_5_rev)
            }
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
    target_store: &TargetStore,
) -> std::io::Result<()> {
    let q_name = query_registry.get_name(hit.query_idx);
    let t_name = target_store.get_name(hit.target_idx);
    let q_seq = query_registry.get(hit.query_idx).sequence();
    let t_idx = hit.target_idx as usize;
    let (_, t_fwd, t_rc, _) = target_store
        .target_seqs(t_idx)
        .map_err(std::io::Error::other)?;
    write_hit_names(
        bufs,
        hit,
        format,
        writer,
        q_name,
        t_name,
        q_seq,
        SeqView::from(t_fwd),
        SeqView::from(t_rc),
    )
}

/// Write a hit with explicit query/target names.
#[allow(clippy::too_many_arguments)]
pub fn write_hit_names<W: Write + ?Sized>(
    bufs: &mut OutputBuffers,
    hit: &SearchHit,
    format: OutputFormat,
    writer: &mut W,
    query_name: &str,
    target_name: &str,
    q_seq: SeqView<'_>,
    t_fwd: SeqView<'_>,
    t_rc: SeqView<'_>,
) -> std::io::Result<()> {
    bufs.line.clear();
    build_line(
        &mut bufs.line,
        &mut bufs.itoa,
        hit,
        query_name,
        target_name,
        q_seq,
        t_fwd,
        t_rc,
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
        target_store: &TargetStore,
    ) -> std::io::Result<()> {
        let mut bufs = OutputBuffers::new();
        write_hit(&mut bufs, self, format, w, query_registry, target_store)
    }
}
