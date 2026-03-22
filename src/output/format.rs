use std::io::Write;

use crate::alignment::{Alignment, PairClass};
use crate::config::OutputFormat;
use crate::index::store::TargetStore;
use crate::registry::QueryRegistry;
use crate::search::SearchHit;
use crate::types::{Base, Strand};

/// Resolved context for formatting a single hit: names and raw sequences.
///
/// All fields are borrows, so `HitCtx` is `Copy` and free to pass by value.
#[derive(Clone, Copy)]
pub struct HitCtx<'a> {
    pub q_name: &'a str,
    pub q_seq: &'a [Base],
    pub t_name: &'a str,
    pub t_fwd: &'a [Base],
    pub t_rc: &'a [Base],
}

// =============================================================================
// FORMAT SCHEMA
// =============================================================================

#[derive(Clone, Copy)]
enum FieldKind {
    QueryId,
    QStart,
    QEnd,
    TargetId,
    TStart,
    TEnd,
    Strand,
    Energy,
    Pairing,
    TargetSeq,
    Flank5,
    Flank3,
}

#[derive(Clone, Copy, PartialEq)]
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
    FieldKind::Energy,
];

const FIELDS_CIGAR: &[FieldKind] = &[
    FieldKind::QueryId,
    FieldKind::QStart,
    FieldKind::QEnd,
    FieldKind::TargetId,
    FieldKind::TStart,
    FieldKind::TEnd,
    FieldKind::Strand,
    FieldKind::Energy,
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
    FieldKind::Energy,
    FieldKind::Pairing,
    FieldKind::TargetSeq,
    FieldKind::Flank5,
    FieldKind::Flank3,
];

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

// =============================================================================
// HELPERS
// =============================================================================

#[inline]
fn append_score_2dp(buf: &mut Vec<u8>, itoa: &mut itoa::Buffer, score: f64) {
    let scaled = (score * 100.0).round_ties_even() as i64;
    if scaled < 0 {
        buf.push(b'-');
    }
    let v = scaled.unsigned_abs();
    buf.extend_from_slice(itoa.format(v / 100).as_bytes());
    buf.push(b'.');
    let frac = (v % 100) as u8;
    buf.push(b'0' + (frac / 10));
    buf.push(b'0' + (frac % 10));
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
            buf.push(q_bases.get(q_idx).copied().unwrap_or(Base::Gap).to_byte());
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
            buf.push(t_bases.get(t_idx).copied().unwrap_or(Base::Gap).to_byte());
            t_idx += 1;
        } else {
            buf.push(Base::Gap.to_byte());
        }
    }
}

/// Emit the target track in C `-p3` orientation (reverse walk, RNA lowercase complement).
///
/// The internal target track is in opposite-orientation transformed space.
/// C expects the reverse-complemented track (RNA lowercase).
#[inline]
fn push_alignment_target_seq_bindingsite(
    buf: &mut Vec<u8>,
    alignment: &Alignment,
    t_bases: &[Base],
) {
    let mut t_idx = alignment
        .steps()
        .iter()
        .filter(|s| s.consumes_target())
        .count();
    for &step in alignment.steps().iter().rev() {
        if step.consumes_target() {
            t_idx = t_idx.saturating_sub(1);
            buf.push(
                t_bases
                    .get(t_idx)
                    .copied()
                    .unwrap_or(Base::Gap)
                    .complement()
                    .to_byte(),
            );
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
    let end = hit.q_end.saturating_add(1).min(q_seq.len());
    if end < start {
        &q_seq[0..0]
    } else {
        &q_seq[start..end]
    }
}

fn hit_target_bases<'a>(hit: &SearchHit, t_fwd: &'a [Base], t_rc: &'a [Base]) -> &'a [Base] {
    match hit.strand {
        Strand::Forward => {
            let start = hit.t_start.min(t_fwd.len());
            let end = hit.t_end.saturating_add(1).min(t_fwd.len());
            if end < start {
                &t_fwd[0..0]
            } else {
                &t_fwd[start..end]
            }
        }
        Strand::Reverse => {
            let len = t_fwd.len();
            if len == 0 {
                return &t_rc[0..0];
            }
            let start = len
                .saturating_sub(hit.t_end.saturating_add(1))
                .min(t_rc.len());
            let end = len
                .saturating_sub(hit.t_start.saturating_add(1))
                .saturating_add(1)
                .min(t_rc.len());
            if end < start {
                &t_rc[0..0]
            } else {
                &t_rc[start..end]
            }
        }
    }
}

const BINDING_SITE_FLANK_LEN: usize = 20;

/// Returns `(flank_5, flank_5_rev, flank_3, flank_3_rev)` in oriented-strand space.
///
/// C `-p3` naming is swapped vs. internal: output `flank5` = our `flank_3` complemented.
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
        Strand::Forward => (
            t_fwd,
            hit.t_start.min(len),
            hit.t_end.min(len.saturating_sub(1)),
        ),
        Strand::Reverse => {
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

    let right_start = end.saturating_add(1).min(oriented.len());
    let right_end = right_start
        .saturating_add(BINDING_SITE_FLANK_LEN)
        .min(oriented.len());
    let left_end = start;
    let left_start = left_end.saturating_sub(BINDING_SITE_FLANK_LEN);

    (
        &oriented[right_start..right_end],
        false,
        &oriented[left_start..left_end],
        true,
    )
}

// =============================================================================
// CORE FORMATTER
// =============================================================================

/// Format a hit line and append it directly to `out`.
///
/// Writes directly into the caller's buffer with no intermediate copy.
/// Flanks are computed only for `OutputFormat::BindingSite`.
#[inline]
pub fn format_hit_into(
    out: &mut Vec<u8>,
    itoa: &mut itoa::Buffer,
    hit: &SearchHit,
    ctx: HitCtx<'_>,
    format: OutputFormat,
) {
    let HitCtx {
        q_name,
        q_seq,
        t_name,
        t_fwd,
        t_rc,
    } = ctx;

    // Minimal fast path: 8 direct writes, no field loop or match dispatch.
    if format == OutputFormat::Minimal {
        out.reserve(q_name.len() + t_name.len() + 72);
        out.extend_from_slice(q_name.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(itoa.format(hit.q_start + 1).as_bytes());
        out.push(b'\t');
        out.extend_from_slice(itoa.format(hit.q_end + 1).as_bytes());
        out.push(b'\t');
        out.extend_from_slice(t_name.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(itoa.format(hit.t_start + 1).as_bytes());
        out.push(b'\t');
        out.extend_from_slice(itoa.format(hit.t_end + 1).as_bytes());
        out.push(b'\t');
        out.push(char::from(hit.strand) as u8);
        out.push(b'\t');
        append_score_2dp(out, itoa, hit.energy.as_f64());
        out.push(b'\n');
        return;
    }

    let spec = format_spec(format);
    let alignment = hit.alignment.as_ref();
    let steps_len = alignment.map(|a| a.steps().len()).unwrap_or(0);

    // Flanks are only needed (and only computed) for BindingSite format.
    let flanks = if format == OutputFormat::BindingSite {
        Some(hit_target_flanks(hit, t_fwd, t_rc))
    } else {
        None
    };

    let flank_len = flanks
        .map(|(f5, _, f3, _)| f5.len() + f3.len())
        .unwrap_or(0);
    out.reserve(q_name.len() + t_name.len() + (steps_len * 2) + flank_len + 96);

    if spec.prelude == PreludeKind::DetailedAlignment {
        if let Some(align) = alignment {
            let q_bases = hit_query_bases(hit, q_seq);
            let t_bases = hit_target_bases(hit, t_fwd, t_rc);
            push_alignment_query_seq(out, align, q_bases);
            out.push(b'\n');
            push_alignment_line(out, align);
            out.push(b'\n');
            push_alignment_target_seq(out, align, t_bases);
            out.push(b'\n');
        }
    }

    for (idx, field) in spec.fields.iter().enumerate() {
        if idx > 0 {
            out.push(b'\t');
        }
        match field {
            FieldKind::QueryId => out.extend_from_slice(q_name.as_bytes()),
            FieldKind::QStart => out.extend_from_slice(itoa.format(hit.q_start + 1).as_bytes()),
            FieldKind::QEnd => out.extend_from_slice(itoa.format(hit.q_end + 1).as_bytes()),
            FieldKind::TargetId => out.extend_from_slice(t_name.as_bytes()),
            FieldKind::TStart => out.extend_from_slice(itoa.format(hit.t_start + 1).as_bytes()),
            FieldKind::TEnd => out.extend_from_slice(itoa.format(hit.t_end + 1).as_bytes()),
            FieldKind::Strand => out.push(char::from(hit.strand) as u8),
            FieldKind::Energy => append_score_2dp(out, itoa, hit.energy.as_f64()),
            FieldKind::Pairing => {
                if let Some(align) = alignment {
                    push_pairing_string(out, align);
                }
            }
            FieldKind::TargetSeq => {
                if let Some(align) = alignment {
                    push_alignment_target_seq_bindingsite(
                        out,
                        align,
                        hit_target_bases(hit, t_fwd, t_rc),
                    );
                }
            }
            // C `-p3` semantics: output `flank5` = our `flank_3` complemented, and vice-versa.
            FieldKind::Flank5 => {
                if let Some((_, _, flank_3, flank_3_rev)) = flanks {
                    if flank_3_rev {
                        for &b in flank_3.iter().rev() {
                            out.push(b.complement().to_byte());
                        }
                    } else {
                        for &b in flank_3 {
                            out.push(b.complement().to_byte());
                        }
                    }
                }
            }
            FieldKind::Flank3 => {
                if let Some((flank_5, flank_5_rev, _, _)) = flanks {
                    if flank_5_rev {
                        for &b in flank_5.iter().rev() {
                            out.push(b.complement().to_byte());
                        }
                    } else {
                        for &b in flank_5 {
                            out.push(b.complement().to_byte());
                        }
                    }
                }
            }
        }
    }
    out.push(b'\n');
}

// =============================================================================
// WRITE HELPER
// =============================================================================

/// Write a hit to a `Write` impl, resolving names and sequences from the registries.
pub fn write_hit<W: Write + ?Sized>(
    hit: &SearchHit,
    format: OutputFormat,
    writer: &mut W,
    query_registry: &QueryRegistry,
    target_store: &TargetStore,
) -> std::io::Result<()> {
    let t_idx = hit.target_idx as usize;
    let (_, t_fwd, t_rc, _) = target_store
        .target_seqs(t_idx)
        .map_err(std::io::Error::other)?;
    let ctx = HitCtx {
        q_name: query_registry.get_name(hit.query_idx),
        q_seq: query_registry.get(hit.query_idx).sequence().as_slice(),
        t_name: target_store.get_name(hit.target_idx),
        t_fwd,
        t_rc,
    };
    let mut line = Vec::new();
    format_hit_into(&mut line, &mut itoa::Buffer::new(), hit, ctx, format);
    writer.write_all(&line)
}

impl SearchHit {
    pub fn write(
        &self,
        w: &mut dyn Write,
        format: OutputFormat,
        query_registry: &QueryRegistry,
        target_store: &TargetStore,
    ) -> std::io::Result<()> {
        write_hit(self, format, w, query_registry, target_store)
    }
}
