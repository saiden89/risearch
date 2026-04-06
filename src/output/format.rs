use crate::alignment::{Alignment, PairClass};
use crate::config::OutputFormat;
use crate::search::SearchHit;
use crate::types::{Base, Strand};

/// Resolved context for formatting a single hit: names and raw sequences.
#[derive(Clone, Copy)]
pub struct HitCtx<'a> {
    pub q_name: &'a str,
    pub q_seq: &'a [Base],
    pub t_name: &'a str,
    pub t_fwd: &'a [Base],
    pub t_rc: &'a [Base],
}

/// Zero-cost TSV line builder. Auto-inserts tab separators between fields.
struct TsvLine<'a> {
    buf: &'a mut Vec<u8>,
    needs_sep: bool,
}

impl<'a> TsvLine<'a> {
    #[inline(always)]
    fn new(buf: &'a mut Vec<u8>) -> Self {
        Self {
            buf,
            needs_sep: false,
        }
    }

    #[inline(always)]
    fn field(&mut self, bytes: &[u8]) {
        if self.needs_sep {
            self.buf.push(b'\t');
        }
        self.needs_sep = true;
        self.buf.extend_from_slice(bytes);
    }

    #[inline(always)]
    fn field_with(&mut self, f: impl FnOnce(&mut Vec<u8>)) {
        if self.needs_sep {
            self.buf.push(b'\t');
        }
        self.needs_sep = true;
        f(self.buf);
    }

    #[inline(always)]
    fn finish(self) {
        self.buf.push(b'\n');
    }
}

/// Format a hit as a TSV line and append it to `out`.
#[inline]
pub fn format_hit_into(
    out: &mut Vec<u8>,
    itoa: &mut itoa::Buffer,
    hit: &SearchHit,
    ctx: HitCtx<'_>,
    format: OutputFormat,
) {
    let alignment = hit.alignment.as_ref();
    let steps_len = alignment.map(|a| a.steps().len()).unwrap_or(0);
    out.reserve(ctx.q_name.len() + ctx.t_name.len() + (steps_len * 2) + 96);

    if format == OutputFormat::Detailed {
        write_alignment_prelude(out, hit, ctx);
    }

    let mut row = TsvLine::new(out);
    write_base_fields(&mut row, itoa, hit, ctx);
    write_extended_fields(&mut row, itoa, hit, ctx, format);
    row.finish();
}

/// The 8 fields shared by every output format.
#[inline(always)]
fn write_base_fields(row: &mut TsvLine<'_>, itoa: &mut itoa::Buffer, hit: &SearchHit, ctx: HitCtx<'_>) {
    row.field(ctx.q_name.as_bytes());
    row.field(itoa.format(hit.q_start + 1).as_bytes());
    row.field(itoa.format(hit.q_end + 1).as_bytes());
    row.field(ctx.t_name.as_bytes());
    row.field(itoa.format(hit.t_start + 1).as_bytes());
    row.field(itoa.format(hit.t_end + 1).as_bytes());
    row.field(&[char::from(hit.strand) as u8]);
    row.field_with(|buf| append_score_2dp(buf, itoa, hit.energy.as_f64()));
}

/// Format-specific fields beyond the base 8.
#[inline(always)]
fn write_extended_fields(
    row: &mut TsvLine<'_>,
    _itoa: &mut itoa::Buffer,
    hit: &SearchHit,
    ctx: HitCtx<'_>,
    format: OutputFormat,
) {
    let alignment = hit.alignment.as_ref();

    if matches!(format, OutputFormat::Cigar | OutputFormat::BindingSite) {
        row.field_with(|buf| {
            if let Some(align) = alignment {
                push_pairing_string(buf, align);
            }
        });
    }

    if format == OutputFormat::BindingSite {
        row.field_with(|buf| {
            if let Some(align) = alignment {
                push_alignment_target_seq_bindingsite(
                    buf,
                    align,
                    hit_target_bases(hit, ctx.t_fwd, ctx.t_rc),
                );
            }
        });

        // C `-p3` semantics: output `flank5` = our `flank_3` complemented, and vice-versa.
        let (flank_5, flank_5_rev, flank_3, flank_3_rev) =
            hit_target_flanks(hit, ctx.t_fwd, ctx.t_rc);
        row.field_with(|buf| push_flank(buf, flank_3, flank_3_rev));
        row.field_with(|buf| push_flank(buf, flank_5, flank_5_rev));
    }
}

/// 3-line alignment visualization (Detailed format only).
fn write_alignment_prelude(out: &mut Vec<u8>, hit: &SearchHit, ctx: HitCtx<'_>) {
    let Some(align) = hit.alignment.as_ref() else {
        return;
    };
    let q_bases = hit_query_bases(hit, ctx.q_seq);
    let t_bases = hit_target_bases(hit, ctx.t_fwd, ctx.t_rc);

    // Query track
    let mut q_idx = 0usize;
    for &step in align.steps() {
        if step.consumes_query() {
            out.push(q_bases.get(q_idx).copied().unwrap_or(Base::Gap).to_byte());
            q_idx += 1;
        } else {
            out.push(Base::Gap.to_byte());
        }
    }
    out.push(b'\n');

    // Alignment symbols
    for &p in align.steps() {
        out.push(p.alignment_symbol() as u8);
    }
    out.push(b'\n');

    // Target track
    let mut t_idx = 0usize;
    for &step in align.steps() {
        if step.consumes_target() {
            out.push(t_bases.get(t_idx).copied().unwrap_or(Base::Gap).to_byte());
            t_idx += 1;
        } else {
            out.push(Base::Gap.to_byte());
        }
    }
    out.push(b'\n');
}

// ---------------------------------------------------------------------------
// Helpers: sequence slicing, flanks, score formatting
// ---------------------------------------------------------------------------

fn hit_query_bases<'a>(hit: &SearchHit, q_seq: &'a [Base]) -> &'a [Base] {
    let start = hit.q_start.min(q_seq.len());
    let end = hit.q_end.saturating_add(1).min(q_seq.len());
    if end < start { &q_seq[0..0] } else { &q_seq[start..end] }
}

fn hit_target_bases<'a>(hit: &SearchHit, t_fwd: &'a [Base], t_rc: &'a [Base]) -> &'a [Base] {
    match hit.strand {
        Strand::Forward => {
            let start = hit.t_start.min(t_fwd.len());
            let end = hit.t_end.saturating_add(1).min(t_fwd.len());
            if end < start { &t_fwd[0..0] } else { &t_fwd[start..end] }
        }
        Strand::Reverse => {
            let len = t_fwd.len();
            if len == 0 {
                return &t_rc[0..0];
            }
            let start = len.saturating_sub(hit.t_end.saturating_add(1)).min(t_rc.len());
            let end = len
                .saturating_sub(hit.t_start.saturating_add(1))
                .saturating_add(1)
                .min(t_rc.len());
            if end < start { &t_rc[0..0] } else { &t_rc[start..end] }
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
            let start = len.saturating_sub(hit.t_end.saturating_add(1)).min(t_rc.len());
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
    let right_end = right_start.saturating_add(BINDING_SITE_FLANK_LEN).min(oriented.len());
    let left_end = start;
    let left_start = left_end.saturating_sub(BINDING_SITE_FLANK_LEN);

    (
        &oriented[right_start..right_end],
        false,
        &oriented[left_start..left_end],
        true,
    )
}

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

fn push_pairing_string(buf: &mut Vec<u8>, alignment: &Alignment) {
    for &p in alignment.steps() {
        buf.push(p.symbol() as u8);
    }
}

/// Target track in C `-p3` orientation (reverse walk, RNA lowercase complement).
fn push_alignment_target_seq_bindingsite(
    buf: &mut Vec<u8>,
    alignment: &Alignment,
    t_bases: &[Base],
) {
    let mut t_idx = alignment.steps().iter().filter(|s| s.consumes_target()).count();
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

fn push_flank(buf: &mut Vec<u8>, bases: &[Base], rev: bool) {
    if rev {
        for &b in bases.iter().rev() {
            buf.push(b.complement().to_byte());
        }
    } else {
        for &b in bases {
            buf.push(b.complement().to_byte());
        }
    }
}
