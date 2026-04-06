use crate::alignment::Alignment;
use crate::config::OutputFormat;
use crate::search::SearchHit;
use crate::types::Base;

const BINDING_SITE_FLANK_LEN: usize = 20;

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

#[inline]
#[allow(clippy::too_many_arguments)]
pub fn format_hit_into(
    out: &mut Vec<u8>,
    itoa: &mut itoa::Buffer,
    hit: &SearchHit,
    q_name: &str,
    q_seq: &[Base],
    t_name: &str,
    t_fwd: &[Base],
    t_rc: &[Base],
    format: OutputFormat,
) {
    let alignment = hit.alignment.as_ref();
    let steps_len = alignment.map(|a| a.steps().len()).unwrap_or(0);
    let flank_reserve = if format == OutputFormat::BindingSite {
        2 * BINDING_SITE_FLANK_LEN
    } else {
        0
    };
    out.reserve(q_name.len() + t_name.len() + (steps_len * 2) + flank_reserve + 96);

    if format == OutputFormat::Detailed {
        write_alignment_prelude(out, hit, q_seq, t_fwd, t_rc);
    }

    let mut row = TsvLine::new(out);
    write_base_fields(&mut row, itoa, hit, q_name, t_name);
    write_extended_fields(&mut row, hit, t_fwd, t_rc, format);
    row.finish();
}

#[inline(always)]
fn write_base_fields(
    row: &mut TsvLine<'_>,
    itoa: &mut itoa::Buffer,
    hit: &SearchHit,
    q_name: &str,
    t_name: &str,
) {
    row.field(q_name.as_bytes());
    row.field(itoa.format(hit.q_start + 1).as_bytes());
    row.field(itoa.format(hit.q_end + 1).as_bytes());
    row.field(t_name.as_bytes());
    row.field(itoa.format(hit.t_start + 1).as_bytes());
    row.field(itoa.format(hit.t_end + 1).as_bytes());
    row.field(&[char::from(hit.strand) as u8]);
    row.field_with(|buf| append_score_2dp(buf, itoa, hit.energy.as_f64()));
}

#[inline(always)]
fn write_extended_fields(
    row: &mut TsvLine<'_>,
    hit: &SearchHit,
    t_fwd: &[Base],
    t_rc: &[Base],
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
                push_alignment_target_seq_bindingsite(buf, align, hit.target_bases(t_fwd, t_rc));
            }
        });

        let (flank_right, flank_left) = hit.target_flanks(t_fwd, t_rc);
        row.field_with(|buf| push_flank_rev(buf, flank_left));
        row.field_with(|buf| push_flank_fwd(buf, flank_right));
    }
}

fn write_alignment_prelude(
    out: &mut Vec<u8>,
    hit: &SearchHit,
    q_seq: &[Base],
    t_fwd: &[Base],
    t_rc: &[Base],
) {
    let Some(align) = hit.alignment.as_ref() else {
        return;
    };
    let q_bases = hit.query_bases(q_seq);
    let t_bases = hit.target_bases(t_fwd, t_rc);

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

    for &p in align.steps() {
        out.push(p.alignment_symbol() as u8);
    }
    out.push(b'\n');

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

fn push_alignment_target_seq_bindingsite(
    buf: &mut Vec<u8>,
    alignment: &Alignment,
    t_bases: &[Base],
) {
    let mut t_idx = t_bases.len();
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

fn push_flank_fwd(buf: &mut Vec<u8>, bases: &[Base]) {
    for &b in bases {
        buf.push(b.complement().to_byte());
    }
}

fn push_flank_rev(buf: &mut Vec<u8>, bases: &[Base]) {
    for &b in bases.iter().rev() {
        buf.push(b.complement().to_byte());
    }
}
