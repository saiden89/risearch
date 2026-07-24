use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use risearch::alignment::{Alignment, PairClass};
use risearch::config::OutputFormat;
use risearch::output::format::format_hit_into;
use risearch::search::SearchHit;
use risearch::types::{Base, Energy, Strand};

fn setup_hit() -> (SearchHit, Vec<Base>, Vec<Base>, Vec<Base>) {
    let q_bases = vec![Base::U, Base::G, Base::C, Base::U, Base::G, Base::C];
    let t_fwd = vec![Base::G, Base::U, Base::C, Base::G, Base::U, Base::C];
    let t_rc = vec![Base::G, Base::U, Base::C, Base::G, Base::U, Base::C];

    let seed: Vec<PairClass> = (0..6).map(|_| PairClass::Canonical).collect();
    let alignment = Alignment::from_parts(&[], &seed, &[]);

    let hit = SearchHit {
        query_idx: 0,
        target_idx: 0,
        q_start: 0,
        q_end: 5,
        t_start: 0,
        t_end: 5,
        strand: Strand::Forward,
        energy: Energy::from_kcal(-10.5),
        seed_start: Some(0),
        seed_end: Some(6),
        alignment: Some(alignment),
    };

    (hit, q_bases, t_fwd, t_rc)
}

fn bench_format_hit(c: &mut Criterion) {
    let (hit, q_seq, t_fwd, t_rc) = setup_hit();

    let mut group = c.benchmark_group("hit_formatting");

    for format in [
        OutputFormat::Minimal,
        OutputFormat::Detailed,
        OutputFormat::Cigar,
        OutputFormat::BindingSite,
    ] {
        let mut itoa = itoa::Buffer::new();
        let mut out = Vec::with_capacity(1024);

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{:?}", format)),
            &format,
            |b, &fmt| {
                b.iter(|| {
                    out.clear();
                    format_hit_into(
                        &mut out,
                        &mut itoa,
                        &hit,
                        "test_query",
                        &q_seq,
                        "test_target",
                        &t_fwd,
                        &t_rc,
                        fmt,
                    );
                    black_box(&out);
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_format_hit);
criterion_main!(benches);
