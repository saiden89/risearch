# RIsearch (Rust)

Energy-based RNA-RNA interaction prediction with a fast suffix-array backend.

Status: alpha. The Rust port is actively optimized and validated against the
legacy C implementation. Some flags and features are still marked TODO or
experimental.

## Features

- Suffix-array (SA) index backend (default).
- Seed-and-extend search with configurable seed length, mismatches, and energy
  thresholds.
- Multi-threaded search via Rayon.
- Streaming output to avoid large memory spikes.
- FASTA/FASTQ input (including .gz) via needletail.

## Build

```bash
cargo build --release
```

Optional (work-in-progress) FM-index backend:

```bash
cargo build --release --features fm-index
```

Install the CLI locally:

```bash
cargo install --path .
```

## Quick start

Index a target FASTA:

```bash
./target/release/risearch index chr22.fa chr22.idx
```

Search for interactions:

```bash
./target/release/risearch search \
  -q queries.fa \
  -i chr22.idx \
  -o results.out
```

Or run directly via Cargo:

```bash
cargo run --release -- search \
  -q queries.fa \
  -i chr22.idx \
  -o results.out
```

Quick sanity check (no data needed):

```bash
cargo run --release -- --help
```

Common tuning flags:

```bash
./target/release/risearch search \
  -q queries.fa \
  -i chr22.idx \
  -o results.out \
  -s 6            \
  -l 20           \
  -e -20          \
  -m 1:3          \
  -f=cigar        \
  -t 8
```

## CLI overview

Global options:

- `-t, --threads <N>`: number of worker threads.
- `-b, --backend <sa|fm>`: index backend (SA default; FM is WIP).
- `-v/-vv/-vvv`: increase logging verbosity.

Subcommands:

- `index <INPUT> <OUTPUT>`: build an index from a FASTA/FASTQ file.
- `search -q <QUERY> -i <INDEX> -o <OUTPUT>`: run the search pipeline.

Search options (selected):

- `-s, --seed <len|start:end>`: seed length or length range (default 6).
- `-m, --mismatch <c:p>`: allow up to `c` mismatches with `p` consecutive
  matches at seed ends.
- `-l, --extension <L>`: max extension length on each side (default 20).
- `-e, --energy <dG>`: energy threshold in kcal/mol (default -20.0).
- `-p, --penalty <P>`: per-nucleotide extension penalty.
- `-f, --format <detailed|cigar|bindingsite>`: output format.

Run `risearch --help` or `risearch search --help` for the full list.

## Input normalization

Sequences are normalized during indexing/search:

- Lowercased and kept as RNA bases.
- `T`/`U` are normalized to `T`.
- Ambiguous bases map to `N`.
- Gaps (`-`/`.`) are removed.
- Duplicate FASTA IDs are rejected.

## Tests and benches

Run the test suite:

```bash
cargo test
```

Benchmarks and profiling utilities live in `bench.sh`, `bench_band.sh`, and
`docs/`. These are evolving; treat them as developer tooling.

## Docs

Additional technical notes and design docs are in `docs/`:

- `docs/implementation_plan.md`
- `docs/strategies.md`
- `docs/parity_divergences.md`
- `docs/SIMD_ARCHITECTURE.md`

## License

See `LICENSE`.
