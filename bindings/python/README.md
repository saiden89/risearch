# risearch Python bindings

## Quick start

From the repo root:

```bash
cd bindings/python
uv sync --locked
uv run --locked maturin develop
uv run --locked python
```

`uv` manages the locked development environment; Maturin builds and installs
the local extension into that environment. Run `maturin develop` again after
changing Rust binding code.

Import the public package as `risearch`. The compiled module is installed as
`risearch._native` and is a private implementation detail.

This Python package does not install the `risearch` command-line program. From
the repository root, install only the CLI with
`cargo install --locked --path .`.

## Legacy CPU / Rosetta install

For older x86-64 CPUs or x86-64 Python running under Rosetta on Apple Silicon,
install the compatibility Polars runtime through the `lts-cpu` extra:

```bash
uv sync --locked --extra lts-cpu
uv run --locked maturin develop
```

For published wheels, the equivalent pip form is:

```bash
pip install "risearch[lts-cpu]"
```

This extra uses Polars' `rtcompat` runtime, so the Python module is still
imported as `polars`.

## What the binding exposes

The Python package is very small:

- `risearch.index(fasta, output)` builds a binary target index
- `risearch.TargetRegistry.open(path)` opens that index for reuse
- `risearch.search(query_fasta, store, **kwargs)` runs the search and returns a Polars `DataFrame`

The result schema is:

| Column | Polars type | Meaning |
| --- | --- | --- |
| `query_idx`, `target_idx` | `UInt64` | Registry positions |
| `query_name`, `target_name` | `String` | FASTA identifiers |
| `q_start`, `q_end`, `t_start`, `t_end` | `UInt64` | Zero-based, inclusive coordinates |
| `strand` | `String` | `+` or `-` |
| `energy` | `Float64` | Free energy in kcal/mol |
| `alignment` | `String` (nullable) | Pairing fingerprint |

The `search()` kwargs map to the canonical Rust-facing options:

- `seed_length`, `seed_start`, `seed_end`
- `mismatches`, `mismatch_prefix`, `mismatch_suffix`
- `seed_wobble`
- `matrix`, `penalty`, `temperature`
- `max_extension`
- `energy_threshold`, `seed_energy`, `no_max_prune`, `no_dedup`

## Smoke test

Run this from `bindings/python/`:

```bash
uv run --locked python - <<'PY'
from pathlib import Path
import tempfile
import risearch

root = Path.cwd().parents[1]
target_fa = root / "tests" / "data" / "target.fa"
query_fa = root / "tests" / "data" / "query.fa"

with tempfile.TemporaryDirectory() as tmp:
    idx = Path(tmp) / "RHOC.idx"
    risearch.index(target_fa, idx)
    store = risearch.TargetRegistry.open(idx)
    df = risearch.search(query_fa, store, seed_length=8, energy_threshold=-10.0)
    print(df.shape)
    print(df.columns)
PY
```

## Development check

Rust-side build check:

```bash
cargo test -p risearch-python
```

Build the extension and run the Python suite:

```bash
uv sync --locked
uv run --locked maturin develop
uv run --locked pytest -q
```

Build a wheel directly:

```bash
uv run --python 3.10 --locked maturin build --out ../../dist
```

Build the source distribution and verify it by rebuilding a wheel from the
unpacked archive:

```bash
uv run --python 3.10 --locked maturin build --sdist --out ../../dist
```
