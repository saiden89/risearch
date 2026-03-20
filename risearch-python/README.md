# risearch Python bindings

## Quick start

From the repo root:

```bash
cd risearch-python
uv sync
uv run python
```

`uv run` will build/install the local `risearch` package into the project venv as needed.

## What the binding exposes

The Python package is very small:

- `risearch.index(fasta, output)` builds a binary target index
- `risearch.TargetStore.open(path)` opens that index for reuse
- `risearch.search(query_fasta, store, **kwargs)` runs the search and returns a Polars `DataFrame`

The `search()` kwargs map to the canonical Rust-facing options:

- `seed_length`, `seed_start`, `seed_end`
- `mismatches`, `mismatch_prefix`, `mismatch_suffix`
- `seed_pairing`
- `matrix`, `penalty`
- `max_extension`
- `energy_threshold`, `seed_energy`, `no_max_prune`

## Smoke test

Run this from `risearch-python/`:

```bash
uv run python - <<'PY'
from pathlib import Path
import tempfile
import risearch

root = Path.cwd().parent
suite = root / "legacy_c" / "RIsearch2" / "test_suite"
target_fa = suite / "RHOC.fa"
query_fa = suite / "mirnas.fa"

with tempfile.TemporaryDirectory() as tmp:
    idx = Path(tmp) / "RHOC.idx"
    risearch.index(target_fa, idx)
    store = risearch.TargetStore.open(idx)
    df = risearch.search(query_fa, store, seed_length=8, energy_threshold=-10.0)
    print(df.shape)
    print(df.columns)
PY
```

Expected shape today is `(81, 9)`.

## Development check

Rust-side build check:

```bash
cargo test -p risearch-python
```

There is also a Python smoke suite in `tests/test_bindings.py`, but `pytest` is not currently declared in `pyproject.toml`, so it will not run until `pytest` is installed in the environment.
