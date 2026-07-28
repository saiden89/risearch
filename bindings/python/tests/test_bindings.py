"""Smoke tests for the risearch Python bindings.

Requires the extension to be built first:
    cd bindings/python && uv run --locked maturin develop
"""

import inspect
from pathlib import Path

import polars as pl
import pytest
import risearch
import risearch._native as native
from polars.testing import assert_frame_equal

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

DATA = Path(__file__).resolve().parents[3] / "tests" / "data"
TARGET_FA = DATA / "target.fa"
QUERY_FA = DATA / "query.fa"


@pytest.fixture(scope="session")
def target_index(tmp_path_factory):
    idx = tmp_path_factory.mktemp("idx") / "RHOC.idx"
    risearch.index(TARGET_FA, idx)
    return idx


@pytest.fixture(scope="session")
def store(target_index):
    return risearch.TargetRegistry.open(target_index)


# ---------------------------------------------------------------------------
# index / TargetRegistry.open
# ---------------------------------------------------------------------------


def test_index_creates_file(tmp_path):
    idx = tmp_path / "out.idx"
    risearch.index(TARGET_FA, idx)
    assert idx.exists() and idx.stat().st_size > 0


def test_target_registry_repr_shows_count(store):
    assert "targets=" in repr(store)
    count = int(repr(store).split("targets=")[1].rstrip(")"))
    assert count > 0


# ---------------------------------------------------------------------------
# search — DataFrame shape and schema
# ---------------------------------------------------------------------------

EXPECTED_COLUMNS = {
    "query_idx": pl.UInt64,
    "query_name": pl.String,
    "target_idx": pl.UInt64,
    "target_name": pl.String,
    "q_start": pl.UInt64,
    "q_end": pl.UInt64,
    "t_start": pl.UInt64,
    "t_end": pl.UInt64,
    "strand": pl.String,
    "energy": pl.Float64,
    "alignment": pl.String,
}


def test_search_returns_dataframe(store):
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    assert isinstance(df, pl.DataFrame)


def test_search_schema(store):
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    for col, dtype in EXPECTED_COLUMNS.items():
        assert col in df.columns, f"missing column: {col}"
        assert df[col].dtype == dtype, f"{col}: expected {dtype}, got {df[col].dtype}"


def test_search_returns_hits(store):
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    assert len(df) > 0


def test_search_resolves_registry_names(store):
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    assert df["query_name"].is_not_null().all()
    assert df["target_name"].is_not_null().all()
    assert (df["query_name"].str.len_chars() > 0).all()
    assert (df["target_name"].str.len_chars() > 0).all()


def test_empty_result_has_correct_schema(store):
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-1000.0)
    assert isinstance(df, pl.DataFrame)
    assert len(df) == 0
    for col, dtype in EXPECTED_COLUMNS.items():
        assert col in df.columns
        assert df[col].dtype == dtype


# ---------------------------------------------------------------------------
# search — value correctness
# ---------------------------------------------------------------------------


def test_all_energies_below_threshold(store):
    threshold = -15.0
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=threshold)
    assert (df["energy"] <= threshold).all()


def test_strand_values(store):
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    assert df["strand"].is_in(["+", "-"]).all()


def test_coordinates_non_negative(store):
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    for col in ("q_start", "q_end", "t_start", "t_end"):
        assert (df[col] >= 0).all()


def test_start_lte_end(store):
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    assert (df["q_start"] <= df["q_end"]).all()
    assert (df["t_start"] <= df["t_end"]).all()


# ---------------------------------------------------------------------------
# search — filtering behaviour
# ---------------------------------------------------------------------------


def test_stricter_threshold_returns_fewer_hits(store):
    loose = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-5.0)
    strict = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-20.0)
    assert len(strict) <= len(loose)


def test_longer_seed_returns_fewer_hits(store):
    short = risearch.search(QUERY_FA, store, seed_length=6, energy_threshold=-10.0)
    long_ = risearch.search(QUERY_FA, store, seed_length=12, energy_threshold=-10.0)
    assert len(long_) <= len(short)


# ---------------------------------------------------------------------------
# search — multi-file input
# ---------------------------------------------------------------------------


def test_single_path_and_list_equivalent(store):
    df_str = risearch.search(
        str(QUERY_FA), store, seed_length=8, energy_threshold=-10.0
    )
    df_list = risearch.search([QUERY_FA], store, seed_length=8, energy_threshold=-10.0)
    assert len(df_str) == len(df_list)


def test_split_files_match_full_file(store, tmp_path):
    lines = QUERY_FA.read_text().splitlines()
    records, current = [], []
    for line in lines:
        if line.startswith(">") and current:
            records.append("\n".join(current))
            current = []
        current.append(line)
    if current:
        records.append("\n".join(current))

    mid = len(records) // 2
    f1, f2 = tmp_path / "part1.fa", tmp_path / "part2.fa"
    f1.write_text("\n".join(records[:mid]) + "\n")
    f2.write_text("\n".join(records[mid:]) + "\n")

    df_split = risearch.search([f1, f2], store, seed_length=8, energy_threshold=-10.0)
    df_full = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    assert len(df_split) == len(df_full)


# ---------------------------------------------------------------------------
# search — invalid arguments
# ---------------------------------------------------------------------------


def test_seed_start_without_end_raises(store):
    with pytest.raises(ValueError, match="seed_start"):
        risearch.search(QUERY_FA, store, seed_start=2)


def test_strict_seed_returns_fewer_or_equal_hits(store):
    wobble = risearch.search(
        QUERY_FA, store, seed_length=8, energy_threshold=-10.0, seed_wobble=True
    )
    strict = risearch.search(
        QUERY_FA, store, seed_length=8, energy_threshold=-10.0, seed_wobble=False
    )
    assert len(strict) <= len(wobble)


def test_invalid_matrix_raises(store):
    with pytest.raises(ValueError, match="DSM id"):
        risearch.search(QUERY_FA, store, matrix="t05")


def test_excessive_max_extension_raises(store):
    with pytest.raises(ValueError, match="max extension"):
        risearch.search(QUERY_FA, store, max_extension=257)


# ---------------------------------------------------------------------------
# Public/private package boundary
# ---------------------------------------------------------------------------


def test_public_objects_identify_as_risearch():
    assert risearch.TargetRegistry.__module__ == "risearch"
    assert risearch.index.__module__ == "risearch"
    assert risearch.index.__name__ == "index"
    assert risearch.search.__module__ == "risearch"


def test_public_and_native_search_defaults_match():
    public = inspect.signature(risearch.search)
    private = inspect.signature(native.search)
    assert list(public.parameters) == list(private.parameters)
    for name, parameter in public.parameters.items():
        native_default = private.parameters[name].default
        if native_default is not Ellipsis:
            assert parameter.default == native_default
    assert public.parameters["energy_threshold"].default == -20.0
    assert public.parameters["penalty"].default == 0.0


def test_native_negative_energy_default_matches_public_default(store):
    implicit = pl.from_arrow(native.search([QUERY_FA], store, seed_length=8))
    explicit = pl.from_arrow(
        native.search([QUERY_FA], store, seed_length=8, energy_threshold=-20.0)
    )
    columns = list(EXPECTED_COLUMNS)
    assert_frame_equal(implicit.sort(columns), explicit.sort(columns))
