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
# search — alignment and threads
# ---------------------------------------------------------------------------

SORT_COLS = [c for c in EXPECTED_COLUMNS if c != "alignment"]


def test_alignment_populated_by_default(store):
    df = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    assert df["alignment"].is_not_null().any()


def test_alignment_false_nulls_column_but_keeps_schema(store):
    df = risearch.search(
        QUERY_FA, store, seed_length=8, energy_threshold=-10.0, alignment=False
    )
    assert df["alignment"].is_null().all()
    for col, dtype in EXPECTED_COLUMNS.items():
        assert col in df.columns
        assert df[col].dtype == dtype


def test_alignment_false_preserves_hits(store):
    with_aln = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    without = risearch.search(
        QUERY_FA, store, seed_length=8, energy_threshold=-10.0, alignment=False
    )
    assert_frame_equal(
        with_aln.drop("alignment").sort(SORT_COLS),
        without.drop("alignment").sort(SORT_COLS),
    )


@pytest.mark.parametrize("threads", [1, 2])
def test_threads_does_not_change_results(store, threads):
    default = risearch.search(QUERY_FA, store, seed_length=8, energy_threshold=-10.0)
    pinned = risearch.search(
        QUERY_FA, store, seed_length=8, energy_threshold=-10.0, threads=threads
    )
    columns = list(EXPECTED_COLUMNS)
    assert_frame_equal(default.sort(columns), pinned.sort(columns))


def test_index_accepts_threads(tmp_path):
    idx = tmp_path / "threaded.idx"
    risearch.index(TARGET_FA, idx, threads=2)
    assert idx.exists() and idx.stat().st_size > 0


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
    with pytest.raises(risearch.ModelError, match="DSM id"):
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


def test_public_and_native_search_params_match():
    """The wrapper must forward every native keyword, in the same order."""
    public = inspect.signature(risearch.search)
    private = inspect.signature(native.search)
    assert list(public.parameters) == list(private.parameters)


def test_native_search_declares_no_defaults():
    """Defaults belong to the wrapper alone, so the native side must have none."""
    private = inspect.signature(native.search)
    defaulted = [
        name
        for name, p in private.parameters.items()
        if p.default is not inspect.Parameter.empty
    ]
    assert defaulted == []


def test_wrapper_defaults_match_rust_config_defaults():
    """The wrapper's defaults are the Rust config's defaults, not a copy of them."""
    public = inspect.signature(risearch.search).parameters
    exempt = {"query", "target", "threads", "alignment"}
    assert {n: p.default for n, p in public.items() if n not in exempt} == (
        native._default_options()
    )


def test_documented_defaults_are_pinned():
    """Absolute pins: the relative check above moves with the Rust config."""
    public = inspect.signature(risearch.search).parameters
    assert public["energy_threshold"].default == -20.0
    assert public["penalty"].default == 0.0
    assert public["seed_wobble"].default is False
    assert public["max_extension"].default == 20
    assert public["temperature"].default == 37
    assert public["matrix"].default == "t04"


def test_native_search_is_callable_with_every_kwarg(store):
    """Keeps a runtime call on the native entry point the stub documents."""
    result = native.search(
        [QUERY_FA],
        store,
        **{**native._default_options(), "seed_length": 8, "alignment": True},
        threads=None,
    )
    assert pl.DataFrame(result).height > 0


# ---------------------------------------------------------------------------
# Error mapping
# ---------------------------------------------------------------------------


def test_risearch_exceptions_share_a_base():
    """One `except risearch.RisearchError` catches every risearch-specific error."""
    for exc in (
        risearch.IndexFormatError,
        risearch.InputError,
        risearch.ModelError,
        risearch.SearchError,
    ):
        assert issubclass(exc, risearch.RisearchError)
        assert issubclass(exc, Exception)


def test_missing_index_raises_file_not_found(tmp_path):
    """An absent index is an OS-level miss, not a risearch-specific failure."""
    with pytest.raises(FileNotFoundError):
        risearch.TargetRegistry.open(tmp_path / "absent.idx")


def test_corrupt_index_raises_index_format_error(tmp_path):
    """A file that is not an index is rejected by the header, before the archive."""
    bogus = tmp_path / "bogus.idx"
    bogus.write_bytes(b"definitely not a risearch index")
    with pytest.raises(risearch.IndexFormatError):
        risearch.TargetRegistry.open(bogus)


def test_missing_query_file_raises_file_not_found(store):
    """A query path that does not exist surfaces as FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        risearch.search("no-such-query.fa", store)
