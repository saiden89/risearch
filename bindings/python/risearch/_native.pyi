from os import PathLike


class TargetRegistry:
    """Mmap-backed target index for RNA-RNA interaction search."""

    @staticmethod
    def open(path: str | PathLike[str]) -> TargetRegistry: ...

    def __repr__(self) -> str: ...


class SearchResult:
    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object: ...


def build_index(
    fasta: str | PathLike[str],
    output: str | PathLike[str],
    *,
    threads: int | None = None,
) -> None: ...


def search(
    query: list[str | PathLike[str]],
    target: TargetRegistry,
    *,
    seed_length: int | None,
    seed_start: int | None,
    seed_end: int | None,
    energy_threshold: float,
    mismatches: int,
    mismatch_prefix: int,
    mismatch_suffix: int,
    seed_wobble: bool,
    matrix: str,
    penalty: float,
    temperature: int,
    max_extension: int,
    seed_energy: float,
    no_max_prune: bool,
    no_dedup: bool,
    alignment: bool,
    threads: int | None,
) -> SearchResult: ...


def _default_options() -> dict[str, object]: ...
