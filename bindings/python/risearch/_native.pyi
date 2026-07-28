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
) -> None: ...


def search(
    query_fasta: list[str | PathLike[str]],
    store: TargetRegistry,
    **kwargs: object,
) -> SearchResult: ...
