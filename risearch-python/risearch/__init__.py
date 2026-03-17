import os

import polars as pl

from risearch._risearch import (
    PyTargetStore as TargetStore,
    build_index as index,
    search as _search_impl,
)


def search(query_fasta, store, **kwargs) -> pl.DataFrame:
    if isinstance(query_fasta, (str, os.PathLike)):
        query_fasta = [query_fasta]
    result = _search_impl(query_fasta, store, **kwargs)
    return pl.from_arrow(result)


__all__ = ["TargetStore", "index", "search"]
