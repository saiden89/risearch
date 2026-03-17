import os
from typing import Optional

import polars as pl

class TargetStore:
    """In-memory index of target sequences for RNA-RNA interaction search."""

    @staticmethod
    def open(path: str | os.PathLike) -> "TargetStore":
        """Load a pre-built index from disk."""
        ...

    def __repr__(self) -> str: ...

def index(
    fasta: str | os.PathLike,
    output: str | os.PathLike,
) -> None:
    """Build a binary index from a FASTA file of target sequences.

    Parameters
    ----------
    fasta:
        Input FASTA file containing target sequences.
    output:
        Destination path for the binary index.
    """
    ...

def search(
    query_fasta: str | os.PathLike | list[str | os.PathLike],
    store: TargetStore,
    *,
    seed_length: Optional[int] = None,
    seed_start: Optional[int] = None,
    seed_end: Optional[int] = None,
    energy_threshold: float = -20.0,
    mismatches: int = 0,
    mismatch_prefix: int = 1,
    mismatch_suffix: int = 0,
    seed_pairing: str = "allow_wobble",
    matrix: str = "t04",
    penalty: float = 3.5,
    max_extension: int = 20,
    seed_energy: float = 0.0,
    no_max_prune: bool = False,
) -> pl.DataFrame:
    """Return a Polars DataFrame with one row per hit.

    Columns
    -------
    query_idx : UInt32
    target_idx : UInt32
    q_start : UInt32
    q_end : UInt32
    t_start : UInt32
    t_end : UInt32
    strand : String  — "+" or "-"
    energy : Float64 — kcal/mol
    alignment : String (nullable) — fingerprint e.g. "PPWW..."
    """
    """Search for RNA-RNA interactions between queries and an indexed target set.

    Parameters
    ----------
    query_fasta:
        FASTA file containing query sequences.
    store:
        Pre-loaded target index (from ``TargetStore.open()``).
    seed_length:
        Seed length (default: 6). Use alone for length-only mode, or combine
        with ``seed_start``/``seed_end`` to fix the length within an interval.
    seed_start:
        Seed interval start position (1-based, negative counts from 3' end).
        Must be provided together with ``seed_end``.
    seed_end:
        Seed interval end position (1-based, negative counts from 3' end).
        Must be provided together with ``seed_start``.
    energy_threshold:
        Maximum free energy (kcal/mol) for a hit to be reported.
    mismatches:
        Maximum number of mismatches allowed in the seed.
    mismatch_prefix:
        Minimum consecutive matches required at the seed 5' end.
    mismatch_suffix:
        Minimum consecutive matches required at the seed 3' end.
    seed_pairing:
        Seed pairing mode: ``"wobble"`` (allows G-U pairs) or ``"strict"``.
    matrix:
        Energy parameter set: ``"t04"`` (Turner 2004) or ``"t99"`` (Turner 1999).
    penalty:
        Per-nucleotide penalty in kcal/mol.
    max_extension:
        Maximum number of bases to extend the seed in each direction via DP.
    seed_energy:
        Energy-per-length threshold for filtering seeds before extension.
    no_max_prune:
        Disable the maximality check (keeps redundant sub-maximal seeds).

    Returns
    -------
    list[SearchHit]
        All hits passing the energy threshold. Order is not guaranteed.
    """
    ...
