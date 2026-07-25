import os
from typing import Optional

import polars as pl

class TargetRegistry:
    """Mmap-backed target index for RNA-RNA interaction search."""

    @staticmethod
    def open(path: str | os.PathLike) -> "TargetRegistry":
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
    store: TargetRegistry,
    *,
    seed_length: Optional[int] = None,
    seed_start: Optional[int] = None,
    seed_end: Optional[int] = None,
    energy_threshold: float = -20.0,
    mismatches: int = 0,
    mismatch_prefix: int = 1,
    mismatch_suffix: int = 0,
    seed_wobble: bool = True,
    matrix: str = "t04",
    penalty: float = 3.5,
    temperature: int = 37,
    max_extension: int = 20,
    seed_energy: float = 0.0,
    no_max_prune: bool = False,
    no_dedup: bool = False,
) -> pl.DataFrame:
    """Search for RNA-RNA interactions between queries and an indexed target set.

    Parameters
    ----------
    query_fasta:
        FASTA file (or list of files) containing query sequences.
    store:
        Pre-loaded target index (from ``TargetRegistry.open()``).
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
    seed_wobble:
        Allow G·U wobble pairs in the seed (default: True). Set False for
        strict Watson-Crick seed pairing.
    matrix:
        Energy parameter set: ``"t04"`` (Turner 2004) or ``"t99"`` (Turner 1999).
    penalty:
        Per-nucleotide penalty in kcal/mol.
    temperature:
        Folding temperature in °C (default: 37).
    max_extension:
        Maximum number of bases to extend the seed in each direction via DP.
        Negative means unlimited (extend across the whole query).
    seed_energy:
        Energy-per-length threshold for filtering seeds before extension.
    no_max_prune:
        Disable the maximality check (keeps redundant sub-maximal seeds).
    no_dedup:
        Disable bounding-box deduplication of overlapping hits.

    Returns
    -------
    pl.DataFrame
        One row per hit. Row order is not guaranteed and varies between runs on
        identical input; sort explicitly if you need a stable order. Columns:
        ``query_idx`` (UInt32), ``target_idx`` (UInt32),
        ``q_start`` / ``q_end`` / ``t_start`` / ``t_end`` (UInt32),
        ``strand`` (String, "+"/"-"), ``energy`` (Float64, kcal/mol),
        ``alignment`` (String, nullable — fingerprint e.g. "PPWW...").
    """
    ...
