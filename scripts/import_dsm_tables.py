#!/usr/bin/env python3
"""Import selected RIsearch3 DSM tables into the canonical data/dsm catalog."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
LEGACY_ROOT = ROOT / "legacy_c" / "risearch3_src" / "tables"
DATA_ROOT = ROOT / "data" / "dsm"

CANONICAL_BASES = ("-", "A", "C", "G", "N", "U")
CANONICAL_BASE_IDX = {base: idx for idx, base in enumerate(CANONICAL_BASES)}
INVALID_TRANSITION = "20.0000"
RISEARCH2_BASES = ("-", "A", "G", "C", "U", "N")
CANONICAL_TO_RISEARCH2 = tuple(RISEARCH2_BASES.index(base) for base in CANONICAL_BASES)
TEMPERATURES = (
    ("273.15", "0"),
    ("298.15", "25"),
    ("310.15", "37"),
    ("315.15", "42"),
    ("323.15", "50"),
)

# RIsearch3 passes matversion=2 in seed-and-extend mode. This is the exact
# index transform from legacy_c/risearch3_src/src/dsm.c::to_version.
RI3_TO_CANONICAL = (1, 3, 2, 4, 5, 0)


@dataclass(frozen=True)
class MatrixSpec:
    matrix_id: str
    family: str
    query: str
    target: str
    publication: str
    doi: str
    legacy_qt: str
    legacy_name: str
    output_dir: str
    legacy_c_table: str | None = None
    # None = use the global TEMPERATURES list; otherwise only these (legacy_path, temperature) pairs
    temperatures: tuple[tuple[str, str], ...] | None = None


MATRICES = (
    MatrixSpec(
        matrix_id="t04",
        family="turner-2004",
        query="rna",
        target="rna",
        publication="Mathews et al. 2004",
        doi="10.1073/pnas.0401799101",
        legacy_qt="RNA/RNA",
        legacy_name="t04.v4",
        output_dir="t04",
    ),
    MatrixSpec(
        matrix_id="slh04",
        family="santalucia-hicks-2004",
        query="dna",
        target="dna",
        publication="SantaLucia and Hicks 2004",
        doi="10.1146/annurev.biophys.32.110601.141800",
        legacy_qt="DNA/DNA",
        legacy_name="slh04.v4",
        output_dir="slh04",
    ),
    MatrixSpec(
        matrix_id="s95-rna-dna",
        family="sugimoto-1995",
        query="rna",
        target="dna",
        publication="Sugimoto et al. 1995",
        doi="10.1021/bi00035a029",
        legacy_qt="RNA/DNA",
        legacy_name="su95.v4",
        output_dir="s95",
    ),
    MatrixSpec(
        matrix_id="s95-dna-rna",
        family="sugimoto-1995",
        query="dna",
        target="rna",
        publication="Sugimoto et al. 1995",
        doi="10.1021/bi00035a029",
        legacy_qt="RNA/DNA",
        legacy_name="su95.v4",
        output_dir="s95",
    ),
    MatrixSpec(
        matrix_id="t99",
        family="turner-1999",
        query="rna",
        target="rna",
        publication="Mathews et al. 1999",
        doi="10.1006/jmbi.1999.2700",
        legacy_qt="RNA/RNA",
        legacy_name="t99.v2",
        output_dir="t99",
        legacy_c_table="dsm-t99-guugle.c",
        temperatures=(("310.15", "37"),),
    ),
)


def read_legacy_tsv(path: Path) -> tuple[str, dict[tuple[int, int, int, int], str]]:
    tokens = path.read_text().split()
    if len(tokens) != 1297:
        raise ValueError(f"{path}: expected 1297 values, found {len(tokens)}")

    initiation = tokens[0]
    table: dict[tuple[int, int, int, int], str] = {}
    for flat, value in enumerate(tokens[1:]):
        t2_raw = flat % 6
        t1_raw = (flat // 6) % 6
        q2_raw = (flat // 36) % 6
        q1_raw = (flat // 216) % 6
        q1 = CANONICAL_BASE_IDX[RISEARCH2_BASES[RI3_TO_CANONICAL[q1_raw]]]
        q2 = CANONICAL_BASE_IDX[RISEARCH2_BASES[RI3_TO_CANONICAL[q2_raw]]]
        t1 = CANONICAL_BASE_IDX[RISEARCH2_BASES[RI3_TO_CANONICAL[t1_raw]]]
        t2 = CANONICAL_BASE_IDX[RISEARCH2_BASES[RI3_TO_CANONICAL[t2_raw]]]
        key = (q1, q2, t1, t2)
        if key in table:
            raise ValueError(f"{path}: duplicate canonical DSM coordinate {key}")
        table[key] = value

    if len(table) != 1296:
        raise ValueError(f"{path}: expected 1296 canonical DSM entries, found {len(table)}")
    return initiation, table


def read_risearch2_c_table(path: Path) -> dict[tuple[int, int, int, int], str]:
    text = path.read_text()
    start = text.index("= {")
    end = text.index("};", start)
    values = [int(value) for value in re.findall(r"-?\d+", text[start:end])]
    if len(values) != 1296:
        raise ValueError(f"{path}: expected 1296 DSM integers, found {len(values)}")

    table: dict[tuple[int, int, int, int], str] = {}
    for q1 in range(6):
        for q2 in range(6):
            for t1 in range(6):
                for t2 in range(6):
                    c_q1 = CANONICAL_TO_RISEARCH2[q1]
                    c_q2 = CANONICAL_TO_RISEARCH2[q2]
                    c_t1 = CANONICAL_TO_RISEARCH2[t1]
                    c_t2 = CANONICAL_TO_RISEARCH2[t2]
                    raw = values[c_q1 * 216 + c_q2 * 36 + c_t1 * 6 + c_t2]
                    table[(q1, q2, t1, t2)] = f"{-raw / 100:.4f}"
    return table


def pretty_float(raw: str) -> str:
    text = f"{float(raw):.4f}".rstrip("0").rstrip(".")
    return text or "0"


def write_canonical_tsv(path: Path, table: dict[tuple[int, int, int, int], str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("q1\tq2\tt1\tt2\tdelta_g_kcal_per_mol\n")
        for q1 in range(6):
            for q2 in range(6):
                for t1 in range(6):
                    for t2 in range(6):
                        value = table[(q1, q2, t1, t2)]
                        if pretty_float(value) == pretty_float(INVALID_TRANSITION):
                            continue
                        out.write(
                            f"{CANONICAL_BASES[q1]}\t{CANONICAL_BASES[q2]}\t{CANONICAL_BASES[t1]}\t{CANONICAL_BASES[t2]}\t{value}\n"
                        )


def reverse_swap(key: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    q1, q2, t1, t2 = key
    return t2, t1, q2, q1


def validate_mixed_orientation(legacy_path: str, rna_dna: dict[tuple[int, int, int, int], str]) -> None:
    src = LEGACY_ROOT / "DNA/RNA" / legacy_path / "su95.v4.tsv"
    _, dna_rna = read_legacy_tsv(src)
    for key, value in rna_dna.items():
        if value != dna_rna[reverse_swap(key)]:
            raise ValueError(
                f"{src}: DNA/RNA table is not reverse-swap equivalent to RNA/DNA at {key}"
            )


def manifest_file_for(spec: MatrixSpec, file_name: str) -> str:
    if spec.query == "dna" and spec.target == "rna" and spec.matrix_id == "s95-dna-rna":
        return f"s95/{file_name}"
    return f"{spec.output_dir}/{file_name}"


def manifest_orientation(spec: MatrixSpec) -> str:
    if spec.query == "dna" and spec.target == "rna" and spec.matrix_id == "s95-dna-rna":
        return "reverse-swap"
    return "identity"


def write_manifest(entries: list[tuple[MatrixSpec, list[tuple[str, str, str]]]]) -> None:
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    with (DATA_ROOT / "manifest.toml").open("w") as out:
        for spec_idx, (spec, temps) in enumerate(entries):
            if spec_idx:
                out.write("\n")
            out.write("[[dsm]]\n")
            out.write(f'id = "{spec.matrix_id}"\n')
            out.write(f'family = "{spec.family}"\n')
            out.write(f'query = "{spec.query}"\n')
            out.write(f'target = "{spec.target}"\n')
            out.write(f'publication = "{spec.publication}"\n')
            out.write(f'doi = "{spec.doi}"\n')
            out.write(f"invalid_transition_kcal = {pretty_float(INVALID_TRANSITION)}\n")
            out.write(f'orientation = "{manifest_orientation(spec)}"\n')
            out.write("default_temperature = 37\n")
            out.write("temperatures = [\n")
            for temperature, file_name, initiation in temps:
                out.write(
                    f'  {{ temperature = {temperature}, file = "{manifest_file_for(spec, file_name)}", '
                    f"initiation_kcal = {initiation} }},\n"
                )
            out.write("]\n")


def main() -> None:
    manifest_entries: list[tuple[MatrixSpec, list[tuple[str, str, str]]]] = []
    for spec in MATRICES:
        temp_entries: list[tuple[str, str, str]] = []
        for legacy_path, temperature in (spec.temperatures or TEMPERATURES):
            if spec.legacy_c_table is None:
                src = LEGACY_ROOT / spec.legacy_qt / legacy_path / f"{spec.legacy_name}.tsv"
                initiation, table = read_legacy_tsv(src)
            else:
                initiation = "5.59"
                table = read_risearch2_c_table(
                    ROOT / "legacy_c" / "RIsearch2" / "src" / spec.legacy_c_table
                )
            file_name = f"{temperature}.tsv"
            if manifest_orientation(spec) == "identity":
                dst = DATA_ROOT / spec.output_dir / file_name
                write_canonical_tsv(dst, table)
                if spec.query == "rna" and spec.target == "dna" and spec.matrix_id == "s95-rna-dna":
                    validate_mixed_orientation(legacy_path, table)
            temp_entries.append((temperature, file_name, pretty_float(initiation)))
        manifest_entries.append((spec, temp_entries))
    write_manifest(manifest_entries)


if __name__ == "__main__":
    main()
