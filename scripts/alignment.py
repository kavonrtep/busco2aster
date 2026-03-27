"""Helpers for retained-locus FASTA export and protein alignment paths."""

from __future__ import annotations

import csv
from pathlib import Path

from .locus_matrix import parse_fasta_records


def busco_sort_key(locus_id: str) -> tuple[int, str]:
    prefix, _, suffix = locus_id.partition("at")
    if prefix.isdigit():
        return int(prefix), suffix
    return 0, locus_id


def locus_output_paths(locus_id: str) -> dict[str, str]:
    raw_fasta = Path("results") / "loci" / "raw_fastas" / f"{locus_id}.faa"
    alignment = Path("results") / "loci" / "alignments" / f"{locus_id}.aln.faa"
    log = Path("results") / "loci" / "logs" / "mafft" / f"{locus_id}.log"
    return {
        "raw_fasta": raw_fasta.as_posix(),
        "alignment": alignment.as_posix(),
        "log": log.as_posix(),
    }


def _resolve_repo_path(repo_root: Path, path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else repo_root / path


def _load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def load_retained_locus_ids(path: Path) -> list[str]:
    rows = _load_tsv_rows(path)
    retained = [row["locus_id"] for row in rows if row.get("decision") == "retain"]
    return sorted(retained, key=busco_sort_key)


def _load_retained_locus_row(path: Path, locus_id: str) -> dict[str, str]:
    matches = [row for row in _load_tsv_rows(path) if row.get("locus_id") == locus_id]
    if not matches:
        raise ValueError(f"Locus {locus_id!r} was not found in retained loci table {path}.")
    if len(matches) > 1:
        raise ValueError(f"Locus {locus_id!r} appears multiple times in retained loci table {path}.")
    row = matches[0]
    if row.get("decision") != "retain":
        raise ValueError(f"Locus {locus_id!r} is not retained and cannot be exported.")
    return row


def _read_single_sequence(faa_path: str, repo_root: Path) -> str:
    resolved_path = _resolve_repo_path(repo_root, faa_path)
    if not resolved_path.is_file():
        raise ValueError(f"Expected locus sequence FASTA does not exist: {resolved_path}")
    records = parse_fasta_records(resolved_path)
    if len(records) != 1:
        raise ValueError(f"Expected exactly one FASTA record in {resolved_path}, found {len(records)}.")
    _, sequence = records[0]
    return sequence.strip().upper()


def export_locus_fasta(
    locus_id: str,
    matrix_path: Path,
    retained_path: Path,
    output_path: Path,
    repo_root: Path,
) -> None:
    retained_row = _load_retained_locus_row(retained_path, locus_id)
    matrix_rows = [
        row
        for row in _load_tsv_rows(matrix_path)
        if row.get("locus_id") == locus_id and row.get("include_in_occupancy") == "true"
    ]
    if not matrix_rows:
        raise ValueError(f"Locus {locus_id!r} has no retained matrix rows in {matrix_path}.")

    matrix_rows.sort(key=lambda row: row["sanitized_taxon_id"])
    headers = [row["sanitized_taxon_id"] for row in matrix_rows]
    if ",".join(headers) != retained_row.get("retained_sanitized_taxon_ids", ""):
        raise ValueError(
            f"Locus {locus_id!r} retained taxon IDs do not match between matrix and retained table."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for row in matrix_rows:
            sequence = _read_single_sequence(row["faa_path"], repo_root)
            handle.write(f">{row['sanitized_taxon_id']}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
