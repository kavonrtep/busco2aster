"""Helpers for retained-locus FASTA export and protein alignment paths."""

from __future__ import annotations

import csv
from pathlib import Path

from .locus_matrix import parse_fasta_records
from .manifest import write_tsv


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
    rows = load_retained_locus_rows(path)
    return [row["locus_id"] for row in rows]


def load_retained_locus_rows(path: Path) -> list[dict[str, str]]:
    rows = _load_tsv_rows(path)
    retained = [row for row in rows if row.get("decision") == "retain"]
    return sorted(retained, key=lambda row: busco_sort_key(row["locus_id"]))


def _build_retained_locus_map(path: Path) -> dict[str, dict[str, str]]:
    retained_rows = load_retained_locus_rows(path)
    retained_map = {row["locus_id"]: row for row in retained_rows}
    if len(retained_map) != len(retained_rows):
        raise ValueError(f"Retained loci table {path} contains duplicate locus IDs.")
    return retained_map


def _load_retained_locus_row(path: Path, locus_id: str) -> dict[str, str]:
    retained_map = _build_retained_locus_map(path)
    row = retained_map.get(locus_id)
    if row is None:
        raise ValueError(f"Locus {locus_id!r} was not found in retained loci table {path}.")
    return row


def _build_matrix_index(
    matrix_path: Path, retained_locus_ids: set[str]
) -> dict[str, list[dict[str, str]]]:
    index = {locus_id: [] for locus_id in retained_locus_ids}
    with matrix_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            locus_id = row.get("locus_id", "")
            if locus_id not in index:
                continue
            if row.get("include_in_occupancy") != "true":
                continue
            index[locus_id].append(dict(row))
    return index


def _read_single_sequence(faa_path: str, repo_root: Path) -> str:
    resolved_path = _resolve_repo_path(repo_root, faa_path)
    if not resolved_path.is_file():
        raise ValueError(f"Expected locus sequence FASTA does not exist: {resolved_path}")
    records = parse_fasta_records(resolved_path)
    if len(records) != 1:
        raise ValueError(f"Expected exactly one FASTA record in {resolved_path}, found {len(records)}.")
    _, sequence = records[0]
    return sequence.strip().upper()


def _build_locus_fasta_records(
    locus_id: str,
    matrix_rows: list[dict[str, str]],
    retained_row: dict[str, str],
    repo_root: Path,
) -> list[tuple[str, str]]:
    if not matrix_rows:
        raise ValueError(f"Locus {locus_id!r} has no retained matrix rows.")

    matrix_rows.sort(key=lambda row: row["sanitized_taxon_id"])
    headers = [row["sanitized_taxon_id"] for row in matrix_rows]
    if ",".join(headers) != retained_row.get("retained_sanitized_taxon_ids", ""):
        raise ValueError(
            f"Locus {locus_id!r} retained taxon IDs do not match between matrix and retained table."
        )

    return [
        (row["sanitized_taxon_id"], _read_single_sequence(row["faa_path"], repo_root))
        for row in matrix_rows
    ]


def _write_fasta_records(output_path: Path, records: list[tuple[str, str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def export_locus_fasta(
    locus_id: str,
    matrix_path: Path,
    retained_path: Path,
    output_path: Path,
    repo_root: Path,
) -> None:
    retained_row = _load_retained_locus_row(retained_path, locus_id)
    matrix_index = _build_matrix_index(matrix_path, {locus_id})
    matrix_rows = matrix_index.get(locus_id, [])
    records = _build_locus_fasta_records(locus_id, matrix_rows, retained_row, repo_root)
    _write_fasta_records(output_path, records)


def export_retained_fastas(
    matrix_path: Path,
    retained_path: Path,
    output_dir: Path,
    manifest_path: Path,
    repo_root: Path,
) -> None:
    retained_map = _build_retained_locus_map(retained_path)
    retained_locus_ids = set(retained_map)
    matrix_index = _build_matrix_index(matrix_path, retained_locus_ids)

    output_dir.mkdir(parents=True, exist_ok=True)
    expected_files = set()
    manifest_rows = []
    for locus_id in sorted(retained_locus_ids, key=busco_sort_key):
        output_path = output_dir / f"{locus_id}.faa"
        records = _build_locus_fasta_records(
            locus_id=locus_id,
            matrix_rows=matrix_index.get(locus_id, []),
            retained_row=retained_map[locus_id],
            repo_root=repo_root,
        )
        _write_fasta_records(output_path, records)
        expected_files.add(output_path.name)
        manifest_rows.append(
            {
                "locus_id": locus_id,
                "raw_fasta": output_path.as_posix(),
                "taxon_count": str(len(records)),
            }
        )

    for path in output_dir.glob("*.faa"):
        if path.name not in expected_files:
            path.unlink()

    write_tsv(
        manifest_path,
        manifest_rows,
        ["locus_id", "raw_fasta", "taxon_count"],
    )
