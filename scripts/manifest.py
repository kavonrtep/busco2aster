"""Helpers for sample manifest normalization and validation."""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Iterable


REQUIRED_SAMPLE_COLUMNS = ("sample_id", "taxon_id", "assembly_fasta")
REFERENCE_SAMPLE_COLUMNS = ("name", "fasta")


class ManifestValidationError(ValueError):
    """Raised when the sample manifest is invalid."""


def sanitize_taxon_label(value: str) -> str:
    sanitized = re.sub(r"[^a-z0-9]+", "_", value.strip().lower())
    sanitized = re.sub(r"_+", "_", sanitized).strip("_")
    if not sanitized:
        raise ManifestValidationError("Taxon label becomes empty after sanitization.")
    return sanitized


def _read_delimited_rows(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fieldnames = tuple(reader.fieldnames or ())
        rows = [dict(row) for row in reader]
    return fieldnames, rows


def load_samples_manifest(path: Path) -> list[dict[str, str]]:
    fieldnames, rows = _read_delimited_rows(path, delimiter="\t")
    missing = [column for column in REQUIRED_SAMPLE_COLUMNS if column not in fieldnames]
    if missing:
        raise ManifestValidationError(
            f"Sample manifest is missing required columns: {', '.join(missing)}"
        )
    unexpected = [column for column in fieldnames if column not in REQUIRED_SAMPLE_COLUMNS]
    if unexpected:
        raise ManifestValidationError(
            f"Sample manifest has unexpected columns: {', '.join(unexpected)}"
        )
    return rows


def load_reference_manifest(path: Path) -> list[dict[str, str]]:
    fieldnames, rows = _read_delimited_rows(path, delimiter="\t")
    missing = [column for column in REFERENCE_SAMPLE_COLUMNS if column not in fieldnames]
    if missing:
        raise ManifestValidationError(
            f"Reference manifest is missing required columns: {', '.join(missing)}"
        )
    return rows


def normalize_reference_manifest_rows(rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    normalized = []
    for row in rows:
        taxon_id = (row.get("name") or "").strip()
        fasta = (row.get("fasta") or "").strip()
        if not taxon_id or not fasta:
            raise ManifestValidationError("Reference manifest contains an empty name or fasta value.")
        normalized.append(
            {
                "sample_id": sanitize_taxon_label(taxon_id),
                "taxon_id": taxon_id,
                "assembly_fasta": Path("test_data", fasta.removeprefix("./")).as_posix(),
            }
        )
    return normalized


def _normalize_manifest_path(path_value: str) -> str:
    path_text = path_value.strip()
    if not path_text:
        raise ManifestValidationError("Assembly path cannot be empty.")
    path_obj = Path(path_text)
    if path_obj.is_absolute():
        return path_obj.as_posix()
    return Path(path_text.removeprefix("./")).as_posix()


def validate_manifest_rows(
    rows: Iterable[dict[str, str]],
    repo_root: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    validated_rows: list[dict[str, str]] = []
    taxon_rows: list[dict[str, str]] = []
    seen_sample_ids: set[str] = set()
    seen_taxon_ids: dict[str, str] = {}

    for row in rows:
        sample_id = (row.get("sample_id") or "").strip()
        taxon_id = (row.get("taxon_id") or "").strip()
        assembly_fasta = _normalize_manifest_path(row.get("assembly_fasta") or "")

        if not sample_id:
            raise ManifestValidationError("Sample manifest contains an empty sample_id.")
        if not taxon_id:
            raise ManifestValidationError("Sample manifest contains an empty taxon_id.")
        if sample_id in seen_sample_ids:
            raise ManifestValidationError(f"Duplicate sample_id found: {sample_id}")
        seen_sample_ids.add(sample_id)

        sanitized_taxon_id = sanitize_taxon_label(taxon_id)
        previous_taxon = seen_taxon_ids.get(sanitized_taxon_id)
        if previous_taxon:
            raise ManifestValidationError(
                "Duplicate taxon identifier after sanitization: "
                f"{previous_taxon!r} and {taxon_id!r} both map to {sanitized_taxon_id!r}"
            )
        seen_taxon_ids[sanitized_taxon_id] = taxon_id

        assembly_path = Path(assembly_fasta)
        resolved_path = assembly_path if assembly_path.is_absolute() else repo_root / assembly_path
        if not resolved_path.is_file():
            raise ManifestValidationError(
                f"Assembly file does not exist for sample {sample_id}: {assembly_fasta}"
            )

        validated_rows.append(
            {
                "sample_id": sample_id,
                "taxon_id": taxon_id,
                "sanitized_taxon_id": sanitized_taxon_id,
                "assembly_fasta": assembly_fasta,
            }
        )
        taxon_rows.append(
            {
                "taxon_id": taxon_id,
                "sanitized_taxon_id": sanitized_taxon_id,
            }
        )

    taxon_rows.sort(key=lambda row: row["sanitized_taxon_id"])
    return validated_rows, taxon_rows


def write_tsv(path: Path, rows: Iterable[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
