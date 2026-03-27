"""Helpers for BUSCO-derived locus matrix construction and selection."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Iterable

from .manifest import write_tsv


AMINO_ACID_ALPHABET = set("ABCDEFGHIKLMNPQRSTVWXYZJUO")
MATRIX_FIELDNAMES = [
    "locus_id",
    "sample_id",
    "taxon_id",
    "sanitized_taxon_id",
    "assembly_fasta",
    "status",
    "record_count",
    "sequence_count",
    "sequence_count_mismatch",
    "full_table_sequence_ids",
    "gene_starts",
    "gene_ends",
    "strands",
    "scores",
    "match_lengths_aa",
    "faa_path",
    "gff_path",
    "protein_lengths_aa",
    "representative_protein_length_aa",
    "has_terminal_stop_codon",
    "terminal_stop_sequence_count",
    "has_internal_stop_codon",
    "has_invalid_amino_acid",
    "invalid_amino_acids",
    "is_complete_single_copy",
    "include_in_occupancy",
    "qc_notes",
    "locus_total_taxa",
    "locus_complete_single_copy_taxa",
    "locus_duplicated_taxa",
    "locus_fragmented_taxa",
    "locus_missing_taxa",
    "locus_stop_codon_taxa",
    "locus_invalid_amino_acid_taxa",
    "locus_occupancy",
]
RETAINED_FIELDNAMES = [
    "locus_id",
    "total_taxa",
    "complete_single_copy_taxa",
    "duplicated_taxa",
    "fragmented_taxa",
    "missing_taxa",
    "stop_codon_taxa",
    "invalid_amino_acid_taxa",
    "occupancy",
    "occupancy_threshold",
    "passes_occupancy",
    "decision",
    "failure_reasons",
    "qc_warnings",
    "retained_sample_ids",
    "retained_sanitized_taxon_ids",
    "retained_protein_lengths_aa",
    "min_protein_length_aa",
    "median_protein_length_aa",
    "max_protein_length_aa",
    "length_span_aa",
    "length_ratio",
    "length_dispersion_observed",
]


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def load_busco_summary_rows(path: Path) -> list[dict[str, str]]:
    rows = _read_tsv_rows(path)
    required = {"sample_id", "taxon_id", "sanitized_taxon_id", "assembly_fasta"}
    missing = required.difference(rows[0].keys() if rows else set())
    if missing:
        raise ValueError(f"BUSCO summary is missing required columns: {', '.join(sorted(missing))}")
    return rows


def load_busco_record_rows(path: Path) -> list[dict[str, str]]:
    rows = _read_tsv_rows(path)
    required = {"sample_id", "busco_id", "status", "faa_path", "gff_path"}
    missing = required.difference(rows[0].keys() if rows else set())
    if missing:
        raise ValueError(f"BUSCO record table is missing required columns: {', '.join(sorted(missing))}")
    return rows


def _resolve_repo_path(repo_root: Path, path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else repo_root / path


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


def _as_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def _format_occupancy(value: float) -> str:
    return f"{value:.4f}"


def _format_stat(value: float | int | str) -> str:
    if value == "":
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    if value.is_integer():
        return str(int(value))
    return f"{value:.4f}"


def _join(values: Iterable[str]) -> str:
    return ",".join(value for value in values if value)


def _busco_sort_key(busco_id: str) -> tuple[int, str]:
    prefix, _, suffix = busco_id.partition("at")
    if prefix.isdigit():
        return int(prefix), suffix
    return 0, busco_id


def parse_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence_lines: list[str] = []

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                records.append((header, "".join(sequence_lines)))
            header = line[1:]
            sequence_lines = []
            continue
        sequence_lines.append(line)

    if header:
        records.append((header, "".join(sequence_lines)))
    return records


def analyze_protein_sequences(faa_path: str, repo_root: Path) -> dict[str, object]:
    if not faa_path:
        return {
            "sequence_count": 0,
            "protein_lengths_aa": "",
            "representative_protein_length_aa": "",
            "has_terminal_stop_codon": "false",
            "terminal_stop_sequence_count": 0,
            "has_internal_stop_codon": "false",
            "has_invalid_amino_acid": "false",
            "invalid_amino_acids": "",
        }

    resolved_path = _resolve_repo_path(repo_root, faa_path)
    if not resolved_path.is_file():
        raise ValueError(f"Expected FASTA file does not exist: {resolved_path}")

    records = parse_fasta_records(resolved_path)
    if not records:
        raise ValueError(f"Expected at least one FASTA record in {resolved_path}")

    lengths: list[int] = []
    invalid_chars: set[str] = set()
    terminal_stop_sequence_count = 0
    has_internal_stop_codon = False

    for _, sequence in records:
        normalized = sequence.strip().upper()
        if normalized.endswith("*"):
            terminal_stop_sequence_count += 1
            core = normalized[:-1]
        else:
            core = normalized
        if "*" in core:
            has_internal_stop_codon = True
        invalid_chars.update(char for char in core if char not in AMINO_ACID_ALPHABET and char != "*")
        lengths.append(len(normalized.replace("*", "")))

    return {
        "sequence_count": len(records),
        "protein_lengths_aa": _join(str(length) for length in lengths),
        "representative_protein_length_aa": str(lengths[0]) if len(lengths) == 1 else "",
        "has_terminal_stop_codon": _bool_text(terminal_stop_sequence_count > 0),
        "terminal_stop_sequence_count": terminal_stop_sequence_count,
        "has_internal_stop_codon": _bool_text(has_internal_stop_codon),
        "has_invalid_amino_acid": _bool_text(bool(invalid_chars)),
        "invalid_amino_acids": "".join(sorted(invalid_chars)),
    }


def _notes_for_cell(
    status: str,
    has_internal_stop_codon: bool,
    has_invalid_amino_acid: bool,
    sequence_count_mismatch: bool,
) -> str:
    notes: list[str] = []
    if status == "Duplicated":
        notes.append("duplicated_status")
    elif status == "Fragmented":
        notes.append("fragmented_status")
    elif status == "Missing":
        notes.append("missing_status")
    elif status == "NotReported":
        notes.append("not_reported")
    elif status == "Conflict":
        notes.append("conflicting_status")

    if has_internal_stop_codon:
        notes.append("internal_stop_codon")
    if has_invalid_amino_acid:
        notes.append("invalid_amino_acid")
    if sequence_count_mismatch:
        notes.append("sequence_count_mismatch")
    return ",".join(notes)


def summarize_locus_sample_cell(
    locus_id: str,
    sample_row: dict[str, str],
    record_rows: list[dict[str, str]],
    repo_root: Path,
) -> dict[str, object]:
    if not record_rows:
        return {
            "locus_id": locus_id,
            **sample_row,
            "status": "NotReported",
            "record_count": 0,
            "sequence_count": 0,
            "sequence_count_mismatch": "false",
            "full_table_sequence_ids": "",
            "gene_starts": "",
            "gene_ends": "",
            "strands": "",
            "scores": "",
            "match_lengths_aa": "",
            "faa_path": "",
            "gff_path": "",
            "protein_lengths_aa": "",
            "representative_protein_length_aa": "",
            "has_terminal_stop_codon": "false",
            "terminal_stop_sequence_count": 0,
            "has_internal_stop_codon": "false",
            "has_invalid_amino_acid": "false",
            "invalid_amino_acids": "",
            "is_complete_single_copy": "false",
            "include_in_occupancy": "false",
            "qc_notes": "not_reported",
        }

    statuses = sorted({row["status"] for row in record_rows if row["status"]})
    status = statuses[0] if len(statuses) == 1 else "Conflict"

    faa_paths = sorted({row["faa_path"] for row in record_rows if row["faa_path"]})
    gff_paths = sorted({row["gff_path"] for row in record_rows if row["gff_path"]})
    if len(faa_paths) > 1 or len(gff_paths) > 1:
        raise ValueError(f"Sample {sample_row['sample_id']} locus {locus_id} has inconsistent BUSCO paths.")

    faa_path = faa_paths[0] if faa_paths else ""
    gff_path = gff_paths[0] if gff_paths else ""
    protein_stats = analyze_protein_sequences(faa_path, repo_root) if faa_path else analyze_protein_sequences("", repo_root)

    record_count = len(record_rows)
    sequence_count = int(protein_stats["sequence_count"])
    sequence_count_mismatch = status not in {"Missing", "NotReported"} and record_count != sequence_count
    is_complete_single_copy = status == "Complete" and record_count == 1 and sequence_count == 1

    return {
        "locus_id": locus_id,
        **sample_row,
        "status": status,
        "record_count": record_count,
        "sequence_count": sequence_count,
        "sequence_count_mismatch": _bool_text(sequence_count_mismatch),
        "full_table_sequence_ids": _join(row["sequence_id"] for row in record_rows),
        "gene_starts": _join(row["gene_start"] for row in record_rows),
        "gene_ends": _join(row["gene_end"] for row in record_rows),
        "strands": _join(row["strand"] for row in record_rows),
        "scores": _join(row["score"] for row in record_rows),
        "match_lengths_aa": _join(row["length"] for row in record_rows),
        "faa_path": faa_path,
        "gff_path": gff_path,
        **protein_stats,
        "is_complete_single_copy": _bool_text(is_complete_single_copy),
        "include_in_occupancy": _bool_text(is_complete_single_copy),
        "qc_notes": _notes_for_cell(
            status,
            _as_bool(str(protein_stats["has_internal_stop_codon"])),
            _as_bool(str(protein_stats["has_invalid_amino_acid"])),
            sequence_count_mismatch,
        ),
    }


def _locus_metrics(rows: list[dict[str, object]]) -> dict[str, object]:
    status_counter = Counter(str(row["status"]) for row in rows)
    total_taxa = len(rows)
    complete_single_copy_taxa = sum(_as_bool(str(row["include_in_occupancy"])) for row in rows)
    stop_codon_taxa = sum(_as_bool(str(row["has_internal_stop_codon"])) for row in rows)
    invalid_amino_acid_taxa = sum(_as_bool(str(row["has_invalid_amino_acid"])) for row in rows)
    occupancy = complete_single_copy_taxa / total_taxa if total_taxa else 0.0

    return {
        "locus_total_taxa": total_taxa,
        "locus_complete_single_copy_taxa": complete_single_copy_taxa,
        "locus_duplicated_taxa": status_counter["Duplicated"],
        "locus_fragmented_taxa": status_counter["Fragmented"],
        "locus_missing_taxa": status_counter["Missing"] + status_counter["NotReported"],
        "locus_stop_codon_taxa": stop_codon_taxa,
        "locus_invalid_amino_acid_taxa": invalid_amino_acid_taxa,
        "locus_occupancy": _format_occupancy(occupancy),
    }


def build_locus_taxon_matrix_rows(
    summary_rows: list[dict[str, str]],
    record_rows: list[dict[str, str]],
    repo_root: Path,
) -> list[dict[str, object]]:
    sample_rows = [
        {
            "sample_id": row["sample_id"],
            "taxon_id": row["taxon_id"],
            "sanitized_taxon_id": row["sanitized_taxon_id"],
            "assembly_fasta": row["assembly_fasta"],
        }
        for row in summary_rows
    ]
    sample_order = [row["sample_id"] for row in sample_rows]
    sample_map = {row["sample_id"]: row for row in sample_rows}

    grouped_records: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in record_rows:
        grouped_records[(row["busco_id"], row["sample_id"])].append(row)

    locus_ids = sorted({row["busco_id"] for row in record_rows}, key=_busco_sort_key)
    matrix_rows: list[dict[str, object]] = []

    for locus_id in locus_ids:
        locus_rows = [
            summarize_locus_sample_cell(
                locus_id,
                sample_map[sample_id],
                grouped_records.get((locus_id, sample_id), []),
                repo_root,
            )
            for sample_id in sample_order
        ]
        metrics = _locus_metrics(locus_rows)
        for row in locus_rows:
            row.update(metrics)
            matrix_rows.append(row)

    return matrix_rows


def build_retained_loci_rows(
    matrix_rows: list[dict[str, str | object]],
    occupancy_threshold: float,
) -> list[dict[str, object]]:
    grouped_rows: dict[str, list[dict[str, str | object]]] = defaultdict(list)
    for row in matrix_rows:
        grouped_rows[str(row["locus_id"])].append(row)

    retained_rows: list[dict[str, object]] = []
    for locus_id in sorted(grouped_rows, key=_busco_sort_key):
        rows = grouped_rows[locus_id]
        first_row = rows[0]
        total_taxa = int(first_row["locus_total_taxa"])
        complete_single_copy_taxa = int(first_row["locus_complete_single_copy_taxa"])
        duplicated_taxa = int(first_row["locus_duplicated_taxa"])
        fragmented_taxa = int(first_row["locus_fragmented_taxa"])
        missing_taxa = int(first_row["locus_missing_taxa"])
        stop_codon_taxa = int(first_row["locus_stop_codon_taxa"])
        invalid_amino_acid_taxa = int(first_row["locus_invalid_amino_acid_taxa"])
        occupancy = complete_single_copy_taxa / total_taxa if total_taxa else 0.0
        passes_occupancy = occupancy >= occupancy_threshold

        retained_rows_for_locus = sorted(
            [row for row in rows if _as_bool(str(row["include_in_occupancy"]))],
            key=lambda row: str(row["sanitized_taxon_id"]),
        )
        retained_lengths = [
            int(str(row["representative_protein_length_aa"]))
            for row in retained_rows_for_locus
            if str(row["representative_protein_length_aa"])
        ]

        min_length = min(retained_lengths) if retained_lengths else ""
        max_length = max(retained_lengths) if retained_lengths else ""
        median_length = median(retained_lengths) if retained_lengths else ""
        length_span = max_length - min_length if retained_lengths else ""
        length_ratio = (max_length / min_length) if retained_lengths and min_length else ""
        length_dispersion_observed = len(set(retained_lengths)) > 1

        failure_reasons: list[str] = []
        if complete_single_copy_taxa == 0:
            failure_reasons.append("no_complete_single_copy_hits")
        if not passes_occupancy:
            failure_reasons.append("occupancy_below_threshold")

        qc_warnings: list[str] = []
        if duplicated_taxa > 0:
            qc_warnings.append("duplicated_taxa_present")
        if fragmented_taxa > 0:
            qc_warnings.append("fragmented_taxa_present")
        if missing_taxa > 0:
            qc_warnings.append("missing_taxa_present")
        if stop_codon_taxa > 0:
            qc_warnings.append("internal_stop_codon_present")
        if invalid_amino_acid_taxa > 0:
            qc_warnings.append("invalid_amino_acid_present")
        if length_dispersion_observed:
            qc_warnings.append("length_dispersion_observed")

        retained_rows.append(
            {
                "locus_id": locus_id,
                "total_taxa": total_taxa,
                "complete_single_copy_taxa": complete_single_copy_taxa,
                "duplicated_taxa": duplicated_taxa,
                "fragmented_taxa": fragmented_taxa,
                "missing_taxa": missing_taxa,
                "stop_codon_taxa": stop_codon_taxa,
                "invalid_amino_acid_taxa": invalid_amino_acid_taxa,
                "occupancy": _format_occupancy(occupancy),
                "occupancy_threshold": _format_occupancy(occupancy_threshold),
                "passes_occupancy": _bool_text(passes_occupancy),
                "decision": "retain" if passes_occupancy else "exclude",
                "failure_reasons": ",".join(failure_reasons),
                "qc_warnings": ",".join(qc_warnings),
                "retained_sample_ids": _join(str(row["sample_id"]) for row in retained_rows_for_locus),
                "retained_sanitized_taxon_ids": _join(
                    str(row["sanitized_taxon_id"]) for row in retained_rows_for_locus
                ),
                "retained_protein_lengths_aa": _join(str(length) for length in retained_lengths),
                "min_protein_length_aa": _format_stat(min_length),
                "median_protein_length_aa": _format_stat(median_length),
                "max_protein_length_aa": _format_stat(max_length),
                "length_span_aa": _format_stat(length_span),
                "length_ratio": _format_stat(length_ratio),
                "length_dispersion_observed": _bool_text(length_dispersion_observed),
            }
        )

    return retained_rows


def write_locus_taxon_matrix(path: Path, rows: list[dict[str, object]]) -> None:
    write_tsv(path, rows, MATRIX_FIELDNAMES)


def write_retained_loci(path: Path, rows: list[dict[str, object]]) -> None:
    write_tsv(path, rows, RETAINED_FIELDNAMES)
