"""Helpers for extracting retained-locus BUSCO CDS sequences with gffread."""

from __future__ import annotations

import csv
import shlex
import subprocess
from pathlib import Path

from .locus_matrix import parse_fasta_records
from .manifest import write_tsv


DNA_RECORD_FIELDNAMES = [
    "sample_id",
    "locus_id",
    "sanitized_taxon_id",
    "transcript_id",
    "source_gff_path",
    "per_locus_fasta",
    "sequence_length_nt",
]


def busco_sort_key(locus_id: str) -> tuple[int, str]:
    prefix, _, suffix = locus_id.partition("at")
    if prefix.isdigit():
        return int(prefix), suffix
    return 0, locus_id


def sample_dna_output_paths(sample_id: str) -> dict[str, str]:
    root = Path("results") / "dna_sequences" / sample_id
    return {
        "root": root.as_posix(),
        "aggregated_gff": (root / "retained_loci.gff3").as_posix(),
        "extracted_fasta": (root / "retained_loci.fna").as_posix(),
        "records": (root / "retained_loci.records.tsv").as_posix(),
        "command": (root / "gffread.command.sh").as_posix(),
        "log": (root / "gffread.log").as_posix(),
        "per_locus_dir": (root / "per_locus").as_posix(),
        "completion": (root / "run.complete").as_posix(),
    }


def per_locus_dna_sequence_path(sample_id: str, locus_id: str) -> str:
    root = Path(sample_dna_output_paths(sample_id)["per_locus_dir"])
    return (root / f"{locus_id}.fna").as_posix()


def _resolve_repo_path(repo_root: Path, path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else repo_root / path


def _load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def _parse_gff_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for part in text.strip().strip(";").split(";"):
        item = part.strip()
        if not item:
            continue
        if "=" not in item:
            attributes[item] = ""
            continue
        key, value = item.split("=", 1)
        attributes[key] = value
    return attributes


def _format_gff_attributes(attributes: dict[str, str]) -> str:
    parts = []
    for key, value in attributes.items():
        if value == "":
            parts.append(key)
        else:
            parts.append(f"{key}={value}")
    return ";".join(parts)


def _read_retained_locus_ids(retained_path: Path) -> set[str]:
    rows = _load_tsv_rows(retained_path)
    return {row["locus_id"] for row in rows if row.get("decision") == "retain"}


def _select_sample_rows(
    *,
    matrix_path: Path,
    retained_path: Path,
    sample_id: str,
) -> list[dict[str, str]]:
    retained_locus_ids = _read_retained_locus_ids(retained_path)
    rows = [
        row
        for row in _load_tsv_rows(matrix_path)
        if row.get("sample_id") == sample_id
        and row.get("locus_id") in retained_locus_ids
        and row.get("include_in_occupancy", "").strip().lower() == "true"
    ]
    return sorted(rows, key=lambda row: busco_sort_key(row["locus_id"]))


def _rewrite_locus_gff(
    *,
    locus_id: str,
    source_path: Path,
    destination_handle,
    transcript_id: str,
) -> None:
    transcript_seen = False
    for raw_line in source_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) != 9:
            raise ValueError(f"GFF record in {source_path} does not have 9 columns: {raw_line!r}")
        feature_type = columns[2]
        attributes = _parse_gff_attributes(columns[8])

        if feature_type in {"mRNA", "transcript"}:
            attributes["ID"] = transcript_id
            attributes.pop("Parent", None)
            transcript_seen = True
        elif feature_type in {"CDS", "exon", "stop_codon", "start_codon"}:
            attributes["Parent"] = transcript_id
        else:
            continue

        columns[8] = _format_gff_attributes(attributes)
        destination_handle.write("\t".join(columns) + "\n")

    if not transcript_seen:
        raise ValueError(f"BUSCO GFF for locus {locus_id!r} does not contain an mRNA/transcript feature: {source_path}")


def _write_split_fasta(output_path: Path, header: str, sequence: str) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(f">{header}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start : start + 80] + "\n")


def build_gffread_command(
    *,
    gffread_executable: str,
    annotation_path: Path,
    genome_fasta_path: Path,
    output_fasta_path: Path,
) -> list[str]:
    return [
        gffread_executable,
        annotation_path.as_posix(),
        "-g",
        genome_fasta_path.as_posix(),
        "-x",
        output_fasta_path.as_posix(),
    ]


def write_command_script(path: Path, command: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rendered = " ".join(shlex.quote(token) for token in command)
    path.write_text(f"#!/usr/bin/env bash\nset -euo pipefail\n{rendered}\n", encoding="utf-8")
    path.chmod(0o755)


def extract_retained_sample_dna(
    *,
    matrix_path: Path,
    retained_path: Path,
    sample_id: str,
    genome_fasta_path: Path,
    repo_root: Path,
    gffread_executable: str = "gffread",
) -> list[dict[str, str]]:
    output_paths = sample_dna_output_paths(sample_id)
    aggregated_gff_path = Path(output_paths["aggregated_gff"])
    extracted_fasta_path = Path(output_paths["extracted_fasta"])
    records_path = Path(output_paths["records"])
    per_locus_dir = Path(output_paths["per_locus_dir"])
    command_path = Path(output_paths["command"])
    log_path = Path(output_paths["log"])
    completion_path = Path(output_paths["completion"])

    selected_rows = _select_sample_rows(
        matrix_path=matrix_path,
        retained_path=retained_path,
        sample_id=sample_id,
    )
    if not selected_rows:
        raise ValueError(f"No retained complete single-copy loci were found for sample {sample_id!r}.")

    aggregated_gff_path.parent.mkdir(parents=True, exist_ok=True)
    per_locus_dir.mkdir(parents=True, exist_ok=True)
    with aggregated_gff_path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for row in selected_rows:
            locus_id = row["locus_id"]
            gff_path_text = row.get("gff_path", "").strip()
            if not gff_path_text:
                raise ValueError(f"Matrix row for sample {sample_id!r} locus {locus_id!r} is missing gff_path.")
            source_gff_path = _resolve_repo_path(repo_root, gff_path_text)
            if not source_gff_path.is_file():
                raise ValueError(f"BUSCO GFF does not exist for sample {sample_id!r} locus {locus_id!r}: {source_gff_path}")
            _rewrite_locus_gff(
                locus_id=locus_id,
                source_path=source_gff_path,
                destination_handle=handle,
                transcript_id=locus_id,
            )

    command = build_gffread_command(
        gffread_executable=gffread_executable,
        annotation_path=aggregated_gff_path,
        genome_fasta_path=genome_fasta_path,
        output_fasta_path=extracted_fasta_path,
    )
    write_command_script(command_path, command)
    with log_path.open("w", encoding="utf-8") as log_handle:
        subprocess.run(command, check=True, stdout=log_handle, stderr=log_handle)

    sequence_by_transcript = {
        header.split()[0]: sequence.upper()
        for header, sequence in parse_fasta_records(extracted_fasta_path)
    }

    expected_files = set()
    record_rows: list[dict[str, str]] = []
    for row in selected_rows:
        locus_id = row["locus_id"]
        sequence = sequence_by_transcript.get(locus_id)
        if sequence is None:
            raise ValueError(
                f"gffread output for sample {sample_id!r} does not contain expected transcript {locus_id!r}."
            )
        output_fasta_path = Path(per_locus_dna_sequence_path(sample_id, locus_id))
        _write_split_fasta(output_fasta_path, row["sanitized_taxon_id"], sequence)
        expected_files.add(output_fasta_path.name)
        record_rows.append(
            {
                "sample_id": sample_id,
                "locus_id": locus_id,
                "sanitized_taxon_id": row["sanitized_taxon_id"],
                "transcript_id": locus_id,
                "source_gff_path": row["gff_path"],
                "per_locus_fasta": output_fasta_path.as_posix(),
                "sequence_length_nt": str(len(sequence)),
            }
        )

    for path in per_locus_dir.glob("*.fna"):
        if path.name not in expected_files:
            path.unlink()

    write_tsv(records_path, record_rows, DNA_RECORD_FIELDNAMES)
    completion_path.write_text(f"sample_id\t{sample_id}\nrecords\t{len(record_rows)}\n", encoding="utf-8")
    return record_rows
