"""Helpers for BUSCO execution, parsing, and QC summarization."""

from __future__ import annotations

import csv
import json
import re
import shutil
from collections import Counter
from pathlib import Path

from .busco_sequences import (
    CATEGORY_BY_STATUS,
    CONCATENATED_LAYOUT,
    FAA_ARTIFACT_BY_CATEGORY,
    GFF_ARTIFACT_BY_CATEGORY,
    LEGACY_DIR_BY_CATEGORY,
    artifact_layout,
    concatenate_busco_category,
    concatenated_category_paths,
    has_fasta_records_for_id,
    has_gff_block,
    is_concatenated_artifact,
)
from .manifest import write_tsv

VALIDATED_MANIFEST_COLUMNS = ("sample_id", "taxon_id", "sanitized_taxon_id", "assembly_fasta")
BUSCO_SEQUENCE_ARTIFACTS = {
    status: FAA_ARTIFACT_BY_CATEGORY[category] for status, category in CATEGORY_BY_STATUS.items()
}
BUSCO_SUMMARY_FIELDNAMES = [
    "sample_id",
    "taxon_id",
    "sanitized_taxon_id",
    "assembly_fasta",
    "busco_version",
    "lineage_name",
    "lineage_creation_date",
    "lineage_busco_count",
    "lineage_species_count",
    "input_fasta",
    "busco_mode",
    "gene_predictor",
    "one_line_summary",
    "n_markers",
    "complete_buscos",
    "single_copy_buscos",
    "multi_copy_buscos",
    "fragmented_buscos",
    "missing_buscos",
    "complete_percent",
    "single_copy_percent",
    "multi_copy_percent",
    "fragmented_percent",
    "missing_percent",
    "avg_identity",
    "internal_stop_codon_count",
    "internal_stop_codon_percent",
    "full_table_row_count",
    "full_table_single_copy_buscos",
    "full_table_multi_copy_buscos",
    "full_table_fragmented_buscos",
    "full_table_missing_buscos",
    "full_table_duplicated_rows",
    "full_table_consistent",
    "raw_root",
    "paths_tsv",
    "short_summary_path",
    "full_table_path",
    "run_complete_path",
    "sequence_layout",
    "single_copy_sequence_path",
    "multi_copy_sequence_path",
    "fragmented_sequence_path",
]
BUSCO_RECORD_FIELDNAMES = [
    "sample_id",
    "taxon_id",
    "sanitized_taxon_id",
    "assembly_fasta",
    "busco_id",
    "status",
    "sequence_id",
    "gene_start",
    "gene_end",
    "strand",
    "score",
    "length",
    "orthodb_url",
    "description",
    "sequence_category",
    "sequence_dir",
    "faa_path",
    "gff_path",
]
BUSCO_FULL_TABLE_COLUMNS = [
    "busco_id",
    "status",
    "sequence_id",
    "gene_start",
    "gene_end",
    "strand",
    "score",
    "length",
    "orthodb_url",
    "description",
]


class BuscoValidationError(ValueError):
    """Raised when BUSCO metadata is missing or invalid."""


def parse_dataset_names(text: str) -> list[str]:
    dataset_names: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        match = re.search(r"([A-Za-z0-9._-]+_odb\d+)", line)
        if match:
            dataset_names.append(match.group(1))
    return dataset_names


def load_dataset_names(path: Path) -> list[str]:
    return parse_dataset_names(path.read_text(encoding="utf-8"))


def verify_lineage_name(lineage: str, available_datasets: list[str]) -> str:
    normalized = (lineage or "").strip()
    if not normalized:
        raise BuscoValidationError("BUSCO lineage is not configured.")
    if normalized not in available_datasets:
        preview = ", ".join(sorted(available_datasets)[:10])
        raise BuscoValidationError(
            f"Configured BUSCO lineage {normalized!r} was not found in the available dataset "
            f"list. First available entries: {preview}"
        )
    return normalized


def _resolve_repo_path(repo_root: Path, path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else repo_root / path


def _load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def load_validated_manifest_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = tuple(reader.fieldnames or ())
        missing = [name for name in VALIDATED_MANIFEST_COLUMNS if name not in fieldnames]
        if missing:
            raise BuscoValidationError(
                "Validated manifest is missing required columns: " + ", ".join(missing)
            )
        return [dict(row) for row in reader]


def busco_output_paths(sample_id: str) -> dict[str, str]:
    sample_root = Path("results") / "busco" / sample_id
    sequence_root = sample_root / "busco_sequences"
    paths = {
        "sample_root": sample_root.as_posix(),
        "raw_root": (Path("work") / "busco" / sample_id / "raw").as_posix(),
        "sequence_root": sequence_root.as_posix(),
        "command": (sample_root / "command.sh").as_posix(),
        "paths": (sample_root / "paths.tsv").as_posix(),
        "short_summary": (sample_root / "short_summary.txt").as_posix(),
        "full_table": (sample_root / "full_table.tsv").as_posix(),
        "completion": (sample_root / "run.complete").as_posix(),
    }

    for category, legacy_dir in LEGACY_DIR_BY_CATEGORY.items():
        faa_path, gff_path = concatenated_category_paths(sequence_root, category)
        paths[f"{category}_faa"] = faa_path.as_posix()
        paths[f"{category}_gff"] = gff_path.as_posix()
        # Where runs predating the concatenated layout keep the same sequences.
        paths[f"{category}_per_locus_dir"] = (sequence_root / legacy_dir).as_posix()

    return paths


def load_busco_artifact_paths(path: Path) -> dict[str, str]:
    rows = _load_tsv_rows(path)
    if not rows:
        raise BuscoValidationError(f"BUSCO artifact index is empty: {path}")
    artifact_paths: dict[str, str] = {}
    for row in rows:
        artifact = (row.get("artifact") or "").strip()
        artifact_path = (row.get("path") or "").strip()
        if not artifact:
            raise BuscoValidationError(f"BUSCO artifact index contains an empty artifact name: {path}")
        artifact_paths[artifact] = artifact_path
    return artifact_paths


def _require_mapping(payload: dict[str, object], key: str, context: str) -> dict[str, object]:
    value = payload.get(key)
    if not isinstance(value, dict):
        raise BuscoValidationError(f"BUSCO {context} is missing mapping {key!r}.")
    return value


def _require_value(payload: dict[str, object], key: str, context: str) -> object:
    if key not in payload:
        raise BuscoValidationError(f"BUSCO {context} is missing required field {key!r}.")
    return payload[key]


def parse_busco_short_summary(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise BuscoValidationError(f"BUSCO short summary is not valid JSON: {path}") from exc
    if not isinstance(payload, dict):
        raise BuscoValidationError(f"BUSCO short summary must be a JSON object: {path}")

    parameters = _require_mapping(payload, "parameters", "short summary")
    lineage = _require_mapping(payload, "lineage_dataset", "short summary")
    versions = _require_mapping(payload, "versions", "short summary")
    results = _require_mapping(payload, "results", "short summary")

    return {
        "busco_version": str(_require_value(versions, "busco", "short summary.versions")),
        "lineage_name": str(_require_value(lineage, "name", "short summary.lineage_dataset")),
        "lineage_creation_date": str(
            _require_value(lineage, "creation_date", "short summary.lineage_dataset")
        ),
        "lineage_busco_count": int(
            _require_value(lineage, "number_of_buscos", "short summary.lineage_dataset")
        ),
        "lineage_species_count": int(
            _require_value(lineage, "number_of_species", "short summary.lineage_dataset")
        ),
        "input_fasta": str(_require_value(parameters, "in", "short summary.parameters")),
        "busco_mode": str(_require_value(parameters, "mode", "short summary.parameters")),
        "gene_predictor": str(parameters.get("gene_predictor", "")),
        "one_line_summary": str(_require_value(results, "one_line_summary", "short summary.results")),
        "n_markers": int(_require_value(results, "n_markers", "short summary.results")),
        "complete_buscos": int(_require_value(results, "Complete BUSCOs", "short summary.results")),
        "single_copy_buscos": int(
            _require_value(results, "Single copy BUSCOs", "short summary.results")
        ),
        "multi_copy_buscos": int(
            _require_value(results, "Multi copy BUSCOs", "short summary.results")
        ),
        "fragmented_buscos": int(
            _require_value(results, "Fragmented BUSCOs", "short summary.results")
        ),
        "missing_buscos": int(_require_value(results, "Missing BUSCOs", "short summary.results")),
        "complete_percent": float(
            _require_value(results, "Complete percentage", "short summary.results")
        ),
        "single_copy_percent": float(
            _require_value(results, "Single copy percentage", "short summary.results")
        ),
        "multi_copy_percent": float(
            _require_value(results, "Multi copy percentage", "short summary.results")
        ),
        "fragmented_percent": float(
            _require_value(results, "Fragmented percentage", "short summary.results")
        ),
        "missing_percent": float(
            _require_value(results, "Missing percentage", "short summary.results")
        ),
        "avg_identity": float(results.get("avg_identity", 0.0)),
        "internal_stop_codon_count": int(results.get("internal_stop_codon_count", 0)),
        "internal_stop_codon_percent": float(results.get("internal_stop_codon_percent", 0.0)),
    }


def parse_busco_full_table(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    header_seen = False
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.rstrip()
        if not line:
            continue
        if line.startswith("# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength"):
            header_seen = True
            continue
        if line.startswith("#"):
            continue
        parts = line.split("\t", len(BUSCO_FULL_TABLE_COLUMNS) - 1)
        if len(parts) < 2:
            raise BuscoValidationError(f"Malformed BUSCO full_table row in {path}: {line!r}")
        if len(parts) < len(BUSCO_FULL_TABLE_COLUMNS):
            parts.extend([""] * (len(BUSCO_FULL_TABLE_COLUMNS) - len(parts)))
        rows.append(dict(zip(BUSCO_FULL_TABLE_COLUMNS, parts, strict=True)))
    if not header_seen:
        raise BuscoValidationError(f"BUSCO full_table header was not found in {path}")
    return rows


def summarize_busco_full_table(rows: list[dict[str, str]]) -> dict[str, int]:
    status_counts: Counter[str] = Counter()
    status_buscos = {status: set() for status in ("Complete", "Duplicated", "Fragmented", "Missing")}

    for row in rows:
        status = row["status"]
        busco_id = row["busco_id"]
        if status not in status_buscos:
            raise BuscoValidationError(f"Unsupported BUSCO status in full_table: {status!r}")
        status_counts[status] += 1
        status_buscos[status].add(busco_id)

    return {
        "full_table_row_count": len(rows),
        "full_table_single_copy_buscos": len(status_buscos["Complete"]),
        "full_table_multi_copy_buscos": len(status_buscos["Duplicated"]),
        "full_table_fragmented_buscos": len(status_buscos["Fragmented"]),
        "full_table_missing_buscos": len(status_buscos["Missing"]),
        "full_table_duplicated_rows": status_counts["Duplicated"],
    }


def validate_busco_summary_consistency(
    summary: dict[str, object],
    table_counts: dict[str, int],
    sample_id: str,
) -> None:
    checks = {
        "single_copy_buscos": "full_table_single_copy_buscos",
        "multi_copy_buscos": "full_table_multi_copy_buscos",
        "fragmented_buscos": "full_table_fragmented_buscos",
        "missing_buscos": "full_table_missing_buscos",
    }
    for summary_key, table_key in checks.items():
        summary_value = int(summary[summary_key])
        table_value = int(table_counts[table_key])
        if summary_value != table_value:
            raise BuscoValidationError(
                "BUSCO summary/full_table mismatch for sample "
                f"{sample_id!r}: {summary_key}={summary_value} but {table_key}={table_value}"
            )
    complete_buscos = int(summary["complete_buscos"])
    if complete_buscos != int(summary["single_copy_buscos"]) + int(summary["multi_copy_buscos"]):
        raise BuscoValidationError(
            f"BUSCO complete counts are inconsistent for sample {sample_id!r}."
        )


def resolve_busco_sequence_paths(
    status: str,
    busco_id: str,
    artifact_paths: dict[str, str],
    repo_root: Path,
) -> dict[str, str]:
    """Address one locus's sequences, in either the concatenated or per-locus layout."""
    category = CATEGORY_BY_STATUS.get(status, "")
    if not category:
        return {"sequence_category": "", "sequence_dir": "", "faa_path": "", "gff_path": ""}

    artifact_name = FAA_ARTIFACT_BY_CATEGORY[category]
    recorded_path = artifact_paths.get(artifact_name, "")
    if not recorded_path:
        raise BuscoValidationError(
            f"BUSCO artifact index does not provide {artifact_name!r} for status {status!r}."
        )

    if is_concatenated_artifact(recorded_path):
        faa_path = recorded_path
        gff_path = artifact_paths.get(GFF_ARTIFACT_BY_CATEGORY[category], "")
        if not gff_path:
            raise BuscoValidationError(
                f"BUSCO artifact index does not provide {GFF_ARTIFACT_BY_CATEGORY[category]!r} "
                f"for status {status!r}."
            )
        sequence_dir = Path(recorded_path).parent.as_posix()
    else:
        sequence_dir = recorded_path
        faa_path = (Path(recorded_path) / f"{busco_id}.faa").as_posix()
        gff_path = (Path(recorded_path) / f"{busco_id}.gff").as_posix()

    if not has_fasta_records_for_id(faa_path, busco_id, repo_root) or not has_gff_block(
        gff_path, busco_id, repo_root
    ):
        raise BuscoValidationError(
            f"Expected BUSCO sequence records for {busco_id!r} in {faa_path!r} and {gff_path!r}."
        )

    return {
        "sequence_category": artifact_name,
        "sequence_dir": sequence_dir,
        "faa_path": faa_path,
        "gff_path": gff_path,
    }


def hydrate_stable_busco_sequence_paths(
    artifact_paths: dict[str, str],
    stable_paths: dict[str, str],
    repo_root: Path,
) -> dict[str, str]:
    """Prefer the copies under `results/busco/` over the scratch paths BUSCO reported."""
    hydrated = dict(artifact_paths)
    for category in LEGACY_DIR_BY_CATEGORY:
        faa_artifact = FAA_ARTIFACT_BY_CATEGORY[category]
        gff_artifact = GFF_ARTIFACT_BY_CATEGORY[category]
        stable_faa = stable_paths[f"{category}_faa"]
        stable_per_locus_dir = stable_paths[f"{category}_per_locus_dir"]

        if _resolve_repo_path(repo_root, stable_faa).is_file():
            hydrated[faa_artifact] = stable_faa
            hydrated[gff_artifact] = stable_paths[f"{category}_gff"]
        elif _resolve_repo_path(repo_root, stable_per_locus_dir).is_dir():
            hydrated[faa_artifact] = stable_per_locus_dir
            hydrated.pop(gff_artifact, None)
    return hydrated


def _sequence_path_columns(artifact_paths: dict[str, str]) -> dict[str, str]:
    columns: dict[str, str] = {}
    layouts = set()
    for category in LEGACY_DIR_BY_CATEGORY:
        recorded_path = artifact_paths.get(FAA_ARTIFACT_BY_CATEGORY[category], "")
        columns[f"{category}_sequence_path"] = recorded_path
        if recorded_path:
            layouts.add(artifact_layout(recorded_path))
    columns["sequence_layout"] = "/".join(sorted(layouts))
    return columns


def build_busco_tables(validated_manifest_path: Path, repo_root: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    summary_rows: list[dict[str, object]] = []
    record_rows: list[dict[str, object]] = []

    for sample_row in load_validated_manifest_rows(validated_manifest_path):
        sample_id = sample_row["sample_id"]
        stable_paths = busco_output_paths(sample_id)

        required_paths = {
            name: _resolve_repo_path(repo_root, stable_paths[name])
            for name in ("paths", "short_summary", "full_table", "completion")
        }
        for name, path in required_paths.items():
            if not path.is_file():
                raise BuscoValidationError(
                    f"Expected BUSCO {name} file for sample {sample_id!r}: {path}"
                )

        artifact_paths = hydrate_stable_busco_sequence_paths(
            load_busco_artifact_paths(required_paths["paths"]),
            stable_paths,
            repo_root,
        )
        summary = parse_busco_short_summary(required_paths["short_summary"])
        full_table_rows = parse_busco_full_table(required_paths["full_table"])
        table_counts = summarize_busco_full_table(full_table_rows)
        validate_busco_summary_consistency(summary, table_counts, sample_id)

        summary_rows.append(
            {
                **sample_row,
                **summary,
                **table_counts,
                "full_table_consistent": "true",
                "raw_root": artifact_paths.get("raw_root", ""),
                "paths_tsv": stable_paths["paths"],
                "short_summary_path": stable_paths["short_summary"],
                "full_table_path": stable_paths["full_table"],
                "run_complete_path": stable_paths["completion"],
                **_sequence_path_columns(artifact_paths),
            }
        )

        for row in full_table_rows:
            sequence_paths = resolve_busco_sequence_paths(
                row["status"],
                row["busco_id"],
                artifact_paths,
                repo_root,
            )
            record_rows.append({**sample_row, **row, **sequence_paths})

    return summary_rows, record_rows


def _find_single_path(base_dir: Path, pattern: str) -> Path:
    matches = sorted(base_dir.rglob(pattern))
    if not matches:
        raise BuscoValidationError(f"Expected BUSCO artifact {pattern!r} under {base_dir}")
    return matches[0]


def _find_optional_dir(base_dir: Path, name: str) -> str:
    matches = sorted(path for path in base_dir.rglob(name) if path.is_dir())
    return matches[0].as_posix() if matches else ""


def _concatenate_optional_category(
    raw_root: Path,
    sequence_root: Path,
    category: str,
) -> dict[str, str]:
    """Fold one BUSCO output directory into `<category>.faa` / `<category>.gff`.

    BUSCO emits two files per locus, which is two inodes per locus per sample;
    the workflow keeps the same records addressed by BUSCO id instead.
    """
    source_dir_text = _find_optional_dir(raw_root, LEGACY_DIR_BY_CATEGORY[category])
    if not source_dir_text:
        return {"faa_path": "", "gff_path": ""}
    return concatenate_busco_category(Path(source_dir_text), sequence_root, category)


def standardize_busco_run(sample_id: str, sample_dir: Path, raw_root: Path) -> None:
    short_summary_source = _find_single_path(raw_root, "short_summary*")
    full_table_source = _find_single_path(raw_root, "full_table.tsv")

    stable_paths = busco_output_paths(sample_id)
    short_summary_target = Path(stable_paths["short_summary"])
    full_table_target = Path(stable_paths["full_table"])
    paths_target = Path(stable_paths["paths"])
    completion_target = Path(stable_paths["completion"])
    sequence_root = Path(stable_paths["sequence_root"])

    short_summary_target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(short_summary_source, short_summary_target)
    shutil.copy2(full_table_source, full_table_target)
    if sequence_root.exists():
        shutil.rmtree(sequence_root)
    sequence_root.mkdir(parents=True, exist_ok=True)

    path_rows = [
        {"artifact": "raw_root", "path": raw_root.as_posix()},
        {"artifact": "short_summary_source", "path": short_summary_source.as_posix()},
        {"artifact": "full_table_source", "path": full_table_source.as_posix()},
        {"artifact": "sequence_layout", "path": CONCATENATED_LAYOUT},
    ]
    for category in LEGACY_DIR_BY_CATEGORY:
        written = _concatenate_optional_category(raw_root, sequence_root, category)
        path_rows.append(
            {"artifact": FAA_ARTIFACT_BY_CATEGORY[category], "path": written["faa_path"]}
        )
        path_rows.append(
            {"artifact": GFF_ARTIFACT_BY_CATEGORY[category], "path": written["gff_path"]}
        )

    write_tsv(paths_target, path_rows, ["artifact", "path"])
    completion_target.write_text(
        f"sample_id\t{sample_id}\nraw_root\t{raw_root.as_posix()}\n",
        encoding="utf-8",
    )
