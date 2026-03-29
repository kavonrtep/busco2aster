"""Helpers for final workflow QC reporting."""

from __future__ import annotations

import csv
from collections import Counter
from datetime import datetime
from pathlib import Path

from .concordance import summarize_gcf_stat
from .gene_trees import read_single_line_tree


BUSCO_REQUIRED_COLUMNS = {
    "sample_id",
    "busco_version",
    "lineage_name",
    "complete_percent",
    "single_copy_percent",
    "multi_copy_percent",
    "fragmented_percent",
    "missing_percent",
    "internal_stop_codon_count",
}

RETAINED_LOCI_REQUIRED_COLUMNS = {
    "locus_id",
    "decision",
    "occupancy",
    "occupancy_threshold",
    "duplicated_taxa",
    "stop_codon_taxa",
    "failure_reasons",
    "qc_warnings",
}

GENE_TREE_REQUIRED_COLUMNS = {
    "locus_id",
    "selected_model",
    "support_mode",
    "support_values_present",
}

REPORT_SECTION_TITLES = (
    "Run Summary",
    "BUSCO Completeness",
    "Locus Retention",
    "Gene Trees",
    "Species Tree",
    "Concordance Metrics",
    "Output Index",
)


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def _require_columns(rows: list[dict[str, str]], required: set[str], label: str) -> None:
    fieldnames = set(rows[0].keys()) if rows else set()
    missing = sorted(required - fieldnames)
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def _parse_int(value: str) -> int:
    return int(value.strip())


def _parse_float(value: str) -> float:
    return float(value.strip())


def _split_csv_field(value: str) -> list[str]:
    return [item for item in value.split(",") if item]


def _format_percent(value: float, digits: int = 1) -> str:
    return f"{value:.{digits}f}%"


def _format_ratio(numerator: int, denominator: int, digits: int = 1) -> str:
    if denominator == 0:
        return "0 / 0 (0.0%)"
    percent = 100.0 * numerator / denominator
    return f"{numerator} / {denominator} ({percent:.{digits}f}%)"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def parse_species_tree_log(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("Version: "):
                metadata["version"] = line.split(": ", 1)[1]
            elif line.startswith("#Genetrees: "):
                metadata["gene_trees"] = line.split(": ", 1)[1]
            elif line.startswith("#Species: "):
                metadata["species"] = line.split(": ", 1)[1]
            elif line.startswith("#Rounds: "):
                metadata["rounds"] = line.split(": ", 1)[1]
            elif line.startswith("#Samples: "):
                metadata["subsamples"] = line.split(": ", 1)[1]
            elif line.startswith("#Threads: "):
                metadata["threads"] = line.split(": ", 1)[1]
            elif line.startswith("Initial score: "):
                metadata["initial_score"] = line.split(": ", 1)[1]
            elif line.startswith("Score: "):
                metadata["final_score"] = line.split(": ", 1)[1]
            elif line.startswith("Final Tree: "):
                metadata["final_tree"] = line.split(": ", 1)[1]
    return metadata


def summarize_busco_rows(rows: list[dict[str, str]]) -> dict[str, object]:
    if not rows:
        raise ValueError("BUSCO summary table is empty.")
    _require_columns(rows, BUSCO_REQUIRED_COLUMNS, "BUSCO summary table")
    return {
        "sample_count": len(rows),
        "busco_version": rows[0]["busco_version"],
        "lineage_name": rows[0]["lineage_name"],
        "rows": rows,
    }


def summarize_retained_loci(rows: list[dict[str, str]]) -> dict[str, object]:
    if not rows:
        raise ValueError("Retained loci table is empty.")
    _require_columns(rows, RETAINED_LOCI_REQUIRED_COLUMNS, "Retained loci table")

    retained_rows = [row for row in rows if row["decision"] == "retain"]
    excluded_rows = [row for row in rows if row["decision"] != "retain"]
    occupancy_distribution = Counter(row["occupancy"] for row in retained_rows)
    failure_counts = Counter()
    warning_counts = Counter()

    for row in rows:
        for reason in _split_csv_field(row["failure_reasons"]):
            failure_counts[reason] += 1
        for warning in _split_csv_field(row["qc_warnings"]):
            warning_counts[warning] += 1

    duplicated_loci = sum(1 for row in rows if _parse_int(row["duplicated_taxa"]) > 0)
    duplicated_retained = sum(1 for row in retained_rows if _parse_int(row["duplicated_taxa"]) > 0)
    stop_codon_loci = sum(1 for row in rows if _parse_int(row["stop_codon_taxa"]) > 0)

    return {
        "total_loci": len(rows),
        "retained_loci": len(retained_rows),
        "excluded_loci": len(excluded_rows),
        "occupancy_threshold": _parse_float(rows[0]["occupancy_threshold"]),
        "occupancy_distribution": occupancy_distribution,
        "failure_counts": failure_counts,
        "warning_counts": warning_counts,
        "duplicated_loci": duplicated_loci,
        "duplicated_retained": duplicated_retained,
        "stop_codon_loci": stop_codon_loci,
    }


def summarize_gene_trees(rows: list[dict[str, str]]) -> dict[str, object]:
    if not rows:
        raise ValueError("Gene-tree manifest is empty.")
    _require_columns(rows, GENE_TREE_REQUIRED_COLUMNS, "Gene-tree manifest")

    support_modes = Counter(row["support_mode"] for row in rows)
    support_values = Counter(row["support_values_present"] for row in rows)
    selected_models = Counter(row["selected_model"] for row in rows)

    return {
        "gene_tree_count": len(rows),
        "support_modes": support_modes,
        "support_values": support_values,
        "selected_models": selected_models,
    }


def render_report(
    *,
    busco_summary_path: Path,
    retained_loci_path: Path,
    gene_tree_manifest_path: Path,
    species_tree_path: Path,
    species_tree_log_path: Path,
    species_tree_backend: str,
    concordance_paths: list[Path] | None = None,
) -> str:
    busco_rows = read_tsv_rows(busco_summary_path)
    retained_rows = read_tsv_rows(retained_loci_path)
    gene_tree_rows = read_tsv_rows(gene_tree_manifest_path)
    species_tree_text = read_single_line_tree(species_tree_path)
    species_log = parse_species_tree_log(species_tree_log_path)

    busco_summary = summarize_busco_rows(busco_rows)
    retained_summary = summarize_retained_loci(retained_rows)
    gene_tree_summary = summarize_gene_trees(gene_tree_rows)

    valid_concordance_paths = [
        path for path in (concordance_paths or []) if path.is_file()
    ]
    gcf_stat_path = next(
        (path for path in valid_concordance_paths if path.name == "gcf.cf.stat"),
        None,
    )
    gcf_summary = summarize_gcf_stat(gcf_stat_path) if gcf_stat_path else None
    now = datetime.now().astimezone().isoformat(timespec="seconds")

    busco_table = markdown_table(
        [
            "Sample",
            "Complete",
            "Single-copy",
            "Multi-copy",
            "Fragmented",
            "Missing",
            "Internal stops",
        ],
        [
            [
                row["sample_id"],
                _format_percent(_parse_float(row["complete_percent"])),
                _format_percent(_parse_float(row["single_copy_percent"])),
                _format_percent(_parse_float(row["multi_copy_percent"])),
                _format_percent(_parse_float(row["fragmented_percent"])),
                _format_percent(_parse_float(row["missing_percent"])),
                row["internal_stop_codon_count"],
            ]
            for row in busco_summary["rows"]
        ],
    )

    occupancy_rows = [
        [_format_percent(100.0 * _parse_float(occupancy)), str(count)]
        for occupancy, count in sorted(
            retained_summary["occupancy_distribution"].items(),
            key=lambda item: _parse_float(item[0]),
            reverse=True,
        )
    ]
    failure_rows = [
        [reason, str(count)]
        for reason, count in sorted(
            retained_summary["failure_counts"].items(),
            key=lambda item: (-item[1], item[0]),
        )
    ]
    warning_rows = [
        [warning, str(count)]
        for warning, count in sorted(
            retained_summary["warning_counts"].items(),
            key=lambda item: (-item[1], item[0]),
        )[:5]
    ]
    model_rows = [
        [model, str(count)]
        for model, count in gene_tree_summary["selected_models"].most_common(5)
    ]

    lines = [
        "# Workflow Report",
        "",
        f"Generated: `{now}`",
        "",
        "## Run Summary",
        f"- Samples: `{busco_summary['sample_count']}`",
        f"- BUSCO lineage: `{busco_summary['lineage_name']}` with BUSCO `{busco_summary['busco_version']}`",
        f"- Candidate loci: `{retained_summary['total_loci']}`",
        f"- Retained loci: `{_format_ratio(retained_summary['retained_loci'], retained_summary['total_loci'])}`",
        f"- Gene trees inferred: `{gene_tree_summary['gene_tree_count']}`",
        f"- Species-tree backend: `{species_tree_backend}`",
        f"- ASTER version: `{species_log.get('version', 'unknown')}`",
        "",
        "Source artifacts: `results/qc/busco_summary.tsv`, `results/qc/retained_loci.tsv`, "
        "`results/gene_trees/gene_tree_manifest.tsv`, "
        f"`{species_tree_path.as_posix()}`, `{species_tree_log_path.as_posix()}`",
        "",
        "## BUSCO Completeness",
        busco_table,
        "",
        "## Locus Retention",
        "Only occupancy is a hard exclusion in v1. Duplicate hits, internal stops, "
        "and length-dispersion signals are carried forward as QC warnings.",
        "",
        f"- Occupancy threshold: `>= {retained_summary['occupancy_threshold']:.1f}`",
        f"- Excluded loci: `{_format_ratio(retained_summary['excluded_loci'], retained_summary['total_loci'])}`",
        f"- Loci with duplicated taxa: `{_format_ratio(retained_summary['duplicated_loci'], retained_summary['total_loci'])}`",
        f"- Retained loci with duplicated taxa: `{_format_ratio(retained_summary['duplicated_retained'], retained_summary['retained_loci'])}`",
        f"- Loci with internal stop codons: `{_format_ratio(retained_summary['stop_codon_loci'], retained_summary['total_loci'])}`",
        "",
        "Retained occupancy distribution:",
        "",
        markdown_table(["Occupancy", "Retained loci"], occupancy_rows),
        "",
        "Exclusion reasons:",
        "",
        markdown_table(["Reason", "Excluded loci"], failure_rows),
        "",
        "Top QC warnings:",
        "",
        markdown_table(["Warning", "Loci"], warning_rows),
        "",
        "## Gene Trees",
        f"- Support mode: `{next(iter(gene_tree_summary['support_modes']))}`",
        f"- Trees with parsed support values: `{_format_ratio(gene_tree_summary['support_values'].get('true', 0), gene_tree_summary['gene_tree_count'])}`",
        f"- Trees without parsed support values: `{_format_ratio(gene_tree_summary['support_values'].get('false', 0), gene_tree_summary['gene_tree_count'])}`",
        "",
        "Most common selected models:",
        "",
        markdown_table(["Model", "Loci"], model_rows),
        "",
        "## Species Tree",
        f"- Backend: `{species_tree_backend}`",
        f"- Input gene trees: `{species_log.get('gene_trees', 'unknown')}`",
        f"- Species: `{species_log.get('species', 'unknown')}`",
        f"- ASTER rounds: `{species_log.get('rounds', 'unknown')}`",
        f"- ASTER subsamples: `{species_log.get('subsamples', 'unknown')}`",
        f"- Threads: `{species_log.get('threads', 'unknown')}`",
        f"- Final score: `{species_log.get('final_score', 'unknown')}`",
        "",
        "Final tree:",
        "",
        "```newick",
        species_tree_text,
        "```",
        "",
        "## Concordance Metrics",
    ]

    if gcf_summary:
        low_rows = [
            [row["branch_id"], f"{row['gcf']:.2f}", str(row["gn"])]
            for row in gcf_summary["lowest_rows"]
        ]
        lines.extend(
            [
                f"- gCF branches scored: `{gcf_summary['branch_count']}`",
                f"- Mean gCF: `{gcf_summary['mean_gcf']:.2f}`",
                f"- Median gCF: `{gcf_summary['median_gcf']:.2f}`",
                f"- Min gCF: `{gcf_summary['min_gcf']:.2f}`",
                f"- Max gCF: `{gcf_summary['max_gcf']:.2f}`",
                "",
                "Lowest-gCF branches:",
                "",
                markdown_table(["Branch ID", "gCF", "Decisive gene trees"], low_rows),
                "",
            ]
        )

    if valid_concordance_paths:
        lines.append("Available concordance artifacts:")
        lines.append("")
        for path in valid_concordance_paths:
            lines.append(f"- `{path.as_posix()}`")
    else:
        lines.append(
            "Not implemented in v1. No concordance or quartet-score outputs were generated for this run."
        )

    lines.extend(
        [
            "",
            "## Output Index",
            markdown_table(
                ["Artifact", "Path"],
                [
                    ["BUSCO summary table", busco_summary_path.as_posix()],
                    ["Retained loci table", retained_loci_path.as_posix()],
                    ["Gene-tree manifest", gene_tree_manifest_path.as_posix()],
                    ["Species tree", species_tree_path.as_posix()],
                    ["Species-tree log", species_tree_log_path.as_posix()],
                    *[
                        [f"Concordance artifact {index}", path.as_posix()]
                        for index, path in enumerate(valid_concordance_paths, start=1)
                    ],
                ],
            ),
            "",
        ]
    )
    return "\n".join(lines)


def write_report(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
