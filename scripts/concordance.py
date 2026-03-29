"""Helpers for branch concordance analyses."""

from __future__ import annotations

import csv
import shlex
import statistics
from pathlib import Path


GCF_REQUIRED_COLUMNS = {
    "ID",
    "gCF",
    "gCF_N",
    "gDF1",
    "gDF1_N",
    "gDF2",
    "gDF2_N",
    "gDFP",
    "gDFP_N",
    "gN",
    "Label",
    "Length",
}


def concordance_output_paths(name: str) -> dict[str, str]:
    prefix = Path("results") / "concordance" / name
    return {
        "name": name,
        "prefix": prefix.as_posix(),
        "command": Path(f"{prefix}.command.sh").as_posix(),
        "stat": Path(f"{prefix}.cf.stat").as_posix(),
        "tree": Path(f"{prefix}.cf.tree").as_posix(),
        "branch": Path(f"{prefix}.cf.branch").as_posix(),
        "log": Path(f"{prefix}.log").as_posix(),
        "completion": Path(f"{prefix}.complete").as_posix(),
    }


def build_iqtree_gcf_command(
    *,
    executable: str,
    reference_tree_path: str,
    gene_tree_path: str,
    prefix: str,
    threads: int,
) -> list[str]:
    return [
        executable,
        "-t",
        reference_tree_path,
        "--gcf",
        gene_tree_path,
        "--prefix",
        prefix,
        "-T",
        str(threads),
        "--redo",
        "--quiet",
    ]


def write_concordance_command_script(path: Path, command: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rendered = " ".join(shlex.quote(token) for token in command)
    path.write_text(f"#!/usr/bin/env bash\nset -euo pipefail\n{rendered}\n", encoding="utf-8")
    path.chmod(0o755)


def read_cf_stat_rows(path: Path) -> list[dict[str, str]]:
    lines = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            lines.append(line)
    if not lines:
        raise ValueError(f"Concordance statistics file is empty: {path}")
    reader = csv.DictReader(lines, delimiter="\t")
    return [dict(row) for row in reader]


def _parse_optional_float(value: str) -> float | None:
    text = value.strip()
    if text in {"", "NA", "N/A", "nan", "NaN"}:
        return None
    return float(text)


def _parse_optional_int(value: str) -> int | None:
    text = value.strip()
    if text in {"", "NA", "N/A"}:
        return None
    return int(text)


def summarize_gcf_stat(path: Path) -> dict[str, object]:
    rows = read_cf_stat_rows(path)
    fieldnames = set(rows[0].keys()) if rows else set()
    missing = sorted(GCF_REQUIRED_COLUMNS - fieldnames)
    if missing:
        raise ValueError(
            f"gCF statistics file is missing required columns: {', '.join(missing)}"
        )

    scored_rows = [
        {
            "branch_id": row["ID"],
            "gcf": _parse_optional_float(row["gCF"]),
            "gcf_n": _parse_optional_int(row["gCF_N"]),
            "gn": _parse_optional_int(row["gN"]),
        }
        for row in rows
        if _parse_optional_float(row["gCF"]) is not None and _parse_optional_int(row["gN"]) is not None
    ]
    if not scored_rows:
        raise ValueError(f"No scored gCF branches were parsed from {path}")
    gcf_values = [row["gcf"] for row in scored_rows]
    lowest_rows = sorted(scored_rows, key=lambda row: (row["gcf"], row["branch_id"]))[:5]

    return {
        "branch_count": len(scored_rows),
        "mean_gcf": statistics.fmean(gcf_values),
        "median_gcf": statistics.median(gcf_values),
        "min_gcf": min(gcf_values),
        "max_gcf": max(gcf_values),
        "lowest_rows": lowest_rows,
    }
