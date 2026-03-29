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

SCFL_REQUIRED_COLUMNS = {
    "ID",
    "sCF",
    "sCF_N",
    "sDF1",
    "sDF1_N",
    "sDF2",
    "sDF2_N",
    "sN",
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


def build_iqtree_scfl_command(
    *,
    executable: str,
    reference_tree_path: str,
    alignment_dir: str,
    prefix: str,
    threads: int,
    quartets: int,
    seqtype: str,
    model: str,
) -> list[str]:
    return [
        executable,
        "-te",
        reference_tree_path,
        "-p",
        alignment_dir,
        "--scfl",
        str(quartets),
        "--seqtype",
        seqtype,
        "-m",
        model,
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


def _summarize_cf_stat(
    *,
    path: Path,
    required_columns: set[str],
    value_column: str,
    count_column: str,
    value_output_key: str,
    count_output_key: str,
    label: str,
    count_parser,
) -> dict[str, object]:
    rows = read_cf_stat_rows(path)
    fieldnames = set(rows[0].keys()) if rows else set()
    missing = sorted(required_columns - fieldnames)
    if missing:
        raise ValueError(f"{label} statistics file is missing required columns: {', '.join(missing)}")

    scored_rows = [
        {
            "branch_id": row["ID"],
            value_output_key: _parse_optional_float(row[value_column]),
            count_output_key: count_parser(row[count_column]),
        }
        for row in rows
        if _parse_optional_float(row[value_column]) is not None and count_parser(row[count_column]) is not None
    ]
    if not scored_rows:
        raise ValueError(f"No scored {label} branches were parsed from {path}")
    values = [row[value_output_key] for row in scored_rows]
    lowest_rows = sorted(scored_rows, key=lambda row: (row[value_output_key], row["branch_id"]))[:5]

    return {
        "branch_count": len(scored_rows),
        f"mean_{value_output_key}": statistics.fmean(values),
        f"median_{value_output_key}": statistics.median(values),
        f"min_{value_output_key}": min(values),
        f"max_{value_output_key}": max(values),
        "lowest_rows": lowest_rows,
    }


def summarize_gcf_stat(path: Path) -> dict[str, object]:
    return _summarize_cf_stat(
        path=path,
        required_columns=GCF_REQUIRED_COLUMNS,
        value_column="gCF",
        count_column="gN",
        value_output_key="gcf",
        count_output_key="gn",
        label="gCF",
        count_parser=_parse_optional_int,
    )


def summarize_scfl_stat(path: Path) -> dict[str, object]:
    return _summarize_cf_stat(
        path=path,
        required_columns=SCFL_REQUIRED_COLUMNS,
        value_column="sCF",
        count_column="sN",
        value_output_key="scfl",
        count_output_key="sn",
        label="sCFL",
        count_parser=_parse_optional_float,
    )
