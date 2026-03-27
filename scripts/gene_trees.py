"""Helpers for IQ-TREE 3 per-locus gene-tree inference and aggregation."""

from __future__ import annotations

import re
import shlex
from pathlib import Path

from .alignment import busco_sort_key, load_retained_locus_ids
from .manifest import write_tsv


SUPPORTED_SUPPORT_MODES = {"abayes", "ufboot"}
BEST_FIT_MODEL_PATTERN = re.compile(r"^Best-fit model according to BIC:\s+(.+?)\s*$")
MODEL_OF_SUBSTITUTION_PATTERN = re.compile(r"^Model of substitution:\s+(.+?)\s*$")
ABAYES_SUPPORT_PATTERN = re.compile(r"\)/([0-9]+(?:\.[0-9]+)?):")
UFBOOT_SUPPORT_PATTERN = re.compile(r"\)([0-9]+(?:\.[0-9]+)?):")


def gene_tree_output_paths(locus_id: str) -> dict[str, str]:
    locus_dir = Path("results") / "gene_trees" / "per_locus" / locus_id
    prefix = locus_dir / locus_id
    return {
        "dir": locus_dir.as_posix(),
        "prefix": prefix.as_posix(),
        "command": (locus_dir / "command.sh").as_posix(),
        "report": Path(f"{prefix}.iqtree").as_posix(),
        "log": Path(f"{prefix}.log").as_posix(),
        "treefile": Path(f"{prefix}.treefile").as_posix(),
        "completion": (locus_dir / "run.complete").as_posix(),
    }


def build_iqtree_command(
    executable: str,
    alignment_path: str,
    prefix: str,
    threads: int,
    model: str,
    support_mode: str,
    seed: int,
    ufboot_replicates: int | None = None,
) -> list[str]:
    if support_mode not in SUPPORTED_SUPPORT_MODES:
        raise ValueError(
            f"Unsupported IQ-TREE support mode {support_mode!r}; "
            f"expected one of {sorted(SUPPORTED_SUPPORT_MODES)}."
        )

    command = [
        executable,
        "-s",
        alignment_path,
        "--seqtype",
        "AA",
        "-m",
        model,
        "-T",
        str(threads),
        "--prefix",
        prefix,
        "--seed",
        str(seed),
        "--redo",
        "--quiet",
    ]
    if support_mode == "abayes":
        command.append("--abayes")
    else:
        if ufboot_replicates is None:
            raise ValueError("ufboot support mode requires a replicate count.")
        command.extend(["-B", str(ufboot_replicates)])
    return command


def write_command_script(path: Path, command: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rendered = " ".join(shlex.quote(token) for token in command)
    path.write_text(f"#!/usr/bin/env bash\nset -euo pipefail\n{rendered}\n", encoding="utf-8")
    path.chmod(0o755)


def parse_selected_model(report_path: Path) -> str:
    selected_model: str | None = None
    fallback_model: str | None = None
    with report_path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            best_fit_match = BEST_FIT_MODEL_PATTERN.match(stripped)
            if best_fit_match:
                selected_model = best_fit_match.group(1)
                continue
            substitution_match = MODEL_OF_SUBSTITUTION_PATTERN.match(stripped)
            if substitution_match:
                fallback_model = substitution_match.group(1)
    if selected_model:
        return selected_model
    if fallback_model:
        return fallback_model
    raise ValueError(f"Could not parse selected model from IQ-TREE report: {report_path}")


def read_single_line_tree(path: Path) -> str:
    with path.open(encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    if len(lines) != 1:
        raise ValueError(f"Expected exactly one Newick tree in {path}, found {len(lines)}.")
    return lines[0]


def tree_has_support_values(tree_text: str, support_mode: str) -> bool:
    if support_mode == "abayes":
        return bool(ABAYES_SUPPORT_PATTERN.search(tree_text))
    if support_mode == "ufboot":
        return bool(UFBOOT_SUPPORT_PATTERN.search(tree_text))
    raise ValueError(f"Unsupported support mode: {support_mode}")


def _resolve_repo_path(repo_root: Path, path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else repo_root / path


def aggregate_gene_tree_outputs(
    retained_path: Path,
    manifest_path: Path,
    aggregate_path: Path,
    support_mode: str,
    repo_root: Path | None = None,
) -> None:
    resolved_repo_root = (repo_root or Path(".")).resolve()
    locus_ids = load_retained_locus_ids(retained_path)
    manifest_rows = []
    aggregate_lines = []

    for locus_id in locus_ids:
        paths = gene_tree_output_paths(locus_id)
        tree_path = _resolve_repo_path(resolved_repo_root, paths["treefile"])
        report_path = _resolve_repo_path(resolved_repo_root, paths["report"])
        tree_text = read_single_line_tree(tree_path)
        selected_model = parse_selected_model(report_path)
        manifest_rows.append(
            {
                "locus_id": locus_id,
                "treefile": paths["treefile"],
                "report": paths["report"],
                "selected_model": selected_model,
                "support_mode": support_mode,
                "support_values_present": str(tree_has_support_values(tree_text, support_mode)).lower(),
            }
        )
        aggregate_lines.append(tree_text)

    write_tsv(
        manifest_path,
        manifest_rows,
        [
            "locus_id",
            "treefile",
            "report",
            "selected_model",
            "support_mode",
            "support_values_present",
        ],
    )
    aggregate_path.parent.mkdir(parents=True, exist_ok=True)
    aggregate_path.write_text("\n".join(aggregate_lines) + "\n", encoding="utf-8")


def load_gene_tree_locus_ids(path: Path) -> list[str]:
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        import csv

        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    return [row["locus_id"] for row in rows]
