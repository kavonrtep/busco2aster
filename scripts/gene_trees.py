"""Helpers for IQ-TREE 3 gene-tree inference and aggregation."""

from __future__ import annotations

import re
import shlex
from pathlib import Path

from .alignment import busco_sort_key, load_retained_locus_ids
from .manifest import write_tsv
from .sequence_mode import locus_id_from_alignment_filename


SUPPORTED_SUPPORT_MODES = {"abayes", "ufboot"}
BEST_FIT_MODEL_PATTERN = re.compile(r"^Best-fit model according to BIC:\s+(.+?)\s*$")
MODEL_OF_SUBSTITUTION_PATTERN = re.compile(r"^Model of substitution:\s+(.+?)\s*$")
ABAYES_SUPPORT_PATTERN = re.compile(r"\)/([0-9]+(?:\.[0-9]+)?):")
UFBOOT_SUPPORT_PATTERN = re.compile(r"\)([0-9]+(?:\.[0-9]+)?):")
DIRECTORY_MODE_MODEL_LINE_PATTERN = re.compile(r"^\s*(.+?):\s+(.+?)\{[^{}]*\},?\s*$")


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


def gene_tree_directory_output_paths() -> dict[str, str]:
    prefix = Path("results") / "gene_trees" / "gene_trees"
    return {
        "prefix": prefix.as_posix(),
        "command": Path("results/gene_trees/gene_trees.command.sh").as_posix(),
        "report": Path(f"{prefix}.iqtree").as_posix(),
        "log": Path(f"{prefix}.log").as_posix(),
        "treefile": Path(f"{prefix}.treefile").as_posix(),
        "best_model_nex": Path(f"{prefix}.best_model.nex").as_posix(),
        "best_scheme": Path(f"{prefix}.best_scheme").as_posix(),
        "best_scheme_nex": Path(f"{prefix}.best_scheme.nex").as_posix(),
        "model_gz": Path(f"{prefix}.model.gz").as_posix(),
        "checkpoint": Path(f"{prefix}.ckp.gz").as_posix(),
        "parstree": Path(f"{prefix}.parstree").as_posix(),
        "completion": Path("results/gene_trees/gene_trees.complete").as_posix(),
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
    seqtype: str = "AA",
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
        seqtype,
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


def build_iqtree_directory_command(
    executable: str,
    alignment_dir: str,
    prefix: str,
    threads: int,
    model: str,
    support_mode: str,
    seed: int,
    ufboot_replicates: int | None = None,
    seqtype: str = "AA",
) -> list[str]:
    if support_mode not in SUPPORTED_SUPPORT_MODES:
        raise ValueError(
            f"Unsupported IQ-TREE support mode {support_mode!r}; "
            f"expected one of {sorted(SUPPORTED_SUPPORT_MODES)}."
        )

    command = [
        executable,
        "-S",
        alignment_dir,
        "--seqtype",
        seqtype,
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


def _alignment_name_to_locus_id(alignment_name: str) -> str:
    return locus_id_from_alignment_filename(alignment_name)


def parse_directory_mode_selected_models(best_model_nex_path: Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    in_charpartition = False
    with best_model_nex_path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("charpartition mymodels"):
                in_charpartition = True
                line = line.split("=", 1)[1].strip()
                if not line:
                    continue
            if not in_charpartition:
                continue
            if line == "end;":
                break
            terminal = line.endswith(";")
            normalized = line.rstrip(",;").strip()
            if not normalized:
                if terminal:
                    break
                continue
            match = DIRECTORY_MODE_MODEL_LINE_PATTERN.match(normalized)
            if not match:
                raise ValueError(f"Could not parse directory-mode model line: {raw_line.rstrip()}")
            model_text = re.sub(r"\{[^{}]*\}$", "", match.group(1).strip())
            alignment_name = match.group(2).strip()
            entries.append((_alignment_name_to_locus_id(alignment_name), model_text))
            if terminal:
                break
    if not entries:
        raise ValueError(f"No per-locus models were parsed from {best_model_nex_path}")
    return entries


def aggregate_directory_mode_gene_tree_outputs(
    *,
    best_model_nex_path: Path,
    treefile_path: Path,
    report_path: Path,
    manifest_path: Path,
    aggregate_path: Path,
    support_mode: str,
) -> None:
    model_entries = parse_directory_mode_selected_models(best_model_nex_path)
    tree_lines = [line.strip() for line in treefile_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(model_entries) != len(tree_lines):
        raise ValueError(
            f"Directory-mode tree count ({len(tree_lines)}) does not match model count ({len(model_entries)})."
        )

    aggregate_path.parent.mkdir(parents=True, exist_ok=True)
    aggregate_path.write_text("\n".join(tree_lines) + "\n", encoding="utf-8")

    manifest_rows = []
    for index, ((locus_id, selected_model), tree_text) in enumerate(zip(model_entries, tree_lines, strict=True), start=1):
        manifest_rows.append(
            {
                "locus_id": locus_id,
                "treefile": aggregate_path.as_posix(),
                "tree_row_index": str(index),
                "report": report_path.as_posix(),
                "selected_model": selected_model,
                "support_mode": support_mode,
                "support_values_present": str(tree_has_support_values(tree_text, support_mode)).lower(),
            }
        )

    write_tsv(
        manifest_path,
        manifest_rows,
        [
            "locus_id",
            "treefile",
            "tree_row_index",
            "report",
            "selected_model",
            "support_mode",
            "support_values_present",
        ],
    )
