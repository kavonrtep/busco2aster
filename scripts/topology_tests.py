"""Helpers for alternative topology hypothesis testing."""

from __future__ import annotations

import shlex
from pathlib import Path


def topology_tests_output_paths() -> dict[str, str]:
    prefix = Path("results") / "topology_tests"
    return {
        "wastral_u2_tree": (prefix / "species_tree.wastral.u2.tre").as_posix(),
        "wastral_u2_command": (prefix / "wastral_u2.command.sh").as_posix(),
        "wastral_u2_log": (prefix / "wastral_u2.log").as_posix(),
        "branch_quartet_support": (prefix / "branch_quartet_support.tsv").as_posix(),
        "contested_branches": (prefix / "contested_branches.tsv").as_posix(),
        "candidate_trees": (prefix / "candidate_trees.tre").as_posix(),
        "candidate_manifest": (prefix / "candidate_trees_manifest.tsv").as_posix(),
        "supermatrix": (prefix / "supermatrix.phy").as_posix(),
        "partitions": (prefix / "supermatrix_partitions.nex").as_posix(),
        "au_prefix": (prefix / "au_test").as_posix(),
        "au_iqtree": (prefix / "au_test.iqtree").as_posix(),
        "au_results": (prefix / "au_test_results.tsv").as_posix(),
        "au_command": (prefix / "au_test.command.sh").as_posix(),
        "au_log": (prefix / "au_test.log").as_posix(),
        "completion": (prefix / "topology_tests.complete").as_posix(),
    }


def build_au_test_command(
    *,
    executable: str,
    supermatrix: str,
    partitions: str,
    candidate_trees: str,
    prefix: str,
    threads: int,
    replicates: int,
) -> list[str]:
    return [
        executable,
        "-s", supermatrix,
        "-p", partitions,
        "--trees", candidate_trees,
        "--test-au",
        "-n", "0",
        "--test", str(replicates),
        "-T", str(threads),
        "--prefix", prefix,
        "--redo",
    ]


def write_topology_tests_command_script(
    path: Path,
    command: list[str],
    stderr_path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rendered = " ".join(shlex.quote(token) for token in command)
    rendered = f"{rendered} 2> {shlex.quote(Path(stderr_path).as_posix())}"
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", rendered]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    path.chmod(0o755)
