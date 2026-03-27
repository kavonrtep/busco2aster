"""Helpers for ASTER species-tree inference."""

from __future__ import annotations

import re
import shlex
from pathlib import Path


SUPPORTED_SPECIES_TREE_BACKENDS = {"astral4", "wastral"}
SUPPORTED_WASTRAL_SUPPORT_MODES = {"abayes", "ufboot"}
ABAYES_INTERNAL_LABEL_PATTERN = re.compile(r"\)/([0-9]+(?:\.[0-9]+)?):")


def species_tree_output_paths(backend: str) -> dict[str, str]:
    if backend not in SUPPORTED_SPECIES_TREE_BACKENDS:
        raise ValueError(
            f"Unsupported species-tree backend {backend!r}; "
            f"expected one of {sorted(SUPPORTED_SPECIES_TREE_BACKENDS)}."
        )

    prefix = Path("results") / "species_tree" / f"species_tree.{backend}"
    return {
        "backend": backend,
        "prefix": prefix.as_posix(),
        "command": Path(f"{prefix}.command.sh").as_posix(),
        "treefile": Path(f"{prefix}.tre").as_posix(),
        "log": Path(f"{prefix}.log").as_posix(),
        "completion": Path(f"{prefix}.complete").as_posix(),
    }


def build_aster_command(
    backend: str,
    executable: str,
    input_path: str,
    output_path: str,
    threads: int,
    support_mode: str | None = None,
    mapping_path: str | None = None,
    root_outgroup: str | None = None,
    annotation_mode: int | None = None,
) -> list[str]:
    if backend not in SUPPORTED_SPECIES_TREE_BACKENDS:
        raise ValueError(
            f"Unsupported species-tree backend {backend!r}; "
            f"expected one of {sorted(SUPPORTED_SPECIES_TREE_BACKENDS)}."
        )

    command = [
        executable,
        "-t",
        str(threads),
        "-o",
        output_path,
        "-i",
        input_path,
    ]
    if backend == "wastral":
        if support_mode not in SUPPORTED_WASTRAL_SUPPORT_MODES:
            raise ValueError(
                f"Unsupported wASTRAL support mode {support_mode!r}; "
                f"expected one of {sorted(SUPPORTED_WASTRAL_SUPPORT_MODES)}."
            )
        if support_mode == "abayes":
            command.append("-B")
        else:
            command.append("-S")
    if mapping_path:
        command.extend(["-a", mapping_path])
    if root_outgroup:
        command.extend(["--root", root_outgroup])
    if annotation_mode is not None:
        command.extend(["-u", str(annotation_mode)])
    return command


def write_species_tree_command_script(path: Path, command: list[str], stderr_path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rendered = " ".join(shlex.quote(token) for token in command)
    rendered = f"{rendered} 2> {shlex.quote(stderr_path.as_posix())}"
    path.write_text(f"#!/usr/bin/env bash\nset -euo pipefail\n{rendered}\n", encoding="utf-8")
    path.chmod(0o755)


def normalize_gene_tree_for_wastral(tree_text: str, support_mode: str) -> str:
    if support_mode == "abayes":
        return ABAYES_INTERNAL_LABEL_PATTERN.sub(r")\1:", tree_text)
    if support_mode == "ufboot":
        return tree_text
    raise ValueError(
        f"Unsupported wASTRAL support mode {support_mode!r}; "
        f"expected one of {sorted(SUPPORTED_WASTRAL_SUPPORT_MODES)}."
    )


def prepare_wastral_gene_tree_input(input_path: Path, output_path: Path, support_mode: str) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    normalized_lines = []
    with input_path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            normalized_lines.append(normalize_gene_tree_for_wastral(stripped, support_mode))
    output_path.write_text("\n".join(normalized_lines) + "\n", encoding="utf-8")
