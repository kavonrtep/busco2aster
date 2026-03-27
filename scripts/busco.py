"""Helpers for BUSCO environment and lineage verification."""

from __future__ import annotations

import re
import shutil
from pathlib import Path

from .manifest import write_tsv


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


def busco_output_paths(sample_id: str) -> dict[str, str]:
    sample_root = Path("results") / "busco" / sample_id
    return {
        "sample_root": sample_root.as_posix(),
        "raw_root": (sample_root / "raw").as_posix(),
        "command": (sample_root / "command.sh").as_posix(),
        "paths": (sample_root / "paths.tsv").as_posix(),
        "short_summary": (sample_root / "short_summary.txt").as_posix(),
        "full_table": (sample_root / "full_table.tsv").as_posix(),
        "completion": (sample_root / "run.complete").as_posix(),
    }


def _find_single_path(base_dir: Path, pattern: str) -> Path:
    matches = sorted(base_dir.rglob(pattern))
    if not matches:
        raise BuscoValidationError(f"Expected BUSCO artifact {pattern!r} under {base_dir}")
    return matches[0]


def _find_optional_dir(base_dir: Path, name: str) -> str:
    matches = sorted(path for path in base_dir.rglob(name) if path.is_dir())
    return matches[0].as_posix() if matches else ""


def standardize_busco_run(sample_id: str, sample_dir: Path, raw_root: Path) -> None:
    short_summary_source = _find_single_path(raw_root, "short_summary*")
    full_table_source = _find_single_path(raw_root, "full_table.tsv")

    stable_paths = busco_output_paths(sample_id)
    short_summary_target = Path(stable_paths["short_summary"])
    full_table_target = Path(stable_paths["full_table"])
    paths_target = Path(stable_paths["paths"])
    completion_target = Path(stable_paths["completion"])

    short_summary_target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(short_summary_source, short_summary_target)
    shutil.copy2(full_table_source, full_table_target)

    write_tsv(
        paths_target,
        [
            {"artifact": "raw_root", "path": raw_root.as_posix()},
            {"artifact": "short_summary_source", "path": short_summary_source.as_posix()},
            {"artifact": "full_table_source", "path": full_table_source.as_posix()},
            {
                "artifact": "single_copy_busco_sequences",
                "path": _find_optional_dir(raw_root, "single_copy_busco_sequences"),
            },
            {
                "artifact": "multi_copy_busco_sequences",
                "path": _find_optional_dir(raw_root, "multi_copy_busco_sequences"),
            },
            {
                "artifact": "fragmented_busco_sequences",
                "path": _find_optional_dir(raw_root, "fragmented_busco_sequences"),
            },
        ],
        ["artifact", "path"],
    )
    completion_target.write_text(
        f"sample_id\t{sample_id}\nraw_root\t{raw_root.as_posix()}\n",
        encoding="utf-8",
    )
