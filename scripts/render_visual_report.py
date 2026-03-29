"""Render the Quarto visual report."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quarto-executable", default="quarto")
    parser.add_argument("--qmd", required=True)
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    qmd_path = Path(args.qmd).resolve()
    data_dir = Path(args.data_dir).resolve()
    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    repo_root = qmd_path.parent.parent.resolve()
    try:
        qmd_arg = qmd_path.relative_to(repo_root).as_posix()
        output_dir_arg = output_path.parent.relative_to(repo_root).as_posix()
    except ValueError:
        qmd_arg = qmd_path.as_posix()
        output_dir_arg = output_path.parent.as_posix()

    command = [
        args.quarto_executable,
        "render",
        qmd_arg,
        "--to",
        "html",
        "--output",
        output_path.name,
        "--output-dir",
        output_dir_arg,
        "-P",
        f"data_dir:{data_dir.as_posix()}",
    ]
    subprocess.run(command, check=True, cwd=repo_root)
    if output_path.is_file():
        return 0

    fallback_candidates = [
        repo_root / "results" / output_path.name,
        repo_root / output_path.name,
    ]
    for candidate in fallback_candidates:
        if candidate.is_file():
            output_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(candidate.as_posix(), output_path.as_posix())
            return 0

    raise FileNotFoundError(f"Quarto render completed but output was not found at {output_path}.")


if __name__ == "__main__":
    raise SystemExit(main())
