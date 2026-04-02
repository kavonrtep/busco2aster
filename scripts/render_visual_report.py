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

    command = [
        args.quarto_executable,
        "render",
        qmd_path.as_posix(),
        "--to",
        "html",
        "--output",
        output_path.name,
        "-P",
        f"data_dir:{data_dir.as_posix()}",
    ]
    subprocess.run(command, check=True, cwd=output_path.parent)

    source_sidecar = qmd_path.parent / f"{qmd_path.stem}_files"
    target_sidecar = output_path.parent / f"{output_path.stem}_files"
    if source_sidecar.is_dir():
        if target_sidecar.exists():
            shutil.rmtree(target_sidecar)
        shutil.copytree(source_sidecar, target_sidecar)
        if source_sidecar.resolve() != target_sidecar.resolve():
            shutil.rmtree(source_sidecar)

    if output_path.is_file():
        return 0

    fallback_candidates = [
        output_path.parent / output_path.name,
        output_path.parent.parent / output_path.name,
        qmd_path.parent / output_path.name,
    ]
    for candidate in fallback_candidates:
        if candidate.is_file():
            output_path.parent.mkdir(parents=True, exist_ok=True)
            if candidate.resolve() != output_path:
                shutil.move(candidate.as_posix(), output_path.as_posix())
            return 0

    raise FileNotFoundError(f"Quarto render completed but output was not found at {output_path}.")


if __name__ == "__main__":
    raise SystemExit(main())
