"""Render the Quarto visual report."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
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

    # Quarto writes scratch files (and may rewrite the source .qmd) next to
    # the input. Inside the Apptainer image the pipeline tree is read-only,
    # so render from a writable scratch copy instead.
    with tempfile.TemporaryDirectory(prefix="quarto_render_", dir=output_path.parent) as tmpdir:
        scratch = Path(tmpdir)
        scratch_qmd = scratch / qmd_path.name
        shutil.copy2(qmd_path, scratch_qmd)

        command = [
            args.quarto_executable,
            "render",
            scratch_qmd.name,
            "--to",
            "html",
            "--output",
            output_path.name,
            "-P",
            f"data_dir:{data_dir.as_posix()}",
        ]
        subprocess.run(command, check=True, cwd=scratch)

        source_sidecar = scratch / f"{scratch_qmd.stem}_files"
        target_sidecar = output_path.parent / f"{output_path.stem}_files"
        if source_sidecar.is_dir():
            if target_sidecar.exists():
                shutil.rmtree(target_sidecar)
            shutil.move(source_sidecar.as_posix(), target_sidecar.as_posix())

        rendered_html = scratch / output_path.name
        if rendered_html.is_file():
            shutil.move(rendered_html.as_posix(), output_path.as_posix())
            return 0

        # Fallback: some Quarto versions emit into the working directory with
        # the stem matching the input filename rather than --output.
        fallback = scratch / f"{scratch_qmd.stem}.html"
        if fallback.is_file():
            shutil.move(fallback.as_posix(), output_path.as_posix())
            return 0

    raise FileNotFoundError(f"Quarto render completed but output was not found at {output_path}.")


if __name__ == "__main__":
    raise SystemExit(main())
