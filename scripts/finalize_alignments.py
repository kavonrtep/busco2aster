"""Finalize batched alignment outputs by removing stale files."""

from __future__ import annotations

import argparse
from pathlib import Path

from .alignment import sync_alignment_outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-dir", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sync_alignment_outputs(
        manifest_path=Path(args.manifest),
        output_dir=Path(args.output_dir),
        log_dir=Path(args.log_dir),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
