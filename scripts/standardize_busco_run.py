"""Create stable per-sample BUSCO outputs after a BUSCO run completes."""

from __future__ import annotations

import argparse
from pathlib import Path

from .busco import standardize_busco_run


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True, help="Normalized sample ID.")
    parser.add_argument("--sample-dir", required=True, help="Stable sample output directory.")
    parser.add_argument("--raw-root", required=True, help="BUSCO raw output directory.")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    standardize_busco_run(
        sample_id=args.sample_id,
        sample_dir=Path(args.sample_dir),
        raw_root=Path(args.raw_root),
    )


if __name__ == "__main__":
    main()
