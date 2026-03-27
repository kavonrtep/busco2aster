"""Convert the reference sample sheet into the normalized internal TSV."""

from __future__ import annotations

import argparse
from pathlib import Path

from .manifest import (
    REQUIRED_SAMPLE_COLUMNS,
    load_reference_manifest,
    normalize_reference_manifest_rows,
    write_tsv,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Reference sample sheet path.")
    parser.add_argument("--output", required=True, help="Normalized TSV output path.")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    rows = load_reference_manifest(input_path)
    normalized_rows = normalize_reference_manifest_rows(rows)
    write_tsv(output_path, normalized_rows, list(REQUIRED_SAMPLE_COLUMNS))


if __name__ == "__main__":
    main()

