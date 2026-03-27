"""Verify that the configured BUSCO lineage exists in the BUSCO dataset list."""

from __future__ import annotations

import argparse
from pathlib import Path

from .busco import load_dataset_names, verify_lineage_name
from .manifest import write_tsv


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-list", required=True, help="Path to BUSCO dataset list output.")
    parser.add_argument("--lineage", required=True, help="Configured BUSCO lineage string.")
    parser.add_argument("--output", required=True, help="Verification TSV output path.")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    dataset_list_path = Path(args.dataset_list)
    output_path = Path(args.output)
    available_datasets = load_dataset_names(dataset_list_path)
    lineage = verify_lineage_name(args.lineage, available_datasets)

    write_tsv(
        output_path,
        [
            {
                "configured_lineage": lineage,
                "matched_lineage": lineage,
                "available_dataset_count": str(len(available_datasets)),
            }
        ],
        ["configured_lineage", "matched_lineage", "available_dataset_count"],
    )


if __name__ == "__main__":
    main()

