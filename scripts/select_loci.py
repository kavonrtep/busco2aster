"""Select retained loci from the long-form locus matrix."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from .locus_matrix import build_retained_loci_rows, write_retained_loci


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, help="Path to results/qc/locus_taxon_matrix.tsv.")
    parser.add_argument("--output", required=True, help="Path to results/qc/retained_loci.tsv.")
    parser.add_argument(
        "--occupancy-threshold",
        required=True,
        type=float,
        help="Minimum fraction of taxa with complete single-copy hits to retain a locus.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    matrix_path = Path(args.matrix)
    with matrix_path.open(newline="", encoding="utf-8") as handle:
        matrix_rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]

    retained_rows = build_retained_loci_rows(matrix_rows, args.occupancy_threshold)
    write_retained_loci(Path(args.output), retained_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
