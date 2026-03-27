"""Build the long-form BUSCO locus-by-taxon matrix."""

from __future__ import annotations

import argparse
from pathlib import Path

from .locus_matrix import (
    build_locus_taxon_matrix_rows,
    load_busco_record_rows,
    load_busco_summary_rows,
    write_locus_taxon_matrix,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, help="Path to results/qc/busco_summary.tsv.")
    parser.add_argument("--records", required=True, help="Path to results/qc/busco_records.tsv.")
    parser.add_argument("--output", required=True, help="Path to results/qc/locus_taxon_matrix.tsv.")
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve relative BUSCO sequence paths.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    summary_rows = load_busco_summary_rows(Path(args.summary))
    record_rows = load_busco_record_rows(Path(args.records))
    matrix_rows = build_locus_taxon_matrix_rows(summary_rows, record_rows, Path(args.repo_root).resolve())
    write_locus_taxon_matrix(Path(args.output), matrix_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
