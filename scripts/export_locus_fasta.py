"""Export one retained-locus protein FASTA with sanitized taxon headers."""

from __future__ import annotations

import argparse
from pathlib import Path

from .alignment import export_locus_fasta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--locus-id", required=True, help="Retained locus identifier to export.")
    parser.add_argument("--matrix", required=True, help="Path to results/qc/locus_taxon_matrix.tsv.")
    parser.add_argument("--retained", required=True, help="Path to results/qc/retained_loci.tsv.")
    parser.add_argument("--output", required=True, help="Path to output locus FASTA.")
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve relative source FASTA paths.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    export_locus_fasta(
        locus_id=args.locus_id,
        matrix_path=Path(args.matrix),
        retained_path=Path(args.retained),
        output_path=Path(args.output),
        repo_root=Path(args.repo_root).resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
