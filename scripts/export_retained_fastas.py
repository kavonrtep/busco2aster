"""Export all retained-locus protein FASTAs in one batch."""

from __future__ import annotations

import argparse
from pathlib import Path

from .alignment import export_retained_fastas


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, help="Path to results/qc/locus_taxon_matrix.tsv.")
    parser.add_argument("--retained", required=True, help="Path to results/qc/retained_loci.tsv.")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where one FASTA per retained locus will be written.",
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="TSV manifest listing the retained-locus FASTA files generated in this batch.",
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve relative source FASTA paths.",
    )
    parser.add_argument(
        "--sequence-type",
        default="protein",
        help="Sequence type to export: protein or dna. Default: protein.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    export_retained_fastas(
        matrix_path=Path(args.matrix),
        retained_path=Path(args.retained),
        output_dir=Path(args.output_dir),
        manifest_path=Path(args.manifest),
        repo_root=Path(args.repo_root).resolve(),
        sequence_type=args.sequence_type,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
