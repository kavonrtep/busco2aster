"""Extract retained BUSCO CDS sequences for one sample with gffread."""

from __future__ import annotations

import argparse
from pathlib import Path

from .dna_extract import extract_retained_sample_dna


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, help="Path to results/qc/locus_taxon_matrix.tsv.")
    parser.add_argument("--retained", required=True, help="Path to results/qc/retained_loci.tsv.")
    parser.add_argument("--sample-id", required=True, help="Validated sample ID.")
    parser.add_argument("--genome-fasta", required=True, help="Prepared plain assembly FASTA for gffread.")
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve relative BUSCO GFF paths.",
    )
    parser.add_argument(
        "--gffread-executable",
        default="gffread",
        help="gffread executable to use. Default: gffread.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    extract_retained_sample_dna(
        matrix_path=Path(args.matrix),
        retained_path=Path(args.retained),
        sample_id=args.sample_id,
        genome_fasta_path=Path(args.genome_fasta),
        repo_root=Path(args.repo_root).resolve(),
        gffread_executable=args.gffread_executable,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
