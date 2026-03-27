"""Aggregate per-locus IQ-TREE outputs into stable workflow artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path

from .gene_trees import aggregate_gene_tree_outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retained", required=True, help="Path to results/qc/retained_loci.tsv.")
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to results/gene_trees/gene_tree_manifest.tsv.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to results/gene_trees/gene_trees.raw.tre.",
    )
    parser.add_argument(
        "--support-mode",
        required=True,
        choices=("abayes", "ufboot"),
        help="Support format expected in the per-locus gene trees.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    aggregate_gene_tree_outputs(
        retained_path=Path(args.retained),
        manifest_path=Path(args.manifest),
        aggregate_path=Path(args.output),
        support_mode=args.support_mode,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
