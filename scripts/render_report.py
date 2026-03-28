"""Render the final workflow report."""

from __future__ import annotations

import argparse
from pathlib import Path

from .report import render_report, write_report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--busco-summary", required=True, help="Path to results/qc/busco_summary.tsv.")
    parser.add_argument("--retained-loci", required=True, help="Path to results/qc/retained_loci.tsv.")
    parser.add_argument(
        "--gene-tree-manifest",
        required=True,
        help="Path to results/gene_trees/gene_tree_manifest.tsv.",
    )
    parser.add_argument("--species-tree", required=True, help="Path to the default species-tree Newick file.")
    parser.add_argument("--species-tree-log", required=True, help="Path to the species-tree log file.")
    parser.add_argument("--backend", required=True, help="Species-tree backend label.")
    parser.add_argument("--output", required=True, help="Output Markdown path.")
    parser.add_argument(
        "--concordance-path",
        action="append",
        default=[],
        help="Optional concordance artifact path. Missing files are ignored.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    content = render_report(
        busco_summary_path=Path(args.busco_summary),
        retained_loci_path=Path(args.retained_loci),
        gene_tree_manifest_path=Path(args.gene_tree_manifest),
        species_tree_path=Path(args.species_tree),
        species_tree_log_path=Path(args.species_tree_log),
        species_tree_backend=args.backend,
        concordance_paths=[Path(path) for path in args.concordance_path],
    )
    write_report(Path(args.output), content)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
