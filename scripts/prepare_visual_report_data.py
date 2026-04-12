"""CLI wrapper for preparing the visual report data bundle."""

from __future__ import annotations

import argparse
from pathlib import Path

from .report_data import build_report_data_bundle


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root used to resolve relative paths.")
    parser.add_argument("--busco-summary", required=True)
    parser.add_argument("--locus-matrix", required=True)
    parser.add_argument("--retained-loci", required=True)
    parser.add_argument("--gene-tree-manifest", required=True)
    parser.add_argument("--species-tree", required=True)
    parser.add_argument("--gcf-stat", required=True)
    parser.add_argument("--gcf-branch", required=True)
    parser.add_argument("--scfl-stat", required=True)
    parser.add_argument("--scfl-branch", required=True)
    parser.add_argument("--quartet-freqquad", required=True)
    parser.add_argument("--alignment-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sequence-type", default="protein")
    parser.add_argument(
        "--topology-tests-results-dir",
        default=None,
        help="Directory containing topology test result TSVs to copy into the report data dir.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir)
    build_report_data_bundle(
        repo_root=Path(args.repo_root).resolve(),
        busco_summary_path=Path(args.busco_summary),
        locus_matrix_path=Path(args.locus_matrix),
        retained_loci_path=Path(args.retained_loci),
        gene_tree_manifest_path=Path(args.gene_tree_manifest),
        species_tree_path=Path(args.species_tree),
        gcf_stat_path=Path(args.gcf_stat),
        gcf_branch_path=Path(args.gcf_branch),
        scfl_stat_path=Path(args.scfl_stat),
        scfl_branch_path=Path(args.scfl_branch),
        quartet_freqquad_path=Path(args.quartet_freqquad),
        alignment_dir=Path(args.alignment_dir),
        output_dir=output_dir,
        sequence_type=args.sequence_type,
    )
    if args.topology_tests_results_dir is not None:
        import shutil

        src_dir = Path(args.topology_tests_results_dir)
        _TOPOLOGY_REPORT_FILES = [
            "branch_quartet_support.tsv",
            "contested_branches.tsv",
            "candidate_trees_manifest.tsv",
            "au_test_results.tsv",
        ]
        for fname in _TOPOLOGY_REPORT_FILES:
            src = src_dir / fname
            if src.is_file():
                shutil.copy2(src, output_dir / fname)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
