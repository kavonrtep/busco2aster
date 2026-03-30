"""Finalize IQ-TREE directory-mode outputs into stable workflow artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path

from .gene_trees import aggregate_directory_mode_gene_tree_outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--best-model-nex", required=True)
    parser.add_argument("--treefile", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--support-mode",
        required=True,
        choices=("abayes", "ufboot"),
        help="Support format expected in the directory-mode gene trees.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    aggregate_directory_mode_gene_tree_outputs(
        best_model_nex_path=Path(args.best_model_nex),
        treefile_path=Path(args.treefile),
        report_path=Path(args.report),
        manifest_path=Path(args.manifest),
        aggregate_path=Path(args.output),
        support_mode=args.support_mode,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
