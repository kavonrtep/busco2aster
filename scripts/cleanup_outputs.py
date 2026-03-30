"""Apply retention-tier cleanup policies to workflow outputs."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


RETENTION_MODES = ("debug", "resume_alignments", "resume_gene_trees", "final_report")

RETENTION_PATTERNS: dict[str, tuple[str, ...]] = {
    "debug": (),
    "resume_alignments": (
        "work/busco/*/raw",
        "results/busco/*/raw",
    ),
    "resume_gene_trees": (
        "work/busco/*/raw",
        "results/busco/*/raw",
        "results/loci/raw_fastas",
        "results/loci/raw_fastas.complete",
        "results/loci/logs",
        "results/gene_trees/per_locus",
    ),
    "final_report": (
        "work/busco/*/raw",
        "results/busco/*/raw",
        "results/loci/raw_fastas",
        "results/loci/raw_fastas.complete",
        "results/loci/logs",
        "results/loci/alignments",
        "results/loci/alignments.complete",
        "results/gene_trees/per_locus",
        "results/gene_trees/gene_trees.command.sh",
        "results/gene_trees/gene_trees.iqtree",
        "results/gene_trees/gene_trees.log",
        "results/gene_trees/gene_trees.treefile",
        "results/gene_trees/gene_trees.best_model.nex",
        "results/gene_trees/gene_trees.best_scheme",
        "results/gene_trees/gene_trees.best_scheme.nex",
        "results/gene_trees/gene_trees.model.gz",
        "results/gene_trees/gene_trees.ckp.gz",
        "results/gene_trees/gene_trees.parstree",
        "results/concordance/*.command.sh",
        "results/concordance/*.log",
        "results/concordance/*.iqtree",
        "results/concordance/*.best_model.nex",
        "results/concordance/*.treefile",
        "results/concordance/*.ckp.gz",
        "results/concordance/*.annotated.tre",
    ),
}


def list_cleanup_targets(repo_root: Path, mode: str) -> list[Path]:
    if mode not in RETENTION_PATTERNS:
        raise ValueError(f"Unsupported cleanup mode {mode!r}. Expected one of {RETENTION_MODES}.")

    targets: set[Path] = set()
    for pattern in RETENTION_PATTERNS[mode]:
        targets.update(path for path in repo_root.glob(pattern) if path.exists())
    return sorted(targets, key=lambda path: (len(path.parts), path.as_posix()))


def remove_cleanup_targets(targets: list[Path]) -> list[Path]:
    removed: list[Path] = []
    for path in sorted(targets, key=lambda candidate: len(candidate.parts), reverse=True):
        if not path.exists():
            continue
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink()
        removed.append(path)
    return removed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--mode", required=True, choices=RETENTION_MODES)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    targets = list_cleanup_targets(repo_root, args.mode)

    for path in targets:
        print(path.relative_to(repo_root).as_posix())

    if not args.dry_run:
        remove_cleanup_targets(targets)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
