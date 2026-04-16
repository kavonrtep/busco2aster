"""Parse wASTRAL -u 2 annotated tree and identify contested branches."""

from __future__ import annotations

import csv
import re
from pathlib import Path

from scripts.tree_utils import (
    TreeNode,
    branch_label_map,
    informative_branch_key,
    iter_nodes,
    leaf_labels,
    read_newick,
)

# Matches the full wASTRAL -u 2 bracketed annotation block. wASTRAL emits
# keys like CULength, f1..f3, localPP, pp1..pp3, q1..q3 in a non-deterministic
# order across releases (observed: pp* before q* in wASTRAL 1.19+), so we
# capture the bracketed span and then extract each key independently.
_ANNOTATION_BLOCK_RE = re.compile(r"\[[^\[\]]*\]", re.DOTALL)
_NUMBER_PATTERN = r"[0-9.eE+\-]+"
_KEY_PATTERNS = {
    key: re.compile(rf"(?:^|[;,\s\[]){re.escape(key)}=({_NUMBER_PATTERN})")
    for key in ("q1", "q2", "q3", "pp1", "pp2", "pp3")
}

BRANCH_QUARTET_SUPPORT_COLUMNS = [
    "branch_id",
    "taxa_clade",
    "taxa_complement",
    "q1",
    "q2",
    "q3",
    "pp1",
    "pp2",
    "pp3",
    "contested",
]

CONTESTED_BRANCHES_COLUMNS = BRANCH_QUARTET_SUPPORT_COLUMNS


def _parse_annotation(label: str | None) -> dict[str, float] | None:
    if label is None:
        return None
    block_match = _ANNOTATION_BLOCK_RE.search(label)
    if block_match is None:
        return None
    block = block_match.group(0)
    parsed: dict[str, float] = {}
    for key, pattern in _KEY_PATTERNS.items():
        m = pattern.search(block)
        if m is None:
            return None
        parsed[key] = float(m.group(1))
    return parsed


def parse_wastral_u2_tree(
    path: Path,
    contested_threshold: float = 0.95,
) -> list[dict]:
    """Parse wASTRAL -u 2 annotated Newick and return per-branch quartet support rows.

    Each internal informative branch appears exactly once in the output.  When
    both nodes sharing a bipartition carry an annotation (which can happen in
    rooted representations with a bifurcating root), the first annotated node
    encountered in a pre-order traversal is used.
    """
    root = read_newick(path)
    total_taxa = leaf_labels(root)
    label_map = branch_label_map(root)
    seen_branch_ids: set[str] = set()
    rows: list[dict] = []

    for node in iter_nodes(root):
        branch_key = informative_branch_key(node, total_taxa)
        if branch_key is None:
            continue
        branch_id = label_map[branch_key]
        if branch_id in seen_branch_ids:
            continue
        clade_taxa, complement_taxa = branch_key.split("||")

        annotation = _parse_annotation(node.label)
        if annotation is None:
            # Branch exists but has no wASTRAL annotation (e.g., poorly annotated node).
            continue

        seen_branch_ids.add(branch_id)
        contested = annotation["pp1"] < contested_threshold
        rows.append(
            {
                "branch_id": branch_id,
                "taxa_clade": clade_taxa,
                "taxa_complement": complement_taxa,
                "q1": annotation["q1"],
                "q2": annotation["q2"],
                "q3": annotation["q3"],
                "pp1": annotation["pp1"],
                "pp2": annotation["pp2"],
                "pp3": annotation["pp3"],
                "contested": contested,
            }
        )

    rows.sort(key=lambda r: r["branch_id"])
    return rows


def filter_contested_branches(rows: list[dict]) -> list[dict]:
    return [r for r in rows if r["contested"]]


def _write_tsv(rows: list[dict], columns: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_branch_quartet_support(rows: list[dict], path: Path) -> None:
    _write_tsv(rows, BRANCH_QUARTET_SUPPORT_COLUMNS, path)


def write_contested_branches(rows: list[dict], path: Path) -> None:
    _write_tsv(rows, CONTESTED_BRANCHES_COLUMNS, path)
