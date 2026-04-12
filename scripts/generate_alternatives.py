"""Generate alternative species tree topologies via NNI swaps at contested branches."""

from __future__ import annotations

import csv
import itertools
from pathlib import Path

from scripts.tree_utils import (
    TreeNode,
    branch_label_map,
    informative_branch_key,
    iter_nodes,
    leaf_labels,
    parse_newick,
    read_newick,
    render_newick,
    write_newick,
)


def _find_node_by_branch_id(root: TreeNode, branch_id: str) -> TreeNode | None:
    """Return the internal node whose stable branch label matches branch_id."""
    total_taxa = leaf_labels(root)
    label_map = branch_label_map(root)
    for node in iter_nodes(root):
        branch_key = informative_branch_id_key(node, total_taxa)
        if branch_key is not None and label_map.get(branch_key) == branch_id:
            return node
    return None


def informative_branch_id_key(node: TreeNode, total_taxa: set[str]) -> str | None:
    return informative_branch_key(node, total_taxa)


def apply_nni_swap(root: TreeNode, branch_id: str, child_index: int = 1) -> TreeNode:
    """Return a new tree with one NNI swap applied at the branch identified by branch_id.

    The swap exchanges node.children[child_index] with the sibling of node (the other
    child of node.parent). child_index=0 and child_index=1 produce the two distinct
    NNI alternatives for that branch.

    Raises ValueError if the branch is not found or NNI is not applicable (e.g., root).
    """
    cloned = root.clone()
    target = _find_node_by_branch_id(cloned, branch_id)
    if target is None:
        raise ValueError(f"Branch {branch_id!r} not found in tree.")
    parent = target.parent
    if parent is None:
        raise ValueError(f"Branch {branch_id!r} is the root; NNI not applicable.")
    if len(target.children) < 2:
        raise ValueError(
            f"Node for branch {branch_id!r} has fewer than 2 children; NNI not applicable."
        )
    if len(parent.children) < 2:
        raise ValueError(f"Parent of branch {branch_id!r} has fewer than 2 children.")

    # Identify the sibling: the other child of parent.
    siblings = [c for c in parent.children if c is not target]
    if not siblings:
        raise ValueError(f"No sibling found for branch {branch_id!r}.")
    sibling = siblings[0]

    # The child of target being swapped out.
    swap_child = target.children[child_index % len(target.children)]

    # Perform the swap in-place on the clone.
    # target.children[child_index] ← sibling
    # parent replaces sibling with swap_child
    target.children[child_index % len(target.children)] = sibling
    sibling.parent = target

    parent_idx = parent.children.index(sibling)
    parent.children[parent_idx] = swap_child
    swap_child.parent = parent

    return cloned


def _load_hypothesis_trees(path: Path) -> list[TreeNode]:
    """Read a multi-tree Newick file (one tree per non-empty line)."""
    trees = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            trees.append(parse_newick(line))
    return trees


def generate_candidate_trees(
    species_tree_path: Path,
    contested_rows: list[dict],
    max_contested_branches: int = 4,
    hypothesis_path: Path | None = None,
) -> tuple[list[TreeNode], list[str], list[str]]:
    """Build the candidate tree set for AU testing.

    Returns (trees, descriptions, sources) where:
      - trees[0] is always the original wASTRAL tree
      - subsequent entries are NNI alternatives (combinations of swaps)
      - user-supplied hypothesis trees follow, if any
    """
    root = read_newick(species_tree_path)

    # Rank contested branches by ascending pp1 (most contested first).
    ranked = sorted(contested_rows, key=lambda r: float(r["pp1"]))
    selected = ranked[:max_contested_branches]

    trees: list[TreeNode] = [root.clone()]
    descriptions: list[str] = ["wASTRAL best tree"]
    sources: list[str] = ["original"]

    # Generate all 2^k - 1 non-empty subsets of selected branches.
    # For each subset, apply NNI swaps at each branch in the subset.
    # We use child_index=1 as the canonical "best alternative" swap.
    branch_ids = [r["branch_id"] for r in selected]
    for size in range(1, len(branch_ids) + 1):
        for combo in itertools.combinations(range(len(branch_ids)), size):
            swapped_root = root.clone()
            desc_parts = []
            for idx in combo:
                bid = branch_ids[idx]
                try:
                    swapped_root = apply_nni_swap(swapped_root, bid, child_index=1)
                    desc_parts.append(f"{bid} swapped")
                except ValueError:
                    # Skip this combination if the swap fails (e.g., branch not found after
                    # previous swap altered structure — rare in practice for independent branches).
                    break
            else:
                trees.append(swapped_root)
                descriptions.append("; ".join(desc_parts))
                sources.append("generated")

    # Append user-supplied hypothesis trees.
    if hypothesis_path is not None and hypothesis_path.is_file():
        for hyp_tree in _load_hypothesis_trees(hypothesis_path):
            trees.append(hyp_tree)
            descriptions.append(f"user hypothesis ({hypothesis_path.name})")
            sources.append("user_hypothesis")

    return trees, descriptions, sources


def write_candidate_trees(trees: list[TreeNode], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [render_newick(t) + ";" for t in trees]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_candidate_manifest(
    descriptions: list[str],
    sources: list[str],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["tree_index", "description", "source"],
            delimiter="\t",
        )
        writer.writeheader()
        for i, (desc, src) in enumerate(zip(descriptions, sources)):
            writer.writerow({"tree_index": i + 1, "description": desc, "source": src})
