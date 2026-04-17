"""Minimal Newick parsing helpers for reporting and concordance joins."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class TreeNode:
    label: str | None = None
    length: str | None = None
    children: list["TreeNode"] = field(default_factory=list)
    parent: "TreeNode | None" = field(default=None, repr=False)

    def is_leaf(self) -> bool:
        return not self.children

    def clone(self) -> "TreeNode":
        cloned_children = [child.clone() for child in self.children]
        node = TreeNode(label=self.label, length=self.length, children=cloned_children)
        for child in cloned_children:
            child.parent = node
        return node


def _parse_label(text: str, start: int) -> tuple[str | None, int]:
    if start >= len(text):
        return None, start
    if text[start] == "'":
        cursor = start + 1
        parts: list[str] = []
        while cursor < len(text):
            char = text[cursor]
            if char == "'":
                if cursor + 1 < len(text) and text[cursor + 1] == "'":
                    parts.append("'")
                    cursor += 2
                    continue
                return "".join(parts), cursor + 1
            parts.append(char)
            cursor += 1
        raise ValueError("Unterminated quoted label in Newick string.")

    cursor = start
    while cursor < len(text) and text[cursor] not in ",():;":
        cursor += 1
    if cursor == start:
        return None, start
    return text[start:cursor], cursor


def _parse_length(text: str, start: int) -> tuple[str | None, int]:
    if start >= len(text) or text[start] != ":":
        return None, start
    cursor = start + 1
    while cursor < len(text) and text[cursor] not in ",);":
        cursor += 1
    return text[start + 1:cursor], cursor


def parse_newick(text: str) -> TreeNode:
    stripped = text.strip()
    if not stripped.endswith(";"):
        raise ValueError("Newick text must end with ';'.")

    def parse_subtree(start: int) -> tuple[TreeNode, int]:
        if stripped[start] == "(":
            cursor = start + 1
            children: list[TreeNode] = []
            while True:
                child, cursor = parse_subtree(cursor)
                children.append(child)
                if stripped[cursor] == ",":
                    cursor += 1
                    continue
                if stripped[cursor] == ")":
                    cursor += 1
                    break
                raise ValueError(f"Unexpected token {stripped[cursor]!r} in Newick string.")
            label, cursor = _parse_label(stripped, cursor)
            length, cursor = _parse_length(stripped, cursor)
            node = TreeNode(label=label, length=length, children=children)
            for child in children:
                child.parent = node
            return node, cursor

        label, cursor = _parse_label(stripped, start)
        if label is None:
            raise ValueError(f"Expected label at position {start} in Newick string.")
        length, cursor = _parse_length(stripped, cursor)
        return TreeNode(label=label, length=length), cursor

    root, cursor = parse_subtree(0)
    if stripped[cursor] != ";":
        raise ValueError("Unexpected trailing content in Newick string.")
    return root


def read_newick(path: Path) -> TreeNode:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError(f"Expected exactly one Newick tree in {path}, found {len(lines)}.")
    return parse_newick(lines[0])


def iter_nodes(root: TreeNode) -> list[TreeNode]:
    nodes: list[TreeNode] = []

    def visit(node: TreeNode) -> None:
        nodes.append(node)
        for child in node.children:
            visit(child)

    visit(root)
    return nodes


def leaf_labels(node: TreeNode) -> set[str]:
    if node.is_leaf():
        if node.label is None:
            raise ValueError("Leaf node is missing a label.")
        return {node.label}
    labels: set[str] = set()
    for child in node.children:
        labels.update(leaf_labels(child))
    return labels


def _quote_label(label: str) -> str:
    if not label:
        return ""
    safe_chars = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-")
    if all(char in safe_chars for char in label):
        return label
    return "'" + label.replace("'", "''") + "'"


def render_newick(node: TreeNode) -> str:
    if node.is_leaf():
        text = _quote_label(node.label or "")
    else:
        rendered_children = ",".join(render_newick(child) for child in node.children)
        text = f"({rendered_children})"
        if node.label:
            text += _quote_label(node.label)
    if node.length is not None:
        text += f":{node.length}"
    return text


def write_newick(path: Path, root: TreeNode) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_newick(root) + ";\n", encoding="utf-8")


def canonical_branch_key(split_a: set[str], split_b: set[str]) -> str:
    left = tuple(sorted(split_a))
    right = tuple(sorted(split_b))
    if (len(left), left) > (len(right), right):
        left, right = right, left
    return f"{','.join(left)}||{','.join(right)}"


def informative_branch_key(node: TreeNode, total_taxa: set[str]) -> str | None:
    if node.parent is None or node.is_leaf():
        return None
    node_taxa = leaf_labels(node)
    other_taxa = total_taxa - node_taxa
    if len(node_taxa) < 2 or len(other_taxa) < 2:
        return None
    return canonical_branch_key(node_taxa, other_taxa)


def internal_branch_key_map(root: TreeNode) -> dict[str, TreeNode]:
    total_taxa = leaf_labels(root)
    branch_map: dict[str, TreeNode] = {}
    for node in iter_nodes(root):
        branch_key = informative_branch_key(node, total_taxa)
        if branch_key is None:
            continue
        branch_map[branch_key] = node
    return branch_map


def branch_label_map(root: TreeNode) -> dict[str, str]:
    branch_keys = sorted(internal_branch_key_map(root))
    return {branch_key: f"B{index}" for index, branch_key in enumerate(branch_keys, start=1)}


def relabel_tree_with_branch_ids(root: TreeNode) -> TreeNode:
    cloned_root = root.clone()
    label_map = branch_label_map(cloned_root)
    total_taxa = leaf_labels(cloned_root)
    for node in iter_nodes(cloned_root):
        branch_key = informative_branch_key(node, total_taxa)
        if branch_key is None:
            if not node.is_leaf():
                node.label = None
            continue
        node.label = label_map[branch_key]
    return cloned_root


def fill_missing_branch_lengths(root: TreeNode, default: str = "0") -> None:
    for node in iter_nodes(root):
        if node is root:
            continue
        if node.length is None:
            node.length = default


def canonical_topology_key(root: TreeNode) -> str:
    return ";".join(sorted(internal_branch_key_map(root)))
