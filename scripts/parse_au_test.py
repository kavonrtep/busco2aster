"""Parse IQ-TREE AU topology test results from a .iqtree report file."""

from __future__ import annotations

import csv
import re
from pathlib import Path

# Matches one data row in the USER TREES table.
# Format: index logL deltaL value sign [value sign ...]
# The p-AU column is always last.
_USER_TREES_ROW_RE = re.compile(
    r"^\s*(?P<index>\d+)"
    r"\s+(?P<logL>-?[0-9]+(?:\.[0-9]+)?)"
    r"\s+(?P<deltaL>[0-9]+(?:\.[0-9]+)?)"
    # remaining pairs of (value, sign) for bp-RELL, p-KH, p-SH, p-WKH, p-WSH, c-ELW, p-AU
    r"(?P<rest>\s+[0-9.]+\s+[+-].*)"
)

# Names for the sign-marked columns (in order as they appear in the table).
_SIGNED_COLUMNS = ["bp_rell", "p_kh", "p_sh", "p_wkh", "p_wsh", "c_elw", "p_au"]

AU_RESULTS_COLUMNS = [
    "tree_index",
    "description",
    "source",
    "logL",
    "deltaL",
    "bp_rell",
    "p_kh",
    "p_sh",
    "p_wkh",
    "p_wsh",
    "c_elw",
    "p_au",
    "in_confidence_set",
]


def _parse_rest_columns(rest: str) -> dict[str, str]:
    """Parse alternating value/sign tokens from the right portion of a table row."""
    tokens = rest.split()
    result: dict[str, str] = {}
    col_iter = iter(_SIGNED_COLUMNS)
    i = 0
    while i < len(tokens) and True:
        try:
            col = next(col_iter)
        except StopIteration:
            break
        if i >= len(tokens):
            break
        result[col] = tokens[i]
        i += 2  # skip the +/- sign token
    return result


def parse_au_test_iqtree(iqtree_path: Path) -> list[dict]:
    """Return one row per candidate tree from the USER TREES section."""
    lines = iqtree_path.read_text(encoding="utf-8").splitlines()

    # Find the USER TREES header line.
    header_idx = next(
        (i for i, ln in enumerate(lines) if ln.strip() == "USER TREES"),
        None,
    )
    if header_idx is None:
        raise ValueError(f"USER TREES section not found in {iqtree_path}")

    # Find the column separator line (row of dashes) after the header.
    # The column header is between the first separator (after "USER TREES")
    # and the second separator, which precedes the data rows.
    sep_indices = [
        i for i in range(header_idx, len(lines))
        if re.match(r"^[-]+$", lines[i].strip()) and lines[i].strip()
    ]
    if len(sep_indices) < 1:
        raise ValueError(f"Could not find column separator in USER TREES section of {iqtree_path}")

    # Data rows start after the last separator found in the section.
    data_start = sep_indices[-1] + 1

    rows: list[dict] = []
    for line in lines[data_start:]:
        m = _USER_TREES_ROW_RE.match(line)
        if m is None:
            # Empty line or footnote line — stop if we've already collected rows.
            if rows and not line.strip():
                break
            continue
        cols = _parse_rest_columns(m.group("rest"))
        row = {
            "tree_index": int(m.group("index")),
            "logL": float(m.group("logL")),
            "deltaL": float(m.group("deltaL")),
            **{k: float(v) for k, v in cols.items()},
        }
        rows.append(row)

    if not rows:
        raise ValueError(f"No tree rows found in USER TREES section of {iqtree_path}")

    return rows


def join_with_manifest(
    au_rows: list[dict],
    manifest_path: Path,
) -> list[dict]:
    """Attach description and source from the candidate tree manifest."""
    manifest: dict[int, dict] = {}
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            manifest[int(row["tree_index"])] = row

    result: list[dict] = []
    for row in au_rows:
        idx = row["tree_index"]
        meta = manifest.get(idx, {})
        p_au = row.get("p_au")
        in_cs = (p_au is not None and p_au >= 0.05)
        result.append(
            {
                "tree_index": idx,
                "description": meta.get("description", ""),
                "source": meta.get("source", ""),
                "logL": row.get("logL"),
                "deltaL": row.get("deltaL"),
                "bp_rell": row.get("bp_rell"),
                "p_kh": row.get("p_kh"),
                "p_sh": row.get("p_sh"),
                "p_wkh": row.get("p_wkh"),
                "p_wsh": row.get("p_wsh"),
                "c_elw": row.get("c_elw"),
                "p_au": p_au,
                "in_confidence_set": in_cs,
            }
        )

    return result


def write_au_results(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=AU_RESULTS_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
