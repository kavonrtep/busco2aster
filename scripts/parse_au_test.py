"""Parse IQ-TREE AU topology test results from a .iqtree report file."""

from __future__ import annotations

import csv
import re
from pathlib import Path

# Matches one data row in the USER TREES table.
# Format: index logL deltaL [value sign [value sign ...]]
# With >=2 candidate trees IQ-TREE emits bp-RELL, p-KH, p-SH, p-WKH, p-WSH,
# c-ELW, p-AU — each a value plus "+"/"-" sign. With a single candidate tree
# those columns are omitted and only index/logL/deltaL are printed.
_USER_TREES_ROW_RE = re.compile(
    r"^\s*(?P<index>\d+)"
    r"\s+(?P<logL>-?[0-9]+(?:\.[0-9]+)?)"
    r"\s+(?P<deltaL>[0-9]+(?:\.[0-9]+)?)"
    r"(?P<rest>(?:\s+[0-9.]+\s+[+-].*)?)\s*$"
)

# Mapping from IQ-TREE column header tokens to our normalised TSV column names.
# Which signed columns are emitted depends on the test flags: a default
# `--test-au` run emits bp-RELL, p-KH, p-SH, c-ELW, p-AU; adding `--test-weight`
# inserts p-WKH and p-WSH between p-SH and c-ELW.
_HEADER_TOKEN_TO_COLUMN = {
    "bp-RELL": "bp_rell",
    "p-KH": "p_kh",
    "p-SH": "p_sh",
    "p-WKH": "p_wkh",
    "p-WSH": "p_wsh",
    "c-ELW": "c_elw",
    "p-AU": "p_au",
}
_ALL_SIGNED_COLUMNS = list(_HEADER_TOKEN_TO_COLUMN.values())

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


def _parse_signed_columns_header(header_line: str) -> list[str]:
    """Map column header tokens after deltaL onto our TSV column names."""
    tokens = header_line.split()
    try:
        delta_idx = tokens.index("deltaL")
    except ValueError as exc:
        raise ValueError(f"deltaL not in column header: {header_line!r}") from exc
    return [
        _HEADER_TOKEN_TO_COLUMN[t]
        for t in tokens[delta_idx + 1:]
        if t in _HEADER_TOKEN_TO_COLUMN
    ]


def _parse_rest_columns(rest: str, signed_columns: list[str]) -> dict[str, str | None]:
    """Parse alternating value/sign tokens from the right portion of a table row.

    Columns absent from the table (e.g. p-WKH/p-WSH when --test-weight is off)
    are returned as None.
    """
    tokens = rest.split()
    result: dict[str, str | None] = {col: None for col in _ALL_SIGNED_COLUMNS}
    for i, col in enumerate(signed_columns):
        token_idx = 2 * i  # value, sign, value, sign, ...
        if token_idx >= len(tokens):
            break
        result[col] = tokens[token_idx]
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

    # Locate the column-header separator (a row of dashes immediately under
    # the "Tree  logL  deltaL ..." column header). Limit the search to the
    # USER TREES section itself — later sections (footnotes, model report)
    # also use dash separators and would otherwise be matched.
    column_header_idx = next(
        (
            i for i in range(header_idx + 1, len(lines))
            if lines[i].lstrip().startswith("Tree") and "logL" in lines[i]
        ),
        None,
    )
    if column_header_idx is None:
        raise ValueError(
            f"Could not find 'Tree  logL  deltaL' column header in USER TREES "
            f"section of {iqtree_path}"
        )
    if column_header_idx + 1 >= len(lines) or not re.fullmatch(
        r"-+", lines[column_header_idx + 1].strip()
    ):
        raise ValueError(
            f"Expected dash separator under column header in {iqtree_path}"
        )
    signed_columns = _parse_signed_columns_header(lines[column_header_idx])
    data_start = column_header_idx + 2

    rows: list[dict] = []
    for line in lines[data_start:]:
        m = _USER_TREES_ROW_RE.match(line)
        if m is None:
            # Empty line or footnote line — stop if we've already collected rows.
            if rows and not line.strip():
                break
            continue
        cols = _parse_rest_columns(m.group("rest"), signed_columns)
        row = {
            "tree_index": int(m.group("index")),
            "logL": float(m.group("logL")),
            "deltaL": float(m.group("deltaL")),
            **{k: (float(v) if v is not None else None) for k, v in cols.items()},
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
        # Single-candidate runs omit p-AU entirely; the lone tree is
        # vacuously in the confidence set.
        in_cs = True if p_au is None else p_au >= 0.05
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
