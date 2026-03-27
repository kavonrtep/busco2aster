"""Build stable QC tables from completed BUSCO runs."""

from __future__ import annotations

import argparse
from pathlib import Path

from .busco import (
    BUSCO_RECORD_FIELDNAMES,
    BUSCO_SUMMARY_FIELDNAMES,
    build_busco_tables,
)
from .manifest import write_tsv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Validated manifest TSV path.")
    parser.add_argument("--summary-output", required=True, help="Output path for busco_summary.tsv.")
    parser.add_argument("--records-output", required=True, help="Output path for parsed BUSCO records.")
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve stable BUSCO output paths.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    summary_rows, record_rows = build_busco_tables(Path(args.manifest), repo_root)
    write_tsv(Path(args.summary_output), summary_rows, BUSCO_SUMMARY_FIELDNAMES)
    write_tsv(Path(args.records_output), record_rows, BUSCO_RECORD_FIELDNAMES)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
