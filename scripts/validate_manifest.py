"""Validate the normalized sample manifest and write metadata outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from .manifest import load_samples_manifest, validate_manifest_rows, write_tsv


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Normalized sample TSV.")
    parser.add_argument("--validated-output", required=True, help="Validated TSV output path.")
    parser.add_argument("--taxon-map-output", required=True, help="Taxon mapping TSV output path.")
    parser.add_argument(
        "--repo-root",
        required=True,
        help="Repository root used to resolve relative assembly paths.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    manifest_path = Path(args.input)
    validated_output = Path(args.validated_output)
    taxon_map_output = Path(args.taxon_map_output)
    repo_root = Path(args.repo_root).resolve()

    rows = load_samples_manifest(manifest_path)
    validated_rows, taxon_rows = validate_manifest_rows(rows, repo_root=repo_root)

    write_tsv(
        validated_output,
        validated_rows,
        ["sample_id", "taxon_id", "sanitized_taxon_id", "assembly_fasta"],
    )
    write_tsv(
        taxon_map_output,
        taxon_rows,
        ["taxon_id", "sanitized_taxon_id"],
    )


if __name__ == "__main__":
    main()

