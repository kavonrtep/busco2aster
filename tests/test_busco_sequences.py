"""Tests for the concatenated BUSCO sequence layout and its per-locus predecessor."""

import csv
import json
import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.busco import (
    build_busco_tables,
    busco_output_paths,
    load_busco_artifact_paths,
    standardize_busco_run,
)
from scripts.busco_sequences import (
    GFF_MARKER_PREFIX,
    LEGACY_DIR_BY_CATEGORY,
    parse_fasta_records_text,
    read_fasta_records_for_id,
    read_gff_text,
)
from scripts.manifest import write_tsv


REPO_ROOT = Path(__file__).resolve().parents[1]

FULL_TABLE_HEADER = (
    "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\t"
    "OrthoDB url\tDescription\n"
)


def locus_ids(count: int) -> list[str]:
    return [f"{index}at9999" for index in range(1, count + 1)]


def protein_text(busco_id: str, copy: str = "") -> str:
    return f">chr1_{busco_id}{copy}:1-90\nMAKKQ\nWWYYE\n"


def gff_text(busco_id: str) -> str:
    return (
        "##gff-version 3\n"
        f"chr1\tminiprot\tmRNA\t1\t90\t.\t+\t.\tID=MP_{busco_id};Target={busco_id} 1 30\n"
        f"chr1\tminiprot\tCDS\t1\t90\t.\t+\t0\tParent=MP_{busco_id}\n"
    )


def write_raw_busco_run(raw_root: Path, single_copy: list[str], duplicated: list[str]) -> None:
    """Write what BUSCO leaves in its scratch directory: two files per locus."""
    run_root = raw_root / "run_fabales_odb12"
    sequence_root = run_root / "busco_sequences"
    single_copy_dir = sequence_root / "single_copy_busco_sequences"
    multi_copy_dir = sequence_root / "multi_copy_busco_sequences"
    for directory in (single_copy_dir, multi_copy_dir):
        directory.mkdir(parents=True, exist_ok=True)

    for busco_id in single_copy:
        (single_copy_dir / f"{busco_id}.faa").write_text(protein_text(busco_id), encoding="utf-8")
        (single_copy_dir / f"{busco_id}.gff").write_text(gff_text(busco_id), encoding="utf-8")

    for busco_id in duplicated:
        # A duplicated locus carries more than one protein record in one file.
        (multi_copy_dir / f"{busco_id}.faa").write_text(
            protein_text(busco_id, "a") + protein_text(busco_id, "b"), encoding="utf-8"
        )
        (multi_copy_dir / f"{busco_id}.gff").write_text(gff_text(busco_id), encoding="utf-8")

    run_root.joinpath("short_summary.json").write_text(
        json.dumps(
            {
                "parameters": {"in": "assembly.fa.gz", "mode": "genome"},
                "lineage_dataset": {
                    "name": "fabales_odb12",
                    "creation_date": "2026-01-01",
                    "number_of_buscos": len(single_copy) + len(duplicated),
                    "number_of_species": 10,
                },
                "versions": {"busco": "6.0.0"},
                "results": {
                    "one_line_summary": "C:100%",
                    "n_markers": len(single_copy) + len(duplicated),
                    "Complete BUSCOs": len(single_copy) + len(duplicated),
                    "Single copy BUSCOs": len(single_copy),
                    "Multi copy BUSCOs": len(duplicated),
                    "Fragmented BUSCOs": 0,
                    "Missing BUSCOs": 0,
                    "Complete percentage": 100.0,
                    "Single copy percentage": 100.0,
                    "Multi copy percentage": 0.0,
                    "Fragmented percentage": 0.0,
                    "Missing percentage": 0.0,
                },
            }
        ),
        encoding="utf-8",
    )

    rows = [
        f"{busco_id}\tComplete\tchr1\t1\t90\t+\t200.0\t30\thttps://example.org\tprotein"
        for busco_id in single_copy
    ]
    for busco_id in duplicated:
        rows.append(
            f"{busco_id}\tDuplicated\tchr1\t1\t90\t+\t200.0\t30\thttps://example.org\tprotein"
        )
        rows.append(
            f"{busco_id}\tDuplicated\tchr2\t1\t90\t-\t199.0\t30\thttps://example.org\tprotein"
        )
    run_root.joinpath("full_table.tsv").write_text(
        FULL_TABLE_HEADER + "\n".join(rows) + "\n", encoding="utf-8"
    )


def standardize_in(repo_root: Path, sample_id: str, raw_root: Path) -> dict[str, str]:
    """Run the copy-out step with `repo_root` as the working directory."""
    cwd = Path.cwd()
    try:
        os.chdir(repo_root)
        standardize_busco_run(
            sample_id=sample_id,
            sample_dir=Path(busco_output_paths(sample_id)["sample_root"]),
            raw_root=raw_root,
        )
    finally:
        os.chdir(cwd)
    return load_busco_artifact_paths(repo_root / busco_output_paths(sample_id)["paths"])


class BuscoSequenceLayoutTests(unittest.TestCase):
    def test_sequence_file_count_does_not_grow_with_the_number_of_loci(self):
        """The inode cost of a sample is fixed; BUSCO's own output is two files per locus."""
        counts = {}
        for locus_count in (3, 12):
            with TemporaryDirectory() as tmpdir:
                repo_root = Path(tmpdir)
                raw_root = repo_root / "work" / "busco" / "sample1" / "raw"
                write_raw_busco_run(raw_root, locus_ids(locus_count), [])
                standardize_in(repo_root, "sample1", raw_root)
                sequence_root = repo_root / busco_output_paths("sample1")["sequence_root"]
                counts[locus_count] = sorted(path.name for path in sequence_root.rglob("*"))

        self.assertEqual(counts[3], counts[12])
        # Only the categories BUSCO reported; two files each, whatever the locus count.
        self.assertEqual(
            counts[3],
            ["multi_copy.faa", "multi_copy.gff", "single_copy.faa", "single_copy.gff"],
        )

    def test_every_locus_round_trips_through_the_concatenated_files(self):
        with TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            raw_root = repo_root / "work" / "busco" / "sample1" / "raw"
            single_copy = locus_ids(4)
            duplicated = ["77at9999"]
            write_raw_busco_run(raw_root, single_copy, duplicated)
            artifact_paths = standardize_in(repo_root, "sample1", raw_root)

            for busco_id in single_copy:
                records = read_fasta_records_for_id(
                    artifact_paths["single_copy_busco_sequences"], busco_id, repo_root
                )
                self.assertEqual(
                    [sequence for _, sequence in records],
                    [sequence for _, sequence in parse_fasta_records_text(protein_text(busco_id))],
                )
                self.assertEqual(
                    read_gff_text(artifact_paths["single_copy_busco_gff"], busco_id, repo_root),
                    gff_text(busco_id),
                )

            duplicated_records = read_fasta_records_for_id(
                artifact_paths["multi_copy_busco_sequences"], "77at9999", repo_root
            )
            self.assertEqual(len(duplicated_records), 2)
            self.assertTrue(
                all(header.startswith("77at9999 ") for header, _ in duplicated_records)
            )

    def test_absent_locus_is_reported_rather_than_silently_matched(self):
        with TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            raw_root = repo_root / "work" / "busco" / "sample1" / "raw"
            write_raw_busco_run(raw_root, locus_ids(3), [])
            artifact_paths = standardize_in(repo_root, "sample1", raw_root)

            self.assertEqual(
                read_fasta_records_for_id(
                    artifact_paths["single_copy_busco_sequences"], "9999at9999", repo_root
                ),
                [],
            )

    def test_tables_from_both_layouts_describe_the_same_sequences(self):
        """A per-locus tree and a concatenated one produce identical QC records."""
        with TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            raw_root = repo_root / "work" / "busco" / "sample1" / "raw"
            single_copy = locus_ids(3)
            write_raw_busco_run(raw_root, single_copy, ["77at9999"])
            standardize_in(repo_root, "sample1", raw_root)

            manifest_path = repo_root / "samples.validated.tsv"
            write_tsv(
                manifest_path,
                [
                    {
                        "sample_id": "sample1",
                        "taxon_id": "Sample one",
                        "sanitized_taxon_id": "sample_one",
                        "assembly_fasta": "assemblies/sample1.fa.gz",
                    }
                ],
                ["sample_id", "taxon_id", "sanitized_taxon_id", "assembly_fasta"],
            )
            concatenated_summary, concatenated_records = build_busco_tables(
                manifest_path, repo_root
            )

            _rewrite_as_per_locus_layout(repo_root, "sample1")
            legacy_summary, legacy_records = build_busco_tables(manifest_path, repo_root)

        self.assertEqual(
            [(row["busco_id"], row["status"], row["sequence_category"]) for row in legacy_records],
            [
                (row["busco_id"], row["status"], row["sequence_category"])
                for row in concatenated_records
            ],
        )
        self.assertEqual(concatenated_summary[0]["sequence_layout"], "concatenated")
        self.assertEqual(legacy_summary[0]["sequence_layout"], "per_locus")
        self.assertTrue(legacy_records[0]["faa_path"].endswith(f"{legacy_records[0]['busco_id']}.faa"))


def _rewrite_as_per_locus_layout(repo_root: Path, sample_id: str) -> None:
    """Turn a concatenated sample back into what an older release wrote."""
    stable_paths = busco_output_paths(sample_id)
    path_rows = [
        dict(row)
        for row in _read_tsv(repo_root / stable_paths["paths"])
        if not row["artifact"].endswith("_busco_gff")
    ]

    for category, legacy_dir_name in LEGACY_DIR_BY_CATEGORY.items():
        faa_path = repo_root / stable_paths[f"{category}_faa"]
        gff_path = repo_root / stable_paths[f"{category}_gff"]
        if not faa_path.is_file():
            continue
        legacy_dir = repo_root / stable_paths[f"{category}_per_locus_dir"]
        legacy_dir.mkdir(parents=True, exist_ok=True)

        for header, sequence in parse_fasta_records_text(faa_path.read_text(encoding="utf-8")):
            busco_id, _, original = header.partition(" ")
            target = legacy_dir / f"{busco_id}.faa"
            with target.open("a", encoding="utf-8") as handle:
                handle.write(f">{original or busco_id}\n{sequence}\n")

        block_id = ""
        for raw_line in gff_path.read_text(encoding="utf-8").splitlines(keepends=True):
            if raw_line.startswith(GFF_MARKER_PREFIX):
                block_id = raw_line[len(GFF_MARKER_PREFIX) :].strip()
                continue
            if block_id:
                with (legacy_dir / f"{block_id}.gff").open("a", encoding="utf-8") as handle:
                    handle.write(raw_line)

        faa_path.unlink()
        gff_path.unlink()
        for row in path_rows:
            if row["artifact"] == f"{category}_busco_sequences":
                row["path"] = stable_paths[f"{category}_per_locus_dir"]

    for row in path_rows:
        if row["artifact"] == "sequence_layout":
            row["path"] = "per_locus"
    write_tsv(repo_root / stable_paths["paths"], path_rows, ["artifact", "path"])


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


if __name__ == "__main__":
    unittest.main()
