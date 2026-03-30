import csv
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

from scripts.busco import (
    BuscoValidationError,
    busco_output_paths,
    parse_busco_full_table,
    parse_busco_short_summary,
    summarize_busco_full_table,
)
from scripts.manifest import write_tsv


REPO_ROOT = Path(__file__).resolve().parents[1]


def make_short_summary_payload(
    sample_id: str,
    input_fasta: str,
    single_copy_buscos: int = 1,
    multi_copy_buscos: int = 1,
    fragmented_buscos: int = 1,
    missing_buscos: int = 1,
) -> dict[str, object]:
    complete_buscos = single_copy_buscos + multi_copy_buscos
    n_markers = complete_buscos + fragmented_buscos + missing_buscos
    return {
        "parameters": {
            "in": input_fasta,
            "mode": "euk_genome_min",
            "gene_predictor": "miniprot",
            "out": sample_id,
        },
        "lineage_dataset": {
            "name": "solanales_odb12",
            "creation_date": "2025-07-01",
            "number_of_buscos": str(n_markers),
            "number_of_species": "9",
        },
        "versions": {
            "busco": "6.0.0",
        },
        "results": {
            "one_line_summary": "C:50.0%[S:25.0%,D:25.0%],F:25.0%,M:25.0%,n:4,E:0.0%",
            "Complete percentage": 50.0,
            "Complete BUSCOs": complete_buscos,
            "Single copy percentage": 25.0,
            "Single copy BUSCOs": single_copy_buscos,
            "Multi copy percentage": 25.0,
            "Multi copy BUSCOs": multi_copy_buscos,
            "Fragmented percentage": 25.0,
            "Fragmented BUSCOs": fragmented_buscos,
            "Missing percentage": 25.0,
            "Missing BUSCOs": missing_buscos,
            "n_markers": n_markers,
            "avg_identity": 0.98,
            "internal_stop_codon_count": 0,
            "internal_stop_codon_percent": 0.0,
        },
    }


FULL_TABLE_FIXTURE = """# BUSCO version is: 6.0.0
# The lineage dataset is: solanales_odb12 (Creation date: 2025-07-01, number of genomes: 9, number of BUSCOs: 4)
# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\tOrthoDB url\tDescription
busco_a\tComplete\tchr1\t10\t90\t+\t201.0\t80\thttps://example.org/a\tprotein a
busco_b\tDuplicated\tchr2\t20\t120\t+\t301.0\t101\thttps://example.org/b\tprotein b
busco_b\tDuplicated\tchr3\t30\t140\t-\t299.5\t101\thttps://example.org/b\tprotein b copy
busco_c\tFragmented\tchr4\t40\t88\t+\t91.2\t48\thttps://example.org/c\tprotein c
busco_d\tMissing
"""


def write_fixture_busco_run(repo_root: Path, sample_id: str, taxon_id: str) -> None:
    stable_paths = busco_output_paths(sample_id)
    sample_root = repo_root / stable_paths["sample_root"]
    raw_root = repo_root / stable_paths["raw_root"]
    short_summary_path = repo_root / stable_paths["short_summary"]
    full_table_path = repo_root / stable_paths["full_table"]
    paths_tsv_path = repo_root / stable_paths["paths"]
    completion_path = repo_root / stable_paths["completion"]

    sequence_root = repo_root / stable_paths["sequence_root"]
    single_copy_dir = sequence_root / "single_copy_busco_sequences"
    multi_copy_dir = sequence_root / "multi_copy_busco_sequences"
    fragmented_dir = sequence_root / "fragmented_busco_sequences"

    for directory in (sample_root, raw_root, single_copy_dir, multi_copy_dir, fragmented_dir):
        directory.mkdir(parents=True, exist_ok=True)

    short_summary_path.write_text(
        json.dumps(
            make_short_summary_payload(
                sample_id=sample_id,
                input_fasta=f"assemblies/{taxon_id}.fa.gz",
            ),
            indent=2,
        ),
        encoding="utf-8",
    )
    full_table_path.write_text(FULL_TABLE_FIXTURE, encoding="utf-8")
    completion_path.write_text(f"sample_id\t{sample_id}\n", encoding="utf-8")

    for directory, busco_id in (
        (single_copy_dir, "busco_a"),
        (multi_copy_dir, "busco_b"),
        (fragmented_dir, "busco_c"),
    ):
        (directory / f"{busco_id}.faa").write_text(f">{busco_id}\nMA*\n", encoding="utf-8")
        (directory / f"{busco_id}.gff").write_text("##gff-version 3\n", encoding="utf-8")

    write_tsv(
        paths_tsv_path,
        [
            {"artifact": "raw_root", "path": stable_paths["raw_root"]},
            {"artifact": "short_summary_source", "path": stable_paths["short_summary"]},
            {"artifact": "full_table_source", "path": stable_paths["full_table"]},
            {
                "artifact": "single_copy_busco_sequences",
                "path": stable_paths["single_copy_sequence_dir"],
            },
            {
                "artifact": "multi_copy_busco_sequences",
                "path": stable_paths["multi_copy_sequence_dir"],
            },
            {
                "artifact": "fragmented_busco_sequences",
                "path": stable_paths["fragmented_sequence_dir"],
            },
        ],
        ["artifact", "path"],
    )


class BuscoSummaryParserTests(unittest.TestCase):
    def test_parse_busco_short_summary_fixture(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            summary_path = Path(tmpdir) / "short_summary.txt"
            summary_path.write_text(
                json.dumps(
                    make_short_summary_payload(
                        sample_id="sample_one",
                        input_fasta="assemblies/sample_one.fa.gz",
                    )
                ),
                encoding="utf-8",
            )
            parsed = parse_busco_short_summary(summary_path)

        self.assertEqual(parsed["lineage_name"], "solanales_odb12")
        self.assertEqual(parsed["complete_buscos"], 2)
        self.assertEqual(parsed["single_copy_buscos"], 1)
        self.assertEqual(parsed["multi_copy_buscos"], 1)
        self.assertEqual(parsed["missing_buscos"], 1)

    def test_parse_busco_full_table_fixture(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            table_path = Path(tmpdir) / "full_table.tsv"
            table_path.write_text(FULL_TABLE_FIXTURE, encoding="utf-8")
            rows = parse_busco_full_table(table_path)
            counts = summarize_busco_full_table(rows)

        self.assertEqual(len(rows), 5)
        self.assertEqual(rows[0]["busco_id"], "busco_a")
        self.assertEqual(rows[1]["status"], "Duplicated")
        self.assertEqual(rows[-1]["sequence_id"], "")
        self.assertEqual(counts["full_table_single_copy_buscos"], 1)
        self.assertEqual(counts["full_table_multi_copy_buscos"], 1)
        self.assertEqual(counts["full_table_fragmented_buscos"], 1)
        self.assertEqual(counts["full_table_missing_buscos"], 1)
        self.assertEqual(counts["full_table_duplicated_rows"], 2)

    def test_malformed_busco_output_fails_clearly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            summary_path = Path(tmpdir) / "short_summary.txt"
            summary_path.write_text('{"parameters": {}}', encoding="utf-8")
            with self.assertRaises(BuscoValidationError) as context:
                parse_busco_short_summary(summary_path)

        self.assertIn("lineage_dataset", str(context.exception))


class BuscoSummaryIntegrationTests(unittest.TestCase):
    def test_summary_cli_writes_one_row_per_sample(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            manifest_path = repo_root / "validated_manifest.tsv"
            summary_output = repo_root / "results/qc/busco_summary.tsv"
            records_output = repo_root / "results/qc/busco_records.tsv"

            write_tsv(
                manifest_path,
                [
                    {
                        "sample_id": "sample_one",
                        "taxon_id": "Sample one",
                        "sanitized_taxon_id": "sample_one",
                        "assembly_fasta": "assemblies/sample_one.fa.gz",
                    },
                    {
                        "sample_id": "sample_two",
                        "taxon_id": "Sample two",
                        "sanitized_taxon_id": "sample_two",
                        "assembly_fasta": "assemblies/sample_two.fa.gz",
                    },
                ],
                ["sample_id", "taxon_id", "sanitized_taxon_id", "assembly_fasta"],
            )
            write_fixture_busco_run(repo_root, "sample_one", "sample_one")
            write_fixture_busco_run(repo_root, "sample_two", "sample_two")

            subprocess.run(
                [
                    "python3",
                    "-m",
                    "scripts.summarize_busco",
                    "--manifest",
                    str(manifest_path),
                    "--summary-output",
                    str(summary_output),
                    "--records-output",
                    str(records_output),
                    "--repo-root",
                    str(repo_root),
                ],
                cwd=REPO_ROOT,
                check=True,
                capture_output=True,
                text=True,
            )

            with summary_output.open(newline="", encoding="utf-8") as handle:
                summary_rows = list(csv.DictReader(handle, delimiter="\t"))
            with records_output.open(newline="", encoding="utf-8") as handle:
                record_rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(len(summary_rows), 2)
        self.assertEqual(summary_rows[0]["sample_id"], "sample_one")
        self.assertEqual(summary_rows[1]["sample_id"], "sample_two")
        self.assertEqual(summary_rows[0]["full_table_consistent"], "true")
        self.assertEqual(len(record_rows), 10)


class BuscoSummaryWorkflowTests(unittest.TestCase):
    def test_summarize_busco_depends_on_completed_busco_outputs(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "-R",
                "summarize_busco",
                "--cores",
                "4",
                "results/qc/busco_summary.tsv",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("rule summarize_busco:", result.stdout)
        self.assertIn("results/busco/solanum_chilense/run.complete", result.stdout)
        self.assertIn("results/busco/solanum_chilense/short_summary.txt", result.stdout)
        self.assertIn("results/qc/busco_summary.tsv", result.stdout)


if __name__ == "__main__":
    unittest.main()
