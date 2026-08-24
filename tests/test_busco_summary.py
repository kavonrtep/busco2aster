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


LOCUS_BY_CATEGORY = {
    "single_copy": "busco_a",
    "multi_copy": "busco_b",
    "fragmented": "busco_c",
}


def fixture_locus_sequences(busco_id: str) -> tuple[str, str]:
    faa = f">{busco_id}\nMA*\n"
    gff = f"##gff-version 3\nchr1\tbusco\tmRNA\t10\t90\t.\t+\t.\tID={busco_id}\n"
    return faa, gff


def write_fixture_busco_run(
    repo_root: Path,
    sample_id: str,
    taxon_id: str,
    layout: str = "concatenated",
) -> None:
    """Write one sample's BUSCO outputs in either sequence layout.

    `per_locus` is what releases before the concatenated layout produced, and is
    still read so that existing `results/busco/` trees keep working.
    """
    stable_paths = busco_output_paths(sample_id)
    sample_root = repo_root / stable_paths["sample_root"]
    raw_root = repo_root / stable_paths["raw_root"]
    short_summary_path = repo_root / stable_paths["short_summary"]
    full_table_path = repo_root / stable_paths["full_table"]
    paths_tsv_path = repo_root / stable_paths["paths"]
    completion_path = repo_root / stable_paths["completion"]
    sequence_root = repo_root / stable_paths["sequence_root"]

    for directory in (sample_root, raw_root, sequence_root):
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

    path_rows = [
        {"artifact": "raw_root", "path": stable_paths["raw_root"]},
        {"artifact": "short_summary_source", "path": stable_paths["short_summary"]},
        {"artifact": "full_table_source", "path": stable_paths["full_table"]},
        {"artifact": "sequence_layout", "path": layout},
    ]

    for category, busco_id in LOCUS_BY_CATEGORY.items():
        faa_text, gff_text = fixture_locus_sequences(busco_id)
        if layout == "concatenated":
            (repo_root / stable_paths[f"{category}_faa"]).write_text(faa_text, encoding="utf-8")
            (repo_root / stable_paths[f"{category}_gff"]).write_text(
                f"##busco_id {busco_id}\n{gff_text}", encoding="utf-8"
            )
            path_rows.append(
                {
                    "artifact": f"{category}_busco_sequences",
                    "path": stable_paths[f"{category}_faa"],
                }
            )
            path_rows.append(
                {"artifact": f"{category}_busco_gff", "path": stable_paths[f"{category}_gff"]}
            )
            continue

        per_locus_dir = repo_root / stable_paths[f"{category}_per_locus_dir"]
        per_locus_dir.mkdir(parents=True, exist_ok=True)
        (per_locus_dir / f"{busco_id}.faa").write_text(faa_text, encoding="utf-8")
        (per_locus_dir / f"{busco_id}.gff").write_text(gff_text, encoding="utf-8")
        path_rows.append(
            {
                "artifact": f"{category}_busco_sequences",
                "path": stable_paths[f"{category}_per_locus_dir"],
            }
        )

    write_tsv(paths_tsv_path, path_rows, ["artifact", "path"])


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
    def _run_summary_cli(self, layout: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
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
            write_fixture_busco_run(repo_root, "sample_one", "sample_one", layout=layout)
            write_fixture_busco_run(repo_root, "sample_two", "sample_two", layout=layout)

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
        return summary_rows, record_rows

    def test_summary_cli_writes_one_row_per_sample(self):
        summary_rows, record_rows = self._run_summary_cli("concatenated")

        self.assertEqual(len(summary_rows), 2)
        self.assertEqual(summary_rows[0]["sample_id"], "sample_one")
        self.assertEqual(summary_rows[1]["sample_id"], "sample_two")
        self.assertEqual(summary_rows[0]["full_table_consistent"], "true")
        self.assertEqual(len(record_rows), 10)
        self.assertEqual(summary_rows[0]["sequence_layout"], "concatenated")
        self.assertTrue(
            summary_rows[0]["single_copy_sequence_path"].endswith("busco_sequences/single_copy.faa")
        )

    def test_summary_cli_still_reads_per_locus_runs(self):
        """An existing results/busco/ tree keeps working without re-running BUSCO."""
        legacy_summary, legacy_records = self._run_summary_cli("per_locus")
        current_summary, current_records = self._run_summary_cli("concatenated")

        self.assertEqual(legacy_summary[0]["sequence_layout"], "per_locus")
        self.assertEqual(len(legacy_records), len(current_records))
        self.assertEqual(
            [(row["busco_id"], row["status"], row["sequence_category"]) for row in legacy_records],
            [(row["busco_id"], row["status"], row["sequence_category"]) for row in current_records],
        )
        self.assertTrue(
            legacy_records[0]["faa_path"].endswith("single_copy_busco_sequences/busco_a.faa")
        )


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
