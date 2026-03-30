import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

from scripts.manifest import (
    ManifestValidationError,
    load_samples_manifest,
    sanitize_taxon_label,
    validate_manifest_rows,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


class ManifestUnitTests(unittest.TestCase):
    def write_manifest(self, directory: Path, header: str, rows: list[str]) -> Path:
        manifest_path = directory / "samples.tsv"
        manifest_path.write_text("\n".join([header, *rows]) + "\n", encoding="utf-8")
        return manifest_path

    def test_valid_manifest_is_parsed_and_validated(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / "sample.fna.gz"
            fasta.write_text("placeholder", encoding="utf-8")
            manifest = self.write_manifest(
                root,
                "sample_id\ttaxon_id\tassembly_fasta",
                ["sample_a\tSolanum_chilense\tsample.fna.gz"],
            )

            rows = load_samples_manifest(manifest)
            validated_rows, taxon_rows = validate_manifest_rows(rows, repo_root=root)

            self.assertEqual(len(validated_rows), 1)
            self.assertEqual(validated_rows[0]["sample_id"], "sample_a")
            self.assertEqual(validated_rows[0]["sanitized_taxon_id"], "solanum_chilense")
            self.assertEqual(taxon_rows[0]["sanitized_taxon_id"], "solanum_chilense")

    def test_missing_required_columns_fail(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest = self.write_manifest(
                root,
                "sample_id\ttaxon_id",
                ["sample_a\tSolanum_chilense"],
            )
            with self.assertRaises(ManifestValidationError):
                load_samples_manifest(manifest)

    def test_duplicate_sample_ids_fail(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta_a = root / "a.fna.gz"
            fasta_b = root / "b.fna.gz"
            fasta_a.write_text("a", encoding="utf-8")
            fasta_b.write_text("b", encoding="utf-8")
            manifest = self.write_manifest(
                root,
                "sample_id\ttaxon_id\tassembly_fasta",
                [
                    "dup\tSolanum_chilense\ta.fna.gz",
                    "dup\tSolanum_pennellii\tb.fna.gz",
                ],
            )
            rows = load_samples_manifest(manifest)
            with self.assertRaises(ManifestValidationError):
                validate_manifest_rows(rows, repo_root=root)

    def test_duplicate_taxon_ids_fail(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta_a = root / "a.fna.gz"
            fasta_b = root / "b.fna.gz"
            fasta_a.write_text("a", encoding="utf-8")
            fasta_b.write_text("b", encoding="utf-8")
            manifest = self.write_manifest(
                root,
                "sample_id\ttaxon_id\tassembly_fasta",
                [
                    "sample_a\tSolanum_chilense\ta.fna.gz",
                    "sample_b\tSolanum_chilense\tb.fna.gz",
                ],
            )
            rows = load_samples_manifest(manifest)
            with self.assertRaises(ManifestValidationError):
                validate_manifest_rows(rows, repo_root=root)

    def test_nonexistent_assembly_path_fails(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest = self.write_manifest(
                root,
                "sample_id\ttaxon_id\tassembly_fasta",
                ["sample_a\tSolanum_chilense\tmissing.fna.gz"],
            )
            rows = load_samples_manifest(manifest)
            with self.assertRaises(ManifestValidationError):
                validate_manifest_rows(rows, repo_root=root)

    def test_taxon_sanitization_is_deterministic(self):
        taxon = "Solanum_lycopersicum_var._cerasiforme"
        first = sanitize_taxon_label(taxon)
        second = sanitize_taxon_label(taxon)
        self.assertEqual(first, second)
        self.assertEqual(first, "solanum_lycopersicum_var_cerasiforme")


class ManifestWorkflowTests(unittest.TestCase):
    def test_validate_manifest_target_produces_metadata_files(self):
        validated = REPO_ROOT / "results" / "metadata" / "samples.validated.tsv"
        taxon_map = REPO_ROOT / "results" / "metadata" / "taxon_name_map.tsv"
        validated_target = "results/metadata/samples.validated.tsv"
        taxon_map_target = "results/metadata/taxon_name_map.tsv"

        if validated.exists():
            validated.unlink()
        if taxon_map.exists():
            taxon_map.unlink()

        subprocess.run(
            [
                "snakemake",
                "--nolock",
                "--cores",
                "1",
                validated_target,
                taxon_map_target,
            ],
            cwd=REPO_ROOT,
            check=True,
        )

        self.assertTrue(validated.is_file())
        self.assertTrue(taxon_map.is_file())
        self.assertIn("sanitized_taxon_id", validated.read_text(encoding="utf-8").splitlines()[0])
        self.assertIn("sanitized_taxon_id", taxon_map.read_text(encoding="utf-8").splitlines()[0])

    def test_adding_sample_reuses_existing_busco_and_rebuilds_downstream_aggregates(self):
        source_manifest = REPO_ROOT / "config" / "samples.tsv"
        with source_manifest.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertGreaterEqual(len(rows), 1)
        extra_row = {
            "sample_id": "solanum_demo_added",
            "taxon_id": "Solanum demo added",
            "assembly_fasta": rows[0]["assembly_fasta"],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            manifest_path = tmp_path / "samples.tsv"
            config_path = tmp_path / "config.yaml"

            with manifest_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    delimiter="\t",
                    fieldnames=["sample_id", "taxon_id", "assembly_fasta"],
                )
                writer.writeheader()
                writer.writerows([*rows, extra_row])

            config_path.write_text(
                f"samples: {manifest_path.as_posix()}\n"
                "busco_lineage: solanales_odb12\n"
                "threads:\n"
                "  default: 4\n"
                "  alignment: 1\n"
                "  busco: 4\n"
                "  concordance: 4\n"
                "  iqtree: 4\n"
                "  species_tree: 4\n",
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    "snakemake",
                    "-n",
                    "-p",
                    "--configfile",
                    config_path.as_posix(),
                    "--cores",
                    "4",
                    "results/report/report.html",
                ],
                cwd=REPO_ROOT,
                capture_output=True,
                text=True,
                check=True,
            )

        self.assertGreaterEqual(result.stdout.count("rule run_busco:"), 1)
        self.assertIn("wildcards: sample=solanum_demo_added", result.stdout)
        self.assertIn("rule summarize_busco:", result.stdout)
        self.assertIn("rule build_locus_matrix:", result.stdout)
        self.assertIn("checkpoint select_loci:", result.stdout)
        self.assertIn("rule export_retained_fastas:", result.stdout)
        self.assertIn("rule infer_gene_trees:", result.stdout)
        self.assertIn("rule infer_species_tree_wastral:", result.stdout)
        self.assertIn("rule infer_gene_concordance:", result.stdout)
        self.assertIn("rule infer_site_concordance:", result.stdout)
        self.assertIn("rule render_visual_report:", result.stdout)


if __name__ == "__main__":
    unittest.main()
