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


if __name__ == "__main__":
    unittest.main()
