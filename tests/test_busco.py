import subprocess
import tempfile
import unittest
from pathlib import Path

from scripts.busco import (
    BuscoValidationError,
    busco_output_paths,
    parse_dataset_names,
    verify_lineage_name,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


class BuscoUnitTests(unittest.TestCase):
    def test_parse_dataset_names_extracts_first_column(self):
        text = """
        2026-03-27 11:39:43 INFO:\tDownloading information on latest versions of BUSCO data...

        - viridiplantae_odb12 [822]
            - embryophyta_odb12 [2026]
                - solanales_odb12 [7934]
        """
        self.assertEqual(
            parse_dataset_names(text),
            ["viridiplantae_odb12", "embryophyta_odb12", "solanales_odb12"],
        )

    def test_verify_lineage_name_accepts_valid_lineage(self):
        lineage = verify_lineage_name(
            "solanales_odb12",
            ["embryophyta_odb12", "solanales_odb12", "viridiplantae_odb12"],
        )
        self.assertEqual(lineage, "solanales_odb12")

    def test_verify_lineage_name_rejects_invalid_lineage(self):
        with self.assertRaises(BuscoValidationError):
            verify_lineage_name(
                "solanales_odb10",
                ["embryophyta_odb12", "solanales_odb12", "viridiplantae_odb12"],
            )

    def test_verify_lineage_name_rejects_empty_lineage(self):
        with self.assertRaises(BuscoValidationError):
            verify_lineage_name("", ["solanales_odb12"])

    def test_busco_output_paths_are_derived_from_sample_id(self):
        paths = busco_output_paths("solanum_chilense")
        self.assertEqual(paths["sample_root"], "results/busco/solanum_chilense")
        self.assertEqual(paths["raw_root"], "work/busco/solanum_chilense/raw")
        self.assertEqual(
            paths["single_copy_sequence_dir"],
            "results/busco/solanum_chilense/busco_sequences/single_copy_busco_sequences",
        )
        self.assertEqual(paths["completion"], "results/busco/solanum_chilense/run.complete")
        self.assertEqual(paths["short_summary"], "results/busco/solanum_chilense/short_summary.txt")


class BuscoWorkflowTests(unittest.TestCase):
    def test_busco_rule_dry_run_renders_with_four_cpus(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "-R",
                "run_busco",
                "--cores",
                "4",
                "results/busco/solanum_chilense/run.complete",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("--cpu 4", result.stdout)
        self.assertIn("--lineage_dataset solanales_odb12", result.stdout)
        self.assertIn("results/busco/solanum_chilense/run.complete", result.stdout)


if __name__ == "__main__":
    unittest.main()
