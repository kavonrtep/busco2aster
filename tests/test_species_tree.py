import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.gene_trees import read_single_line_tree
from scripts.species_tree import (
    build_aster_command,
    prepare_wastral_gene_tree_input,
    species_tree_output_paths,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
WASTRAL_EXECUTABLE = REPO_ROOT / "work/tools/aster/current/bin/wastral"
ASTRAL4_EXECUTABLE = REPO_ROOT / "work/tools/aster/current/bin/astral4"


class SpeciesTreeUnitTests(unittest.TestCase):
    def test_species_tree_output_paths_are_backend_specific(self):
        wastral_paths = species_tree_output_paths("wastral")
        astral4_paths = species_tree_output_paths("astral4")

        self.assertEqual(wastral_paths["treefile"], "results/species_tree/species_tree.wastral.tre")
        self.assertEqual(astral4_paths["treefile"], "results/species_tree/species_tree.astral4.tre")
        self.assertNotEqual(wastral_paths["command"], astral4_paths["command"])

    def test_build_aster_command_uses_backend_executable_and_threads(self):
        command = build_aster_command(
            backend="wastral",
            executable="wastral",
            input_path="results/gene_trees/gene_trees.raw.tre",
            output_path="results/species_tree/species_tree.wastral.tre",
            threads=4,
            support_mode="abayes",
        )
        self.assertEqual(command[0], "wastral")
        self.assertIn("-B", command)
        self.assertIn("-t", command)
        self.assertIn("4", command)
        self.assertIn("-i", command)
        self.assertIn("results/gene_trees/gene_trees.raw.tre", command)

    def test_prepare_wastral_gene_tree_input_rewrites_abayes_labels(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            input_path = tmp_path / "raw.tre"
            output_path = tmp_path / "wastral.tre"
            input_path.write_text(
                "(a:0.1,b:0.1,(c:0.1,(d:0.1,e:0.1)/1:0.1)/0.333:0.1);\n",
                encoding="utf-8",
            )

            prepare_wastral_gene_tree_input(input_path, output_path, support_mode="abayes")

            self.assertEqual(
                output_path.read_text(encoding="utf-8").strip(),
                "(a:0.1,b:0.1,(c:0.1,(d:0.1,e:0.1)1:0.1)0.333:0.1);",
            )


class SpeciesTreeSmokeTests(unittest.TestCase):
    @unittest.skipUnless(WASTRAL_EXECUTABLE.is_file(), "ASTER wastral binary is not installed under work/tools/aster/current.")
    def test_wastral_smoke_run_writes_species_tree(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            raw_gene_trees_path = tmp_path / "gene_trees.raw.tre"
            gene_trees_path = tmp_path / "gene_trees.wastral.tre"
            output_path = tmp_path / "species_tree.tre"
            raw_gene_trees_path.write_text(
                "(a:0.0308709501,b:0.0302158872,(c:0.0000010000,(d:0.0000010000,e:0.1309274391)/1:0.0991457827)/0.333:0.0000010000);\n"
                "(a:0.0308709501,b:0.0302158872,(c:0.0000010000,(d:0.0000010000,e:0.1309274391)/1:0.0991457827)/0.333:0.0000010000);\n"
                "(a:0.0308709501,b:0.0302158872,(c:0.0000010000,(d:0.0000010000,e:0.1309274391)/1:0.0991457827)/0.333:0.0000010000);\n",
                encoding="utf-8",
            )
            prepare_wastral_gene_tree_input(raw_gene_trees_path, gene_trees_path, support_mode="abayes")
            command = build_aster_command(
                backend="wastral",
                executable=WASTRAL_EXECUTABLE.as_posix(),
                input_path=gene_trees_path.as_posix(),
                output_path=output_path.as_posix(),
                threads=1,
                support_mode="abayes",
            )
            subprocess.run(command, capture_output=True, text=True, check=True)
            tree_text = read_single_line_tree(output_path)

        self.assertIn("a", tree_text)
        self.assertIn("b", tree_text)
        self.assertIn("c", tree_text)
        self.assertIn("d", tree_text)

    @unittest.skipUnless(ASTRAL4_EXECUTABLE.is_file(), "ASTER astral4 binary is not installed under work/tools/aster/current.")
    def test_astral4_smoke_run_writes_species_tree(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            gene_trees_path = tmp_path / "gene_trees.tre"
            output_path = tmp_path / "species_tree.tre"
            gene_trees_path.write_text(
                "((a,b),(c,d));\n"
                "((a,c),(b,d));\n"
                "((a,b),(c,d));\n",
                encoding="utf-8",
            )
            command = build_aster_command(
                backend="astral4",
                executable=ASTRAL4_EXECUTABLE.as_posix(),
                input_path=gene_trees_path.as_posix(),
                output_path=output_path.as_posix(),
                threads=1,
            )
            subprocess.run(command, capture_output=True, text=True, check=True)
            tree_text = read_single_line_tree(output_path)

        self.assertIn("a", tree_text)
        self.assertIn("b", tree_text)
        self.assertIn("c", tree_text)
        self.assertIn("d", tree_text)


class SpeciesTreeWorkflowTests(unittest.TestCase):
    def test_species_tree_dry_run_defaults_to_wastral(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "--cores",
                "4",
                "results/species_tree/species_tree.complete",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertNotIn("rule infer_species_tree_astral4:", result.stdout)
        self.assertNotIn("rule align_locus:", result.stdout)
        if "Nothing to be done" not in result.stdout:
            self.assertIn("rule prepare_wastral_gene_trees:", result.stdout)
            self.assertIn("rule infer_species_tree_wastral:", result.stdout)

    @unittest.skipUnless((REPO_ROOT / "results/species_tree/species_tree.wastral.tre").is_file(), "Real wastral output is not present.")
    def test_real_species_tree_exists_after_workflow_run(self):
        tree_text = read_single_line_tree(REPO_ROOT / "results/species_tree/species_tree.wastral.tre")
        self.assertTrue(tree_text.endswith(";"))


if __name__ == "__main__":
    unittest.main()
