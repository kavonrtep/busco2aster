import csv
import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.gene_trees import (
    aggregate_directory_mode_gene_tree_outputs,
    build_iqtree_command,
    build_iqtree_directory_command,
    gene_tree_directory_output_paths,
    gene_tree_output_paths,
    load_gene_tree_locus_ids,
    parse_directory_mode_selected_models,
    read_single_line_tree,
    tree_has_support_values,
)
from scripts.manifest import write_tsv


REPO_ROOT = Path(__file__).resolve().parents[1]
IQTREE_EXECUTABLE = REPO_ROOT / "work/tools/iqtree3/current/bin/iqtree3"


class GeneTreeUnitTests(unittest.TestCase):
    def test_gene_tree_output_paths_are_deterministic(self):
        paths = gene_tree_output_paths("35at4069")
        self.assertEqual(paths["command"], "results/gene_trees/per_locus/35at4069/command.sh")
        self.assertEqual(paths["report"], "results/gene_trees/per_locus/35at4069/35at4069.iqtree")
        self.assertEqual(paths["log"], "results/gene_trees/per_locus/35at4069/35at4069.log")
        self.assertEqual(paths["treefile"], "results/gene_trees/per_locus/35at4069/35at4069.treefile")
        self.assertEqual(paths["completion"], "results/gene_trees/per_locus/35at4069/run.complete")

    def test_gene_tree_directory_output_paths_are_deterministic(self):
        paths = gene_tree_directory_output_paths()
        self.assertEqual(paths["prefix"], "results/gene_trees/gene_trees")
        self.assertEqual(paths["command"], "results/gene_trees/gene_trees.command.sh")
        self.assertEqual(paths["report"], "results/gene_trees/gene_trees.iqtree")
        self.assertEqual(paths["treefile"], "results/gene_trees/gene_trees.treefile")
        self.assertEqual(paths["best_model_nex"], "results/gene_trees/gene_trees.best_model.nex")

    def test_build_iqtree_command_uses_abayes_support_mode(self):
        command = build_iqtree_command(
            executable="iqtree3",
            alignment_path="results/loci/alignments/35at4069.aln.faa",
            prefix="results/gene_trees/per_locus/35at4069/35at4069",
            threads=1,
            model="MFP",
            support_mode="abayes",
            seed=20260327,
            ufboot_replicates=1000,
        )
        self.assertIn("--abayes", command)
        self.assertIn("--prefix", command)
        self.assertIn("results/gene_trees/per_locus/35at4069/35at4069", command)
        self.assertIn("-m", command)
        self.assertIn("MFP", command)

    def test_build_iqtree_directory_command_uses_separate_tree_mode(self):
        command = build_iqtree_directory_command(
            executable="iqtree3",
            alignment_dir="results/loci/alignments",
            prefix="results/gene_trees/gene_trees",
            threads=4,
            model="MFP",
            support_mode="abayes",
            seed=20260327,
            ufboot_replicates=1000,
        )
        self.assertIn("-S", command)
        self.assertIn("results/loci/alignments", command)
        self.assertIn("--abayes", command)
        self.assertIn("-T", command)
        self.assertIn("4", command)

    def test_parse_directory_mode_selected_models_reads_best_model_nexus(self):
        with TemporaryDirectory() as tmpdir:
            best_model_path = Path(tmpdir) / "loci.best_model.nex"
            best_model_path.write_text(
                "#nexus\n"
                "begin sets;\n"
                "  charpartition mymodels =\n"
                "    FLU: 15853at4069.aln.faa{0.1},\n"
                "    Q.BIRD+I{0.925262}: 10235at4069.aln.faa{0.2};\n"
                "end;\n",
                encoding="utf-8",
            )

            entries = parse_directory_mode_selected_models(best_model_path)

        self.assertEqual(entries, [("15853at4069", "FLU"), ("10235at4069", "Q.BIRD+I")])

    def test_aggregate_directory_mode_gene_tree_outputs_preserves_all_loci_once(self):
        with TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            manifest_path = repo_root / "gene_tree_manifest.tsv"
            aggregate_path = repo_root / "gene_trees.raw.tre"
            report_path = repo_root / "gene_trees.iqtree"
            treefile_path = repo_root / "gene_trees.treefile"
            best_model_path = repo_root / "gene_trees.best_model.nex"

            report_path.write_text("Best-fit model according to BIC: LG+F:2at4069.aln.faa,WAG+F:10at4069.aln.faa\n", encoding="utf-8")
            treefile_path.write_text(
                "(a:0.1,(b:0.1,c:0.1)/0.9:0.1);\n"
                "(d:0.1,(e:0.1,f:0.1)/0.8:0.1);\n",
                encoding="utf-8",
            )
            best_model_path.write_text(
                "#nexus\n"
                "begin sets;\n"
                "  charpartition mymodels =\n"
                "    LG+F: 2at4069.aln.faa{0.013},\n"
                "    WAG+F: 10at4069.aln.faa{0.023};\n"
                "end;\n",
                encoding="utf-8",
            )

            aggregate_directory_mode_gene_tree_outputs(
                best_model_nex_path=best_model_path,
                treefile_path=treefile_path,
                report_path=report_path,
                manifest_path=manifest_path,
                aggregate_path=aggregate_path,
                support_mode="abayes",
            )

            self.assertEqual(load_gene_tree_locus_ids(manifest_path), ["2at4069", "10at4069"])
            with aggregate_path.open(encoding="utf-8") as handle:
                aggregate_lines = [line.strip() for line in handle if line.strip()]
            self.assertEqual(len(aggregate_lines), 2)
            self.assertTrue(aggregate_lines[0].startswith("(a:"))
            self.assertTrue(aggregate_lines[1].startswith("(d:"))


class GeneTreeSmokeTests(unittest.TestCase):
    @unittest.skipUnless(IQTREE_EXECUTABLE.is_file(), "IQ-TREE 3 binary is not installed under work/tools/iqtree3/current.")
    def test_iqtree_smoke_run_writes_support_annotated_tree(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            alignment_path = tmp_path / "input.faa"
            prefix = tmp_path / "run"
            alignment_path.write_text(
                ">a\nMKTAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCC\n"
                ">b\nMKTAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCA\n"
                ">c\nMKTAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCG\n"
                ">d\nMKTAAAGGGGTTTTCCCCAAAAGGGGTTTTGGGG\n"
                ">e\nMKTAAAGGGGTTTTCCCCAAAAGGGGAAAAGGGG\n",
                encoding="utf-8",
            )
            command = build_iqtree_command(
                executable=IQTREE_EXECUTABLE.as_posix(),
                alignment_path=alignment_path.as_posix(),
                prefix=prefix.as_posix(),
                threads=1,
                model="MFP",
                support_mode="abayes",
                seed=20260327,
                ufboot_replicates=1000,
            )
            subprocess.run(command, capture_output=True, text=True, check=True)
            tree_text = read_single_line_tree(prefix.with_suffix(".treefile"))

        self.assertTrue(tree_has_support_values(tree_text, "abayes"))

    @unittest.skipUnless(IQTREE_EXECUTABLE.is_file(), "IQ-TREE 3 binary is not installed under work/tools/iqtree3/current.")
    def test_iqtree_directory_mode_smoke_run_writes_multi_tree_output(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            alignment_dir = tmp_path / "alignments"
            alignment_dir.mkdir()
            (alignment_dir / "l1.aln.faa").write_text(
                ">a\nMMMMMMMMMM\n>b\nMMMMMMMMML\n>c\nLLLLLLLLLL\n>d\nLLLLLLLLLM\n>e\nLLLLLLLLLA\n",
                encoding="utf-8",
            )
            (alignment_dir / "l2.aln.faa").write_text(
                ">a\nMMMMMMMMMM\n>b\nMMMMMMMMLM\n>c\nLLLLLLLLLL\n>d\nLLLLLLLLLM\n>e\nLLLLLLLALL\n",
                encoding="utf-8",
            )
            prefix = tmp_path / "loci"

            command = build_iqtree_directory_command(
                executable=IQTREE_EXECUTABLE.as_posix(),
                alignment_dir=alignment_dir.as_posix(),
                prefix=prefix.as_posix(),
                threads=1,
                model="MFP",
                support_mode="abayes",
                seed=20260327,
                ufboot_replicates=1000,
            )
            subprocess.run(command, capture_output=True, text=True, check=True)
            tree_lines = [line.strip() for line in prefix.with_suffix(".treefile").read_text(encoding="utf-8").splitlines() if line.strip()]

        self.assertEqual(len(tree_lines), 2)
        self.assertTrue(all(tree_has_support_values(tree, "abayes") for tree in tree_lines))


class GeneTreeWorkflowTests(unittest.TestCase):
    def test_gene_tree_rule_dry_run_renders_directory_mode_targets(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "-R",
                "infer_gene_trees",
                "aggregate_gene_trees",
                "--cores",
                "4",
                "results/gene_trees/gene_trees.raw.tre",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        if "Nothing to be done" not in result.stdout:
            self.assertIn("rule infer_gene_trees:", result.stdout)
            self.assertIn("rule aggregate_gene_trees:", result.stdout)
            self.assertIn("results/gene_trees/gene_trees.treefile", result.stdout)
            self.assertIn("results/loci/alignments.complete", result.stdout)
        else:
            self.assertTrue((REPO_ROOT / "results/gene_trees/gene_trees.raw.tre").is_file())

    @unittest.skipUnless((REPO_ROOT / "results/gene_trees/gene_tree_manifest.tsv").is_file(), "Real gene-tree outputs are not present.")
    def test_real_gene_tree_manifest_matches_retained_locus_count(self):
        with (REPO_ROOT / "results/qc/retained_loci.tsv").open(newline="", encoding="utf-8") as handle:
            retained_rows = list(csv.DictReader(handle, delimiter="\t"))
        retained_ids = [row["locus_id"] for row in retained_rows if row["decision"] == "retain"]

        with (REPO_ROOT / "results/gene_trees/gene_tree_manifest.tsv").open(
            newline="", encoding="utf-8"
        ) as handle:
            manifest_rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(len(manifest_rows), len(retained_ids))
        self.assertEqual({row["locus_id"] for row in manifest_rows}, set(retained_ids))
        self.assertEqual(len({row["locus_id"] for row in manifest_rows}), len(manifest_rows))
        self.assertEqual(
            [int(row["tree_row_index"]) for row in manifest_rows],
            list(range(1, len(manifest_rows) + 1)),
        )


if __name__ == "__main__":
    unittest.main()
