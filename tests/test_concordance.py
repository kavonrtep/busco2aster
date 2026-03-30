import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.concordance import (
    build_iqtree_scfl_command,
    build_iqtree_gcf_command,
    concordance_output_paths,
    quartet_output_paths,
    summarize_freqquad,
    summarize_gcf_stat,
    summarize_scfl_stat,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
IQTREE_EXECUTABLE = REPO_ROOT / "work/tools/iqtree3/current/bin/iqtree3"


class ConcordanceUnitTests(unittest.TestCase):
    def test_concordance_output_paths_are_deterministic(self):
        paths = concordance_output_paths("gcf")
        self.assertEqual(paths["prefix"], "results/concordance/gcf")
        self.assertEqual(paths["command"], "results/concordance/gcf.command.sh")
        self.assertEqual(paths["stat"], "results/concordance/gcf.cf.stat")
        self.assertEqual(paths["tree"], "results/concordance/gcf.cf.tree")
        self.assertEqual(paths["branch"], "results/concordance/gcf.cf.branch")
        self.assertEqual(paths["log"], "results/concordance/gcf.log")
        self.assertEqual(paths["completion"], "results/concordance/gcf.complete")

    def test_quartet_output_paths_are_deterministic(self):
        paths = quartet_output_paths()
        self.assertEqual(paths["prefix"], "results/concordance/wastral_quartets")
        self.assertEqual(paths["tree"], "results/concordance/wastral_quartets.annotated.tre")
        self.assertEqual(paths["freqquad"], "results/concordance/wastral_quartets.freqquad.tsv")
        self.assertEqual(paths["log"], "results/concordance/wastral_quartets.log")
        self.assertEqual(paths["completion"], "results/concordance/wastral_quartets.complete")

    def test_build_iqtree_gcf_command_renders_expected_flags(self):
        command = build_iqtree_gcf_command(
            executable="iqtree3",
            reference_tree_path="results/species_tree/species_tree.wastral.tre",
            gene_tree_path="results/gene_trees/gene_trees.raw.tre",
            prefix="results/concordance/gcf",
            threads=4,
        )
        self.assertEqual(command[0], "iqtree3")
        self.assertIn("--gcf", command)
        self.assertIn("-t", command)
        self.assertIn("--prefix", command)
        self.assertIn("-T", command)
        self.assertIn("4", command)
        self.assertIn("--redo", command)
        self.assertIn("--quiet", command)

    def test_build_iqtree_scfl_command_renders_expected_flags(self):
        command = build_iqtree_scfl_command(
            executable="iqtree3",
            reference_tree_path="results/species_tree/species_tree.wastral.tre",
            alignment_dir="results/loci/alignments",
            prefix="results/concordance/scfl",
            threads=4,
            quartets=100,
            seqtype="AA",
            model="LG+G4",
        )
        self.assertEqual(command[0], "iqtree3")
        self.assertIn("--scfl", command)
        self.assertIn("--seqtype", command)
        self.assertIn("AA", command)
        self.assertIn("-m", command)
        self.assertIn("LG+G4", command)
        self.assertIn("-p", command)
        self.assertIn("-te", command)
        self.assertNotIn("-t", command)
        self.assertIn("results/loci/alignments", command)
        self.assertIn("100", command)
        self.assertIn("--prefix", command)
        self.assertIn("-T", command)

    def test_summarize_gcf_stat_parses_realistic_tabular_output(self):
        with TemporaryDirectory() as tmpdir:
            stat_path = Path(tmpdir) / "gcf.cf.stat"
            stat_path.write_text(
                "# Concordance factor statistics\n"
                "ID\tgCF\tgCF_N\tgDF1\tgDF1_N\tgDF2\tgDF2_N\tgDFP\tgDFP_N\tgN\tLabel\tLength\n"
                "6\t66.67\t2\t0\t0\t0\t0\t33.33\t1\t3\t\t-1\n"
                "7\t50.00\t1\t25.00\t1\t0\t0\t25.00\t1\t2\t\t-1\n",
                encoding="utf-8",
            )

            summary = summarize_gcf_stat(stat_path)

        self.assertEqual(summary["branch_count"], 2)
        self.assertAlmostEqual(summary["mean_gcf"], 58.335, places=3)
        self.assertAlmostEqual(summary["median_gcf"], 58.335, places=3)
        self.assertAlmostEqual(summary["min_gcf"], 50.0)
        self.assertAlmostEqual(summary["max_gcf"], 66.67)
        self.assertEqual(summary["lowest_rows"][0]["branch_id"], "7")

    def test_summarize_scfl_stat_parses_realistic_tabular_output(self):
        with TemporaryDirectory() as tmpdir:
            stat_path = Path(tmpdir) / "scfl.cf.stat"
            stat_path.write_text(
                "# Concordance factor statistics\n"
                "ID\tsCF\tsCF_N\tsDF1\tsDF1_N\tsDF2\tsDF2_N\tsN\tLabel\tLength\n"
                "5\t94.44\t17\t0\t0\t5.56\t1\t18\t\t0.5\n"
                "6\t40.00\t8\t30.00\t6\t30.00\t6\t20\t\t0.4\n",
                encoding="utf-8",
            )

            summary = summarize_scfl_stat(stat_path)

        self.assertEqual(summary["branch_count"], 2)
        self.assertAlmostEqual(summary["mean_scfl"], 67.22, places=2)
        self.assertAlmostEqual(summary["median_scfl"], 67.22, places=2)
        self.assertAlmostEqual(summary["min_scfl"], 40.0)
        self.assertAlmostEqual(summary["max_scfl"], 94.44)
        self.assertEqual(summary["lowest_rows"][0]["branch_id"], "6")

    def test_summarize_freqquad_parses_branch_level_quartet_support(self):
        with TemporaryDirectory() as tmpdir:
            freqquad_path = Path(tmpdir) / "wastral_quartets.freqquad.tsv"
            freqquad_path.write_text(
                "N1\tt1\t{d}|{e}#{c}|{a,b}\t0.952381\t3\t3\n"
                "N1\tt2\t{a,b}|{e}#{c}|{d}\t0.0238095\t0\t3\n"
                "N1\tt3\t{d}|{a,b}#{c}|{e}\t0.0238095\t0\t3\n"
                "N2\tt1\t{d,e}|{c}#{b}|{a}\t0.777777\t4\t10\n"
                "N2\tt2\t{a}|{c}#{b}|{d,e}\t0.111111\t3\t10\n"
                "N2\tt3\t{d,e}|{a}#{b}|{c}\t0.111111\t3\t10\n",
                encoding="utf-8",
            )

            summary = summarize_freqquad(freqquad_path)

        self.assertEqual(summary["node_count"], 2)
        self.assertEqual(summary["scored_node_count"], 2)
        self.assertAlmostEqual(summary["mean_best_frequency"], 0.7, places=6)
        self.assertAlmostEqual(summary["median_best_frequency"], 0.7, places=6)
        self.assertAlmostEqual(summary["min_best_frequency"], 0.4, places=6)
        self.assertAlmostEqual(summary["max_best_frequency"], 1.0, places=6)
        self.assertEqual(summary["lowest_rows"][0]["node_id"], "N2")


class ConcordanceSmokeTests(unittest.TestCase):
    @unittest.skipUnless(IQTREE_EXECUTABLE.is_file(), "IQ-TREE 3 binary is not installed under work/tools/iqtree3/current.")
    def test_iqtree_gcf_smoke_run_writes_cf_stat_outputs(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            species_tree_path = tmp_path / "species_tree.tre"
            gene_trees_path = tmp_path / "gene_trees.tre"
            prefix = tmp_path / "gcf"

            species_tree_path.write_text("((a,b),(c,d));\n", encoding="utf-8")
            gene_trees_path.write_text(
                "((a,b),(c,d));\n"
                "((a,b),(c,d));\n"
                "((a,c),(b,d));\n",
                encoding="utf-8",
            )

            command = build_iqtree_gcf_command(
                executable=IQTREE_EXECUTABLE.as_posix(),
                reference_tree_path=species_tree_path.as_posix(),
                gene_tree_path=gene_trees_path.as_posix(),
                prefix=prefix.as_posix(),
                threads=1,
            )
            subprocess.run(command, capture_output=True, text=True, check=True)
            summary = summarize_gcf_stat(prefix.with_suffix(".cf.stat"))

        self.assertGreaterEqual(summary["branch_count"], 1)
        self.assertGreaterEqual(summary["max_gcf"], summary["min_gcf"])

    @unittest.skipUnless(IQTREE_EXECUTABLE.is_file(), "IQ-TREE 3 binary is not installed under work/tools/iqtree3/current.")
    def test_iqtree_scfl_smoke_run_writes_cf_stat_outputs(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            species_tree_path = tmp_path / "species_tree.tre"
            alignment_dir = tmp_path / "alignments"
            prefix = tmp_path / "scfl"

            species_tree_path.write_text("((a,b),(c,d));\n", encoding="utf-8")
            alignment_dir.mkdir()
            (alignment_dir / "l1.faa").write_text(
                ">a\nMMMMMMMMMM\n>b\nMMMMMMMMML\n>c\nLLLLLLLLLL\n>d\nLLLLLLLLLM\n",
                encoding="utf-8",
            )
            (alignment_dir / "l2.faa").write_text(
                ">a\nMMMMMMMMMM\n>b\nMMMMMMMMLM\n>c\nLLLLLLLLLL\n>d\nLLLLLLLLLM\n",
                encoding="utf-8",
            )

            command = build_iqtree_scfl_command(
                executable=IQTREE_EXECUTABLE.as_posix(),
                reference_tree_path=species_tree_path.as_posix(),
                alignment_dir=alignment_dir.as_posix(),
                prefix=prefix.as_posix(),
                threads=1,
                quartets=100,
                seqtype="AA",
                model="LG+G4",
            )
            subprocess.run(command, capture_output=True, text=True, check=True)
            summary = summarize_scfl_stat(prefix.with_suffix(".cf.stat"))

        self.assertGreaterEqual(summary["branch_count"], 1)
        self.assertGreaterEqual(summary["max_scfl"], summary["min_scfl"])


class ConcordanceWorkflowTests(unittest.TestCase):
    def test_concordance_rule_dry_run_renders_expected_rule(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "--cores",
                "4",
                "results/concordance/gcf.cf.stat",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        if "Nothing to be done" not in result.stdout:
            self.assertIn("rule infer_gene_concordance:", result.stdout)
            self.assertIn("results/species_tree/species_tree.wastral.tre", result.stdout)
            self.assertIn("results/gene_trees/gene_trees.raw.tre", result.stdout)
        else:
            self.assertTrue((REPO_ROOT / "results/concordance/gcf.cf.stat").is_file())

    def test_scfl_rule_dry_run_renders_expected_rule(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "--cores",
                "4",
                "results/concordance/scfl.cf.stat",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        if "rule infer_site_concordance:" in result.stdout:
            self.assertIn("rule infer_site_concordance:", result.stdout)
            self.assertIn("results/species_tree/species_tree.wastral.tre", result.stdout)
            self.assertIn("results/loci/alignments.complete", result.stdout)
        else:
            self.assertTrue((REPO_ROOT / "results/concordance/scfl.cf.stat").is_file())

    def test_wastral_quartet_rule_dry_run_renders_expected_rule(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "--cores",
                "4",
                "results/concordance/wastral_quartets.freqquad.tsv",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        if "rule annotate_species_tree_quartets:" in result.stdout:
            self.assertIn("rule annotate_species_tree_quartets:", result.stdout)
            self.assertIn("results/species_tree/species_tree.wastral.tre", result.stdout)
            self.assertIn("results/gene_trees/gene_trees.wastral.tre", result.stdout)
        else:
            self.assertTrue((REPO_ROOT / "results/concordance/wastral_quartets.freqquad.tsv").is_file())

    @unittest.skipUnless((REPO_ROOT / "results/concordance/gcf.cf.stat").is_file(), "Real gCF output is not present.")
    def test_real_gcf_summary_has_scored_branches(self):
        summary = summarize_gcf_stat(REPO_ROOT / "results/concordance/gcf.cf.stat")
        self.assertGreater(summary["branch_count"], 0)
        self.assertGreaterEqual(summary["max_gcf"], summary["min_gcf"])

    @unittest.skipUnless((REPO_ROOT / "results/concordance/scfl.cf.stat").is_file(), "Real sCFL output is not present.")
    def test_real_scfl_summary_has_scored_branches(self):
        summary = summarize_scfl_stat(REPO_ROOT / "results/concordance/scfl.cf.stat")
        self.assertGreater(summary["branch_count"], 0)
        self.assertGreaterEqual(summary["max_scfl"], summary["min_scfl"])

    @unittest.skipUnless((REPO_ROOT / "results/concordance/wastral_quartets.freqquad.tsv").is_file(), "Real ASTER quartet output is not present.")
    def test_real_freqquad_summary_has_scored_branches(self):
        summary = summarize_freqquad(REPO_ROOT / "results/concordance/wastral_quartets.freqquad.tsv")
        self.assertGreater(summary["node_count"], 0)
        self.assertGreaterEqual(summary["max_best_frequency"], summary["min_best_frequency"])


if __name__ == "__main__":
    unittest.main()
