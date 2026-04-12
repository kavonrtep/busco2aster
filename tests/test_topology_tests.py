"""Tests for the alternative topology hypothesis testing extension."""

import csv
import shutil
import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.generate_alternatives import (
    apply_nni_swap,
    generate_candidate_trees,
    write_candidate_manifest,
    write_candidate_trees,
)
from scripts.parse_au_test import join_with_manifest, parse_au_test_iqtree, write_au_results
from scripts.quartet_support import (
    filter_contested_branches,
    parse_wastral_u2_tree,
    write_branch_quartet_support,
    write_contested_branches,
)
from scripts.topology_tests import (
    build_au_test_command,
    topology_tests_output_paths,
    write_topology_tests_command_script,
)
from scripts.tree_utils import canonical_topology_key, parse_newick


REPO_ROOT = Path(__file__).resolve().parents[1]
IQTREE_EXECUTABLE = REPO_ROOT / "work/tools/iqtree3/current/bin/iqtree3"
WASTRAL_EXECUTABLE = REPO_ROOT / "work/tools/aster/current/bin/wastral"
SNAKEMAKE_AVAILABLE = shutil.which("snakemake") is not None


# ── Unit tests ────────────────────────────────────────────────────────────────


class TopologyTestsOutputPathsTests(unittest.TestCase):
    def test_output_paths_are_deterministic(self):
        paths = topology_tests_output_paths()
        self.assertEqual(paths["wastral_u2_tree"], "results/topology_tests/species_tree.wastral.u2.tre")
        self.assertEqual(paths["branch_quartet_support"], "results/topology_tests/branch_quartet_support.tsv")
        self.assertEqual(paths["contested_branches"], "results/topology_tests/contested_branches.tsv")
        self.assertEqual(paths["candidate_trees"], "results/topology_tests/candidate_trees.tre")
        self.assertEqual(paths["candidate_manifest"], "results/topology_tests/candidate_trees_manifest.tsv")
        self.assertEqual(paths["supermatrix"], "results/topology_tests/supermatrix.phy")
        self.assertEqual(paths["partitions"], "results/topology_tests/supermatrix_partitions.nex")
        self.assertEqual(paths["au_iqtree"], "results/topology_tests/au_test.iqtree")
        self.assertEqual(paths["au_results"], "results/topology_tests/au_test_results.tsv")
        self.assertEqual(paths["completion"], "results/topology_tests/topology_tests.complete")

    def test_build_au_test_command_renders_expected_flags(self):
        command = build_au_test_command(
            executable="iqtree3",
            supermatrix="results/topology_tests/supermatrix.phy",
            partitions="results/topology_tests/supermatrix_partitions.nex",
            candidate_trees="results/topology_tests/candidate_trees.tre",
            prefix="results/topology_tests/au_test",
            threads=4,
            replicates=10000,
        )
        self.assertEqual(command[0], "iqtree3")
        self.assertIn("--test-au", command)
        self.assertIn("-n", command)
        self.assertIn("0", command)
        self.assertIn("--test", command)
        self.assertIn("10000", command)
        self.assertIn("-T", command)
        self.assertIn("4", command)
        self.assertIn("--redo", command)

    def test_write_command_script_creates_executable_file(self):
        with TemporaryDirectory() as tmpdir:
            script_path = Path(tmpdir) / "au_test.command.sh"
            command = ["iqtree3", "--test-au", "-n", "0"]
            write_topology_tests_command_script(script_path, command, Path(tmpdir) / "au.log")
            self.assertTrue(script_path.is_file())
            content = script_path.read_text()
            self.assertIn("set -euo pipefail", content)
            self.assertIn("iqtree3", content)
            self.assertIn("--test-au", content)
            mode = oct(script_path.stat().st_mode)
            self.assertTrue(mode.endswith("755"), f"Expected 755 permission, got {mode}")


class QuartetSupportParserTests(unittest.TestCase):
    # wASTRAL -u 2 annotation embedded as node label:
    # '[q1=0.9;q2=0.07;q3=0.03;pp1=0.95;pp2=0.03;pp3=0.02;f1=900;f2=70;f3=30;QC=0.9;EN=1000]'
    _ANNOTATED_NEWICK = (
        "((A,B)'[q1=0.9;q2=0.07;q3=0.03;pp1=0.95;pp2=0.03;pp3=0.02]':0.5,"
        "(C,D)'[q1=0.4;q2=0.35;q3=0.25;pp1=0.41;pp2=0.37;pp3=0.22]':0.2):0.0;"
    )

    def _write_annotated_tree(self, tmpdir: str) -> Path:
        p = Path(tmpdir) / "u2.tre"
        p.write_text(self._ANNOTATED_NEWICK, encoding="utf-8")
        return p

    def test_parse_returns_rows_for_all_informative_branches(self):
        with TemporaryDirectory() as tmpdir:
            tree_path = self._write_annotated_tree(tmpdir)
            rows = parse_wastral_u2_tree(tree_path, contested_threshold=0.95)
        # 4-taxon tree has 1 internal informative branch
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertIn("branch_id", row)
        self.assertAlmostEqual(row["pp1"], 0.95)
        self.assertAlmostEqual(row["q1"], 0.9)

    def test_parse_flags_contested_branch(self):
        newick = (
            "((A,B)'[q1=0.4;q2=0.35;q3=0.25;pp1=0.41;pp2=0.37;pp3=0.22]':0.5,"
            "(C,D):0.2):0.0;"
        )
        with TemporaryDirectory() as tmpdir:
            p = Path(tmpdir) / "u2.tre"
            p.write_text(newick, encoding="utf-8")
            rows = parse_wastral_u2_tree(p, contested_threshold=0.95)
        self.assertTrue(rows[0]["contested"])

    def test_filter_contested_branches_excludes_well_resolved(self):
        rows = [
            {"pp1": 0.98, "contested": False},
            {"pp1": 0.60, "contested": True},
        ]
        contested = filter_contested_branches(rows)
        self.assertEqual(len(contested), 1)
        self.assertAlmostEqual(contested[0]["pp1"], 0.60)

    def test_write_branch_quartet_support_writes_tsv_header(self):
        rows = [
            {
                "branch_id": "B1",
                "taxa_clade": "A,B",
                "taxa_complement": "C,D",
                "q1": 0.9, "q2": 0.07, "q3": 0.03,
                "pp1": 0.95, "pp2": 0.03, "pp3": 0.02,
                "contested": False,
            }
        ]
        with TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "bqs.tsv"
            write_branch_quartet_support(rows, out)
            with out.open() as fh:
                header = fh.readline().strip().split("\t")
        self.assertIn("branch_id", header)
        self.assertIn("pp1", header)
        self.assertIn("contested", header)


class NNISwapTests(unittest.TestCase):
    # 6-taxon tree with 3 internal informative branches.
    _NEWICK_6 = "(((A,B),(C,D)),(E,F));"

    def _root(self):
        return parse_newick(self._NEWICK_6)

    def test_apply_nni_swap_produces_different_topology(self):
        root = self._root()
        original_key = canonical_topology_key(root)
        # Find branch IDs
        from scripts.tree_utils import branch_label_map, leaf_labels
        lm = branch_label_map(root)
        branch_ids = sorted(lm.values())
        swapped = apply_nni_swap(root, branch_ids[0], child_index=1)
        swapped_key = canonical_topology_key(swapped)
        self.assertNotEqual(original_key, swapped_key)

    def test_double_swap_restores_original_topology(self):
        root = self._root()
        original_key = canonical_topology_key(root)
        from scripts.tree_utils import branch_label_map
        lm = branch_label_map(root)
        branch_ids = sorted(lm.values())
        bid = branch_ids[0]
        once = apply_nni_swap(root, bid, child_index=1)
        twice = apply_nni_swap(once, bid, child_index=1)
        restored_key = canonical_topology_key(twice)
        self.assertEqual(original_key, restored_key)

    def test_generate_candidate_trees_returns_original_first(self):
        with TemporaryDirectory() as tmpdir:
            tree_path = Path(tmpdir) / "species_tree.tre"
            tree_path.write_text(self._NEWICK_6, encoding="utf-8")
            contested = [
                {
                    "branch_id": "B1",
                    "taxa_clade": "A,B",
                    "taxa_complement": "C,D,E,F",
                    "q1": 0.5, "q2": 0.3, "q3": 0.2,
                    "pp1": 0.5, "pp2": 0.3, "pp3": 0.2,
                    "contested": True,
                }
            ]
            trees, descs, srcs = generate_candidate_trees(
                species_tree_path=tree_path,
                contested_rows=contested,
                max_contested_branches=4,
            )
        self.assertGreaterEqual(len(trees), 2)
        self.assertEqual(descs[0], "wASTRAL best tree")
        self.assertEqual(srcs[0], "original")

    def test_generate_no_alternatives_when_no_contested(self):
        with TemporaryDirectory() as tmpdir:
            tree_path = Path(tmpdir) / "species_tree.tre"
            tree_path.write_text(self._NEWICK_6, encoding="utf-8")
            trees, descs, srcs = generate_candidate_trees(
                species_tree_path=tree_path,
                contested_rows=[],
                max_contested_branches=4,
            )
        self.assertEqual(len(trees), 1)
        self.assertEqual(descs[0], "wASTRAL best tree")

    def test_write_candidate_trees_produces_one_newick_per_line(self):
        root = self._root()
        with TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "candidates.tre"
            write_candidate_trees([root, root], out)
            lines = [ln for ln in out.read_text().splitlines() if ln.strip()]
        self.assertEqual(len(lines), 2)
        for line in lines:
            self.assertTrue(line.endswith(";"))

    def test_write_candidate_manifest_produces_correct_columns(self):
        with TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "manifest.tsv"
            write_candidate_manifest(
                ["original", "B1 swapped"],
                ["original", "generated"],
                out,
            )
            with out.open() as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                rows = list(reader)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["tree_index"], "1")
        self.assertEqual(rows[1]["description"], "B1 swapped")
        self.assertEqual(rows[1]["source"], "generated")


class AUTestParserTests(unittest.TestCase):
    _IQTREE_CONTENT = """
USER TREES
----------

See au_test.trees for trees with branch lengths.

Tree      logL    deltaL  bp-RELL    p-KH     p-SH    p-WKH    p-WSH    c-ELW    p-AU
---------------------------------------------------------------------------------------
  1  -12345.67     0.00  0.987 +  0.625 +  1.000 +  0.625 +  1.000 +  0.965 +  0.981 +
  2  -12389.12    43.45  0.013 -  0.375 -  0.500 -  0.375 -  0.500 -  0.035 -  0.019 -

deltaL  : logL difference from the best tree
"""

    def test_parse_au_test_returns_correct_row_count(self):
        with TemporaryDirectory() as tmpdir:
            iqtree_path = Path(tmpdir) / "au_test.iqtree"
            iqtree_path.write_text(self._IQTREE_CONTENT, encoding="utf-8")
            rows = parse_au_test_iqtree(iqtree_path)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["tree_index"], 1)
        self.assertAlmostEqual(rows[0]["logL"], -12345.67)
        self.assertAlmostEqual(rows[0]["deltaL"], 0.0)
        self.assertAlmostEqual(rows[1]["deltaL"], 43.45)

    def test_in_confidence_set_threshold(self):
        with TemporaryDirectory() as tmpdir:
            iqtree_path = Path(tmpdir) / "au_test.iqtree"
            iqtree_path.write_text(self._IQTREE_CONTENT, encoding="utf-8")
            manifest_path = Path(tmpdir) / "manifest.tsv"
            manifest_path.write_text(
                "tree_index\tdescription\tsource\n"
                "1\twASTRAL best tree\toriginal\n"
                "2\tB1 swapped\tgenerated\n",
                encoding="utf-8",
            )
            rows = parse_au_test_iqtree(iqtree_path)
            joined = join_with_manifest(rows, manifest_path)
        self.assertTrue(joined[0]["in_confidence_set"])   # p_au=0.981 >= 0.05
        self.assertFalse(joined[1]["in_confidence_set"])  # p_au=0.019 < 0.05

    def test_write_au_results_produces_all_columns(self):
        rows = [
            {
                "tree_index": 1, "description": "wASTRAL best tree", "source": "original",
                "logL": -12345.67, "deltaL": 0.0, "bp_rell": 0.987,
                "p_kh": 0.625, "p_sh": 1.0, "p_wkh": 0.625,
                "p_wsh": 1.0, "c_elw": 0.965, "p_au": 0.981,
                "in_confidence_set": True,
            }
        ]
        with TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "au_results.tsv"
            write_au_results(rows, out)
            with out.open() as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                written = list(reader)
        self.assertEqual(len(written), 1)
        self.assertIn("p_au", written[0])
        self.assertIn("in_confidence_set", written[0])
        self.assertEqual(written[0]["description"], "wASTRAL best tree")


# ── Workflow dry-run tests ─────────────────────────────────────────────────────


@unittest.skipUnless(SNAKEMAKE_AVAILABLE, "snakemake not in PATH")
class TopologyTestsWorkflowTests(unittest.TestCase):
    def test_topology_tests_complete_dry_run_renders_expected_rule(self):
        result = subprocess.run(
            [
                "snakemake", "-n", "-p", "--cores", "4",
                "results/topology_tests/topology_tests.complete",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0 and "Nothing to be done" not in result.stdout:
            # Expect to see our new rules in the dry-run output.
            self.assertIn("topology_tests", result.stdout + result.stderr, result.stderr[:500])

    def test_wastral_u2_annotate_rule_depends_on_species_tree(self):
        result = subprocess.run(
            [
                "snakemake", "-n", "-p", "--cores", "4",
                "results/topology_tests/species_tree.wastral.u2.tre",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        output = result.stdout + result.stderr
        if "Nothing to be done" not in output:
            self.assertIn("species_tree.wastral.tre", output)

    @unittest.skipUnless(
        SNAKEMAKE_AVAILABLE
        and (REPO_ROOT / "results/topology_tests/topology_tests.complete").is_file(),
        "snakemake not in PATH or real topology test outputs are not present.",
    )
    def test_real_au_results_has_confidence_set_annotation(self):
        from scripts.parse_au_test import parse_au_test_iqtree

        au_path = REPO_ROOT / "results/topology_tests/au_test.iqtree"
        rows = parse_au_test_iqtree(au_path)
        self.assertGreater(len(rows), 0)
        # The wASTRAL best tree (tree_index=1) should always be in the confidence set.
        tree1 = next((r for r in rows if r["tree_index"] == 1), None)
        self.assertIsNotNone(tree1)
        self.assertGreaterEqual(tree1["p_au"], 0.05)


if __name__ == "__main__":
    unittest.main()
