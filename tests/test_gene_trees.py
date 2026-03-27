import csv
import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.gene_trees import (
    aggregate_gene_tree_outputs,
    build_iqtree_command,
    gene_tree_output_paths,
    load_gene_tree_locus_ids,
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

    def test_aggregate_gene_tree_outputs_preserves_all_loci_once(self):
        with TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            retained_path = repo_root / "retained.tsv"
            manifest_path = repo_root / "gene_tree_manifest.tsv"
            aggregate_path = repo_root / "gene_trees.raw.tre"

            write_tsv(
                retained_path,
                [
                    {"locus_id": "2at4069", "decision": "retain"},
                    {"locus_id": "10at4069", "decision": "retain"},
                    {"locus_id": "11at4069", "decision": "exclude"},
                ],
                ["locus_id", "decision"],
            )

            for locus_id, model, tree_text in [
                ("2at4069", "LG+F", "(a:0.1,(b:0.1,c:0.1)/0.9:0.1);"),
                ("10at4069", "WAG+F", "(d:0.1,(e:0.1,f:0.1)/0.8:0.1);"),
            ]:
                paths = gene_tree_output_paths(locus_id)
                report_path = repo_root / paths["report"]
                tree_path = repo_root / paths["treefile"]
                report_path.parent.mkdir(parents=True, exist_ok=True)
                report_path.write_text(
                    f"Best-fit model according to BIC: {model}\nModel of substitution: {model}\n",
                    encoding="utf-8",
                )
                tree_path.write_text(tree_text + "\n", encoding="utf-8")

            aggregate_gene_tree_outputs(
                retained_path=retained_path,
                manifest_path=manifest_path,
                aggregate_path=aggregate_path,
                support_mode="abayes",
                repo_root=repo_root,
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


class GeneTreeWorkflowTests(unittest.TestCase):
    def test_gene_tree_rule_dry_run_renders_for_real_locus(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "-R",
                "infer_gene_tree",
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
        self.assertIn("rule infer_gene_tree:", result.stdout)
        self.assertIn("rule aggregate_gene_trees:", result.stdout)
        self.assertIn("results/gene_trees/per_locus/", result.stdout)

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
        self.assertEqual([row["locus_id"] for row in manifest_rows], sorted(retained_ids, key=lambda value: (int(value.split("at", 1)[0]), value.split("at", 1)[1])))


if __name__ == "__main__":
    unittest.main()
