import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.cleanup_outputs import list_cleanup_targets, remove_cleanup_targets


REPO_ROOT = Path(__file__).resolve().parents[1]


class CleanupUnitTests(unittest.TestCase):
    def test_list_cleanup_targets_for_resume_gene_trees_includes_expected_paths(self):
        with TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            targets = [
                repo_root / "work/busco/sample_a/raw",
                repo_root / "results/busco/sample_a/raw",
                repo_root / "results/loci/raw_fastas",
                repo_root / "results/loci/logs",
                repo_root / "results/gene_trees/per_locus",
            ]
            for path in targets:
                path.mkdir(parents=True, exist_ok=True)

            planned = list_cleanup_targets(repo_root, "resume_gene_trees")

        self.assertEqual(
            {path.relative_to(repo_root).as_posix() for path in planned},
            {
                "work/busco/sample_a/raw",
                "results/busco/sample_a/raw",
                "results/loci/raw_fastas",
                "results/loci/logs",
                "results/gene_trees/per_locus",
            },
        )

    def test_remove_cleanup_targets_deletes_files_and_directories(self):
        with TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            file_path = repo_root / "results/concordance/scfl.log"
            dir_path = repo_root / "results/gene_trees/per_locus"
            file_path.parent.mkdir(parents=True, exist_ok=True)
            dir_path.mkdir(parents=True, exist_ok=True)
            file_path.write_text("log\n", encoding="utf-8")

            removed = remove_cleanup_targets([file_path, dir_path])

            self.assertEqual({path.name for path in removed}, {"scfl.log", "per_locus"})
            self.assertFalse(file_path.exists())
            self.assertFalse(dir_path.exists())


class CleanupWorkflowTests(unittest.TestCase):
    def test_cleanup_cli_dry_run_lists_targets_without_deleting(self):
        with TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            raw_root = repo_root / "results/busco/sample_a/raw"
            raw_root.mkdir(parents=True, exist_ok=True)

            result = subprocess.run(
                [
                    "python3",
                    "-m",
                    "scripts.cleanup_outputs",
                    "--repo-root",
                    str(repo_root),
                    "--mode",
                    "resume_alignments",
                    "--dry-run",
                ],
                cwd=REPO_ROOT,
                capture_output=True,
                text=True,
                check=True,
            )

            self.assertIn("results/busco/sample_a/raw", result.stdout)
            self.assertTrue(raw_root.exists())


if __name__ == "__main__":
    unittest.main()
