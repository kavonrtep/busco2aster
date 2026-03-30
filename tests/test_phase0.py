import importlib
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class Phase0RepositoryTests(unittest.TestCase):
    def test_expected_directories_exist(self):
        expected = [
            REPO_ROOT / "config",
            REPO_ROOT / "workflow",
            REPO_ROOT / "workflow" / "rules",
            REPO_ROOT / "workflow" / "envs",
            REPO_ROOT / "scripts",
            REPO_ROOT / "tests",
        ]
        for path in expected:
            with self.subTest(path=path):
                self.assertTrue(path.is_dir(), f"Missing directory: {path}")

    def test_snakefile_exists(self):
        self.assertTrue((REPO_ROOT / "Snakefile").is_file())

    def test_scripts_package_imports(self):
        module = importlib.import_module("scripts")
        self.assertIsNotNone(module)

    def test_snakefile_exports_repo_root_on_pythonpath_for_shell_rules(self):
        snakefile_text = (REPO_ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertIn("shell.prefix(", snakefile_text)
        self.assertIn("export PYTHONPATH=", snakefile_text)
        self.assertIn("${{{{PYTHONPATH-}}}}", snakefile_text)


if __name__ == "__main__":
    unittest.main()
