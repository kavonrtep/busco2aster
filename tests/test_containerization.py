import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import yaml

import run_pipeline


REPO_ROOT = Path(__file__).resolve().parents[1]


class ContainerizationUnitTests(unittest.TestCase):
    def test_collect_input_paths_resolves_manifest_and_assemblies_from_repo_root(self):
        with TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            config_dir = repo_root / "config"
            data_dir = repo_root / "data"
            config_dir.mkdir(parents=True)
            data_dir.mkdir(parents=True)

            assembly_path = data_dir / "sample1.fa"
            assembly_path.write_text(">chr1\nACGT\n", encoding="utf-8")

            samples_path = config_dir / "samples.tsv"
            samples_path.write_text(
                "sample_id\ttaxon_id\tassembly_fasta\n"
                "sample1\tTaxon one\tdata/sample1.fa\n",
                encoding="utf-8",
            )

            config_path = config_dir / "config.yaml"
            config_path.write_text("samples: config/samples.tsv\n", encoding="utf-8")
            config_obj = yaml.safe_load(config_path.read_text(encoding="utf-8"))

            resolved_samples, input_paths = run_pipeline.collect_input_paths(
                config_path=config_path,
                config_obj=config_obj,
                repo_root=repo_root,
            )

        self.assertEqual(resolved_samples, samples_path)
        self.assertEqual(input_paths["samples"], samples_path)
        self.assertEqual(input_paths["assembly:sample1"], assembly_path)

    def test_prepare_runtime_workdir_creates_fallback_symlinks(self):
        with TemporaryDirectory() as tmpdir:
            workdir = Path(tmpdir) / "run"
            pipeline_dir = Path(tmpdir) / "pipeline"
            (pipeline_dir / "config").mkdir(parents=True)
            (pipeline_dir / "scripts").mkdir(parents=True)

            run_pipeline.prepare_runtime_workdir(workdir, pipeline_dir)

            self.assertTrue((workdir / "config").is_symlink())
            self.assertTrue((workdir / "scripts").is_symlink())


class ContainerizationWorkflowTests(unittest.TestCase):
    def test_container_definition_links_cli_tools_into_usr_local_bin(self):
        definition_text = (REPO_ROOT / "busco2aster.def").read_text(encoding="utf-8")
        self.assertIn(
            "ln -sf /opt/tools/iqtree3/current/bin/iqtree3 /usr/local/bin/iqtree3",
            definition_text,
        )
        self.assertIn(
            "ln -sf /opt/tools/aster/current/bin/wastral /usr/local/bin/wastral",
            definition_text,
        )
        self.assertIn(
            "ln -sf /opt/tools/aster/current/bin/astral4 /usr/local/bin/astral4",
            definition_text,
        )

    def test_create_env_helper_workflow_dry_run_renders_expected_rules(self):
        result = subprocess.run(
            ["snakemake", "-s", "workflow/Snakefile_create_envs", "-n", "-p"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("rule create_env_busco:", result.stdout)
        self.assertIn("rule create_env_alignment:", result.stdout)
        self.assertIn("rule create_env_report:", result.stdout)

    def test_create_env_report_rule_validates_quarto_and_r_packages(self):
        helper_text = (REPO_ROOT / "workflow" / "rules" / "_create_envs.smk").read_text(
            encoding="utf-8"
        )
        self.assertIn("quarto --version", helper_text)
        self.assertIn("library(ggtree); library(phangorn)", helper_text)

    def test_github_workflow_smoke_tests_report_environment_inside_container(self):
        workflow_text = (REPO_ROOT / ".github" / "workflows" / "build-sif.yaml").read_text(
            encoding="utf-8"
        )
        self.assertIn("--snakefile /opt/pipeline/workflow/Snakefile_create_envs", workflow_text)
        self.assertIn("/tmp/busco2aster_env_report", workflow_text)

    def test_run_pipeline_wrapper_supports_local_dry_run(self):
        result = subprocess.run(
            [
                "python3",
                "run_pipeline.py",
                "--config",
                "config/config.yaml",
                "--repo-root",
                REPO_ROOT.as_posix(),
                "--directory",
                REPO_ROOT.as_posix(),
                "--target",
                "results/metadata/samples.validated.tsv",
                "--snakemake-args=--dry-run",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("Running:", result.stdout)
        self.assertTrue(
            "rule validate_manifest:" in result.stdout or "Nothing to be done" in result.stdout
        )


if __name__ == "__main__":
    unittest.main()
