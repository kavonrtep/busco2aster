import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import yaml

import run_pipeline


REPO_ROOT = Path(__file__).resolve().parents[1]


class ContainerizationUnitTests(unittest.TestCase):
    def test_compute_thread_overrides_uses_auto_policy(self):
        overrides = run_pipeline.compute_thread_overrides(
            total_cores=40,
            sample_count=25,
            config_obj={},
            snakemake_args="",
        )

        self.assertEqual(overrides["prepare_assembly"], 10)
        self.assertEqual(overrides["run_busco"], 10)
        self.assertEqual(overrides["align_locus_batch"], 1)
        self.assertEqual(overrides["infer_gene_trees"], 40)
        self.assertEqual(overrides["infer_gene_concordance"], 40)
        self.assertEqual(overrides["infer_site_concordance"], 40)
        self.assertEqual(overrides["annotate_species_tree_quartets"], 40)
        self.assertEqual(overrides["infer_species_tree_wastral"], 40)

    def test_compute_thread_overrides_respects_fixed_policy(self):
        overrides = run_pipeline.compute_thread_overrides(
            total_cores=40,
            sample_count=25,
            config_obj={"thread_policy": "fixed"},
            snakemake_args="",
        )
        self.assertEqual(overrides, {})

    def test_compute_thread_overrides_skips_auto_when_set_threads_is_explicit(self):
        overrides = run_pipeline.compute_thread_overrides(
            total_cores=40,
            sample_count=25,
            config_obj={},
            snakemake_args="--set-threads infer_gene_trees=12",
        )
        self.assertEqual(overrides, {})

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
        self.assertEqual(input_paths["samples"].path, samples_path)
        self.assertEqual(input_paths["assembly:sample1"].path, assembly_path)

    def test_collect_input_paths_strips_trailing_whitespace_from_manifest_paths(self):
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
                "sample1\tTaxon one\tdata/sample1.fa \n",
                encoding="utf-8",
            )

            config_path = config_dir / "config.yaml"
            config_path.write_text("samples: config/samples.tsv\n", encoding="utf-8")
            config_obj = yaml.safe_load(config_path.read_text(encoding="utf-8"))

            _, input_paths = run_pipeline.collect_input_paths(
                config_path=config_path,
                config_obj=config_obj,
                repo_root=repo_root,
            )

        self.assertEqual(input_paths["assembly:sample1"].path, assembly_path)
        self.assertEqual(input_paths["assembly:sample1"].raw_text, "data/sample1.fa ")

    def test_missing_path_diagnosis_reports_missing_mount_root(self):
        diagnosis = run_pipeline.missing_path_diagnosis(
            Path("/definitely_missing_mount/example/sample.fa")
        )
        self.assertIn("mount root", diagnosis)

    def test_missing_path_diagnosis_reports_missing_file_inside_visible_tree(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            visible_dir = root / "visible"
            visible_dir.mkdir()

            diagnosis = run_pipeline.missing_path_diagnosis(visible_dir / "missing.fa")

        self.assertIn("visible tree", diagnosis)

    def test_missing_path_diagnosis_reports_whitespace_manifest_value(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            diagnosis = run_pipeline.missing_path_diagnosis(
                root / "sample.fa",
                "sample.fa ",
            )

        self.assertIn("whitespace", diagnosis)

    def test_missing_path_diagnosis_reports_missing_symlink_target_mount(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            link_path = root / "sample.fa"
            link_path.symlink_to("/mnt/not_mounted/data/sample.fa")

            diagnosis = run_pipeline.missing_path_diagnosis(link_path)

        self.assertIn("symlink target", diagnosis)
        self.assertIn("/mnt/not_mounted", diagnosis)

    def test_suggested_bind_dir_uses_mount_root_for_non_visible_absolute_tree(self):
        bind_dir = run_pipeline.suggested_bind_dir(Path("/mnt/not_mounted/example/sample.fa"))
        self.assertEqual(bind_dir, Path("/mnt/not_mounted"))

    def test_suggested_bind_dir_uses_existing_ancestor_for_wrong_location(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            visible_dir = root / "visible"
            visible_dir.mkdir()

            bind_dir = run_pipeline.suggested_bind_dir(visible_dir / "missing" / "sample.fa")

        self.assertEqual(bind_dir, visible_dir)

    def test_suggested_bind_dir_uses_symlink_target_mount_root(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            link_path = root / "sample.fa"
            link_path.symlink_to("/mnt/not_mounted/data/sample.fa")

            bind_dir = run_pipeline.suggested_bind_dir(link_path)

        self.assertEqual(bind_dir, Path("/mnt/not_mounted"))

    def test_prepare_runtime_workdir_creates_fallback_symlinks(self):
        with TemporaryDirectory() as tmpdir:
            workdir = Path(tmpdir) / "run"
            pipeline_dir = Path(tmpdir) / "pipeline"
            (pipeline_dir / "config").mkdir(parents=True)
            (pipeline_dir / "scripts").mkdir(parents=True)
            (pipeline_dir / "reports").mkdir(parents=True)

            run_pipeline.prepare_runtime_workdir(workdir, pipeline_dir)

            self.assertTrue((workdir / "config").is_symlink())
            self.assertTrue((workdir / "scripts").is_symlink())
            self.assertTrue((workdir / "reports").is_symlink())

    def test_build_snakemake_command_defaults_to_all_target(self):
        command = run_pipeline.build_snakemake_command(
            config_path=Path("/tmp/config.yaml"),
            repo_root=Path("/tmp/repo"),
            workdir=Path("/tmp/work"),
            threads=4,
            target=None,
            snakemake_args="--dry-run",
            thread_overrides={"infer_gene_trees": 4},
        )

        self.assertIn("all", command)
        self.assertNotIn("results/report/report.html", command)
        self.assertIn("--set-threads", command)
        self.assertIn("infer_gene_trees=4", command)


class ContainerizationWorkflowTests(unittest.TestCase):
    def test_minimal_config_template_omits_explicit_executable_paths(self):
        config_obj = yaml.safe_load(
            (REPO_ROOT / "config" / "config.template.yaml").read_text(encoding="utf-8")
        )

        self.assertEqual(config_obj["samples"], "config/samples.tsv")
        self.assertEqual(config_obj["busco_lineage"], "solanales_odb12")
        self.assertNotIn("iqtree_executable", config_obj)
        self.assertNotIn("wastral_executable", config_obj)
        self.assertNotIn("astral4_executable", config_obj)

    def test_workflow_dry_run_accepts_template_without_executable_paths(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "-R",
                "infer_species_tree_wastral",
                "--configfile",
                (REPO_ROOT / "config" / "config.template.yaml").as_posix(),
                "--cores",
                "4",
                "results/species_tree/species_tree.wastral.complete",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("rule infer_species_tree_wastral:", result.stdout)

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
        self.assertIn("rule create_env_assembly_prep:", result.stdout)
        self.assertIn("rule create_env_dna_extract:", result.stdout)
        self.assertIn("rule create_env_alignment:", result.stdout)
        self.assertNotIn("create_env_report", result.stdout)

    def test_create_env_assembly_prep_rule_validates_seqkit(self):
        helper_text = (REPO_ROOT / "workflow" / "rules" / "_create_envs.smk").read_text(
            encoding="utf-8"
        )
        self.assertIn("/tmp/busco2aster_env_assembly_prep", helper_text)
        self.assertIn("seqkit version", helper_text)

    def test_create_env_dna_extract_rule_validates_gffread(self):
        helper_text = (REPO_ROOT / "workflow" / "rules" / "_create_envs.smk").read_text(
            encoding="utf-8"
        )
        self.assertIn("/tmp/busco2aster_env_dna_extract", helper_text)
        self.assertIn("gffread --version", helper_text)

    def test_create_envs_does_not_include_report_env(self):
        helper_text = (REPO_ROOT / "workflow" / "rules" / "_create_envs.smk").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("report.yaml", helper_text)
        self.assertNotIn("quarto", helper_text)

    def test_github_workflow_smoke_tests_jinja2_inside_container(self):
        workflow_text = (REPO_ROOT / ".github" / "workflows" / "build-sif.yaml").read_text(
            encoding="utf-8"
        )
        self.assertIn("--snakefile /opt/pipeline/workflow/Snakefile_create_envs", workflow_text)
        self.assertIn("/tmp/busco2aster_env_assembly_prep", workflow_text)
        self.assertIn("/tmp/busco2aster_env_dna_extract", workflow_text)
        self.assertIn("import jinja2", workflow_text)

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
