import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from scripts.busco import (
    BuscoValidationError,
    busco_output_paths,
    hydrate_stable_busco_sequence_paths,
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
            paths["single_copy_faa"],
            "results/busco/solanum_chilense/busco_sequences/single_copy.faa",
        )
        self.assertEqual(
            paths["single_copy_gff"],
            "results/busco/solanum_chilense/busco_sequences/single_copy.gff",
        )
        self.assertEqual(
            paths["single_copy_per_locus_dir"],
            "results/busco/solanum_chilense/busco_sequences/single_copy_busco_sequences",
        )
        self.assertEqual(paths["completion"], "results/busco/solanum_chilense/run.complete")
        self.assertEqual(paths["short_summary"], "results/busco/solanum_chilense/short_summary.txt")

    def test_hydrate_stable_busco_sequence_paths_prefers_concatenated_files(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            repo_root = Path(temp_dir)
            stable_paths = busco_output_paths("sample1")
            stable_faa = repo_root / stable_paths["single_copy_faa"]
            stable_faa.parent.mkdir(parents=True)
            stable_faa.write_text(">busco_a\nMA\n", encoding="utf-8")
            hydrated = hydrate_stable_busco_sequence_paths(
                {
                    "single_copy_busco_sequences": "results/busco/sample1/raw/run/busco_sequences/single_copy_busco_sequences",
                    "multi_copy_busco_sequences": "results/busco/sample1/raw/run/busco_sequences/multi_copy_busco_sequences",
                },
                stable_paths,
                repo_root,
            )
            self.assertEqual(
                hydrated["single_copy_busco_sequences"], stable_paths["single_copy_faa"]
            )
            self.assertEqual(hydrated["single_copy_busco_gff"], stable_paths["single_copy_gff"])
            self.assertEqual(
                hydrated["multi_copy_busco_sequences"],
                "results/busco/sample1/raw/run/busco_sequences/multi_copy_busco_sequences",
            )

    def test_hydrate_stable_busco_sequence_paths_keeps_per_locus_directories(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            repo_root = Path(temp_dir)
            stable_paths = busco_output_paths("sample1")
            (repo_root / stable_paths["single_copy_per_locus_dir"]).mkdir(parents=True)
            hydrated = hydrate_stable_busco_sequence_paths(
                {
                    "single_copy_busco_sequences": "results/busco/sample1/raw/run/busco_sequences/single_copy_busco_sequences",
                },
                stable_paths,
                repo_root,
            )
            self.assertEqual(
                hydrated["single_copy_busco_sequences"],
                stable_paths["single_copy_per_locus_dir"],
            )
            self.assertNotIn("single_copy_busco_gff", hydrated)


def _rule_block(text: str, rule_name: str) -> str:
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith(f"rule {rule_name}:"))
    end = start + 1
    while end < len(lines) and (not lines[end].strip() or lines[end].startswith((" ", "\t"))):
        end += 1
    return "\n".join(lines[start:end])


def _directive_block(rule_text: str, directive: str) -> str:
    lines = rule_text.splitlines()
    start = next(
        (i for i, line in enumerate(lines) if line.strip() == f"{directive}:"),
        None,
    )
    if start is None:
        return ""
    end = start + 1
    while end < len(lines) and lines[end].startswith("        "):
        end += 1
    return "\n".join(lines[start + 1 : end])


class BuscoScratchRootTests(unittest.TestCase):
    """The BUSCO scratch root must stay out of `params` (GitHub issue #1).

    Under a batch scheduler the scratch root is job-scoped (PBS `$SCRATCHDIR`),
    so a `params` entry derived from it changes on every submission, fires
    snakemake's params rerun-trigger, and re-runs BUSCO for every sample even
    though nothing about the data or the method changed.
    """

    def setUp(self):
        self.rule_text = _rule_block(
            (REPO_ROOT / "workflow" / "rules" / "busco.smk").read_text(encoding="utf-8"),
            "run_busco",
        )

    def test_scratch_root_is_not_a_param(self):
        params = _directive_block(self.rule_text, "params")
        self.assertNotIn("BUSCO_WORK_ROOT", params)
        self.assertNotIn("raw_root", params)
        self.assertNotIn("{params.raw_root", self.rule_text)

    def test_scratch_root_is_declared_as_a_resource(self):
        resources = _directive_block(self.rule_text, "resources")
        self.assertIn("raw_root=", resources)
        self.assertIn("BUSCO_WORK_ROOT", resources)

    def test_recorded_command_still_uses_the_real_scratch_out_path(self):
        self.assertIn("--out_path {resources.raw_root:q}", self.rule_text)
        self.assertIn("--raw-root {resources.raw_root:q}", self.rule_text)


class BuscoOfflineLineageTests(unittest.TestCase):
    """`run_busco` must not touch the network (GitHub issue #3).

    Concurrent jobs share one `--download_path`, and an online BUSCO rewrites
    `file_versions.tsv` there in place; a job that reads it mid-rewrite aborts
    with "<lineage> is not a valid option for 'lineages'" even though the
    dataset is fully downloaded.
    """

    def setUp(self):
        self.rule_text = (REPO_ROOT / "workflow" / "rules" / "busco.smk").read_text(
            encoding="utf-8"
        )

    def test_recorded_busco_command_runs_offline(self):
        run_busco = _rule_block(self.rule_text, "run_busco")
        self.assertIn("--offline", run_busco)
        self.assertIn("--download_path {params.download_path:q} --offline", run_busco)

    def test_lineage_download_rules_stay_online(self):
        """The dataset has to be fetched by something; these two rules never run concurrently."""
        for rule_name in ("capture_busco_metadata", "fetch_busco_lineage"):
            self.assertNotIn("--offline", _rule_block(self.rule_text, rule_name))

    def test_run_busco_waits_for_the_downloaded_dataset(self):
        """Offline mode is only safe because the dataset is already an input."""
        run_busco = _rule_block(self.rule_text, "run_busco")
        self.assertIn("lineage_dataset=ancient(BUSCO_LINEAGE_DATASET_CFG)", run_busco)


@unittest.skipUnless(shutil.which("snakemake"), "snakemake is not installed")
class SnakemakeRerunTriggerContractTests(unittest.TestCase):
    """Guard the snakemake behaviour the scratch-root fix relies on.

    A `params` value that changes between runs re-triggers the job; the same
    value expressed as a `resources` entry does not. If a snakemake upgrade
    changes that, `run_busco` silently goes back to re-running every submission.
    """

    SNAKEFILE = """
configfile: "cfg.yaml"
WORK_ROOT = config["work_root"]

rule all:
    input:
        "from_params.txt",
        "from_resources.txt",

rule from_params:
    output:
        "from_params.txt",
    params:
        raw_root=lambda wildcards: f"{WORK_ROOT}/raw",
    shell:
        "echo {params.raw_root} > {output}"

rule from_resources:
    output:
        "from_resources.txt",
    resources:
        raw_root=lambda wildcards: f"{WORK_ROOT}/raw",
    shell:
        "echo {resources.raw_root} > {output}"
"""

    def _snakemake(self, workdir: Path, *args: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            ["snakemake", "--cores", "1", *args],
            cwd=workdir,
            capture_output=True,
            text=True,
            check=True,
        )

    def test_resources_survive_a_changed_scratch_root_but_params_do_not(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            workdir = Path(temp_dir)
            (workdir / "Snakefile").write_text(self.SNAKEFILE, encoding="utf-8")

            (workdir / "cfg.yaml").write_text("work_root: /scratch/job_A\n", encoding="utf-8")
            self._snakemake(workdir)
            self.assertEqual(
                (workdir / "from_resources.txt").read_text(encoding="utf-8").strip(),
                "/scratch/job_A/raw",
            )

            # Same workflow, same data, resubmitted under a new job id.
            (workdir / "cfg.yaml").write_text("work_root: /scratch/job_B\n", encoding="utf-8")
            planned = self._snakemake(workdir, "--dry-run").stdout

            self.assertIn("from_params", planned)
            self.assertNotIn("rule from_resources:", planned)


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
        self.assertTrue(
            "rule prepare_assembly:" in result.stdout
            or "work/assemblies_prepared/solanum_chilense.fa.gz" in result.stdout
        )
        self.assertIn("--cpu 4", result.stdout)
        self.assertIn("--lineage_dataset solanales_odb12", result.stdout)
        self.assertIn("--offline", result.stdout)
        self.assertIn("work/assemblies_prepared/solanum_chilense.fa.gz", result.stdout)
        self.assertIn("results/busco/solanum_chilense/run.complete", result.stdout)
        self.assertIn("results/busco/solanum_chilense/busco_sequences", result.stdout)


if __name__ == "__main__":
    unittest.main()
