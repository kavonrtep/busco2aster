import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.report import REPORT_SECTION_TITLES, render_report


REPO_ROOT = Path(__file__).resolve().parents[1]


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


class ReportUnitTests(unittest.TestCase):
    def test_render_report_includes_required_sections_and_artifacts(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            busco_path = tmp_path / "busco_summary.tsv"
            retained_path = tmp_path / "retained_loci.tsv"
            gene_manifest_path = tmp_path / "gene_tree_manifest.tsv"
            species_tree_path = tmp_path / "species_tree.tre"
            species_log_path = tmp_path / "species_tree.log"

            write_text(
                busco_path,
                "\t".join(
                    [
                        "sample_id",
                        "busco_version",
                        "lineage_name",
                        "complete_percent",
                        "single_copy_percent",
                        "multi_copy_percent",
                        "fragmented_percent",
                        "missing_percent",
                        "internal_stop_codon_count",
                    ]
                )
                + "\n"
                + "\t".join(["sample_a", "6.0.0", "solanales_odb12", "99.6", "99.1", "0.5", "0.0", "0.4", "3"])
                + "\n",
            )
            write_text(
                retained_path,
                "\t".join(
                    [
                        "locus_id",
                        "decision",
                        "occupancy",
                        "occupancy_threshold",
                        "duplicated_taxa",
                        "stop_codon_taxa",
                        "failure_reasons",
                        "qc_warnings",
                    ]
                )
                + "\n"
                + "\t".join(["1at1", "retain", "1.0000", "0.8", "0", "0", "", "length_dispersion_observed"])
                + "\n"
                + "\t".join(["2at1", "exclude", "0.5000", "0.8", "1", "1", "occupancy_below_threshold", "duplicated_taxa_present,internal_stop_codon_present"])
                + "\n",
            )
            write_text(
                gene_manifest_path,
                "\t".join(
                    [
                        "locus_id",
                        "selected_model",
                        "support_mode",
                        "support_values_present",
                    ]
                )
                + "\n"
                + "\t".join(["1at1", "LG+F", "abayes", "true"])
                + "\n"
                + "\t".join(["2at1", "WAG", "abayes", "false"])
                + "\n",
            )
            write_text(species_tree_path, "((a,b),c);\n")
            write_text(
                species_log_path,
                "\n".join(
                    [
                        "Version: v1.24.3.8",
                        "#Genetrees: 2",
                        "#Species: 3",
                        "#Rounds: 4",
                        "#Samples: 4",
                        "#Threads: 4",
                        "Score: 12.5",
                        "Final Tree: ((a,b),c);",
                    ]
                )
                + "\n",
            )

            report_text = render_report(
                busco_summary_path=busco_path,
                retained_loci_path=retained_path,
                gene_tree_manifest_path=gene_manifest_path,
                species_tree_path=species_tree_path,
                species_tree_log_path=species_log_path,
                species_tree_backend="wastral",
            )

        for title in REPORT_SECTION_TITLES:
            self.assertIn(f"## {title}", report_text)
        self.assertIn("results/qc/busco_summary.tsv", report_text)
        self.assertIn("results/qc/retained_loci.tsv", report_text)
        self.assertIn("gene_tree_manifest.tsv", report_text)
        self.assertIn("species_tree.tre", report_text)

    def test_render_report_handles_missing_optional_concordance_paths(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            busco_path = tmp_path / "busco_summary.tsv"
            retained_path = tmp_path / "retained_loci.tsv"
            gene_manifest_path = tmp_path / "gene_tree_manifest.tsv"
            species_tree_path = tmp_path / "species_tree.tre"
            species_log_path = tmp_path / "species_tree.log"

            write_text(
                busco_path,
                "sample_id\tbusco_version\tlineage_name\tcomplete_percent\tsingle_copy_percent\tmulti_copy_percent\tfragmented_percent\tmissing_percent\tinternal_stop_codon_count\n"
                "sample_a\t6.0.0\tsolanales_odb12\t99.6\t99.1\t0.5\t0.0\t0.4\t3\n",
            )
            write_text(
                retained_path,
                "locus_id\tdecision\toccupancy\toccupancy_threshold\tduplicated_taxa\tstop_codon_taxa\tfailure_reasons\tqc_warnings\n"
                "1at1\tretain\t1.0000\t0.8\t0\t0\t\t\n",
            )
            write_text(
                gene_manifest_path,
                "locus_id\tselected_model\tsupport_mode\tsupport_values_present\n"
                "1at1\tLG+F\tabayes\ttrue\n",
            )
            write_text(species_tree_path, "((a,b),c);\n")
            write_text(species_log_path, "Version: v1.24.3.8\n#Genetrees: 1\n#Species: 3\n#Rounds: 4\n#Samples: 4\n#Threads: 4\nScore: 12.5\n")

            report_text = render_report(
                busco_summary_path=busco_path,
                retained_loci_path=retained_path,
                gene_tree_manifest_path=gene_manifest_path,
                species_tree_path=species_tree_path,
                species_tree_log_path=species_log_path,
                species_tree_backend="wastral",
                concordance_paths=[tmp_path / "missing.tsv"],
            )

        self.assertIn("Not implemented in v1", report_text)


class ReportWorkflowTests(unittest.TestCase):
    def test_report_rule_dry_run_uses_expected_upstream_artifacts(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "--cores",
                "4",
                "results/report/report.md",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        if "Nothing to be done" not in result.stdout:
            self.assertIn("rule render_report:", result.stdout)
            self.assertIn("results/qc/busco_summary.tsv", result.stdout)
            self.assertIn("results/qc/retained_loci.tsv", result.stdout)
            self.assertIn("results/gene_trees/gene_tree_manifest.tsv", result.stdout)
            self.assertIn("results/species_tree/species_tree.wastral.tre", result.stdout)
        else:
            self.assertIn("Updating job render_report.", result.stdout)

    @unittest.skipUnless((REPO_ROOT / "results/report/report.md").is_file(), "Real report output is not present.")
    def test_real_report_references_species_tree_output(self):
        report_text = (REPO_ROOT / "results/report/report.md").read_text(encoding="utf-8")
        self.assertIn("results/species_tree/species_tree.wastral.tre", report_text)
        self.assertIn("## Species Tree", report_text)


if __name__ == "__main__":
    unittest.main()
