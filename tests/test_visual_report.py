import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.report_data import build_report_data_bundle
from scripts.tree_utils import canonical_topology_key, parse_newick, relabel_tree_with_branch_ids, render_newick


REPO_ROOT = Path(__file__).resolve().parents[1]


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


class TreeUtilsUnitTests(unittest.TestCase):
    def test_canonical_topology_key_is_child_order_invariant(self):
        tree_a = parse_newick("((a,b),(c,d));")
        tree_b = parse_newick("((d,c),(b,a));")
        self.assertEqual(canonical_topology_key(tree_a), canonical_topology_key(tree_b))

    def test_relabel_tree_with_branch_ids_assigns_stable_labels(self):
        relabeled = relabel_tree_with_branch_ids(parse_newick("(((a,b),c),d,e);"))
        rendered = render_newick(relabeled)
        self.assertIn("B1", rendered)
        self.assertIn("B2", rendered)


class VisualReportDataTests(unittest.TestCase):
    def test_report_data_bundle_writes_joined_branch_metrics(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            align_dir = tmp_path / "alignments"
            align_dir.mkdir()

            write_text(
                tmp_path / "busco_summary.tsv",
                "sample_id\ttaxon_id\tsanitized_taxon_id\tassembly_fasta\tbusco_version\tlineage_name\tcomplete_percent\tsingle_copy_percent\tmulti_copy_percent\tfragmented_percent\tmissing_percent\tinternal_stop_codon_count\n"
                "a\tA\ta\tdata/a.fa\t6.0.0\ttest\t99\t99\t0\t0\t1\t0\n"
                "b\tB\tb\tdata/b.fa\t6.0.0\ttest\t99\t99\t0\t0\t1\t0\n"
                "c\tC\tc\tdata/c.fa\t6.0.0\ttest\t99\t99\t0\t0\t1\t0\n"
                "d\tD\td\tdata/d.fa\t6.0.0\ttest\t99\t99\t0\t0\t1\t0\n",
            )
            write_text(
                tmp_path / "locus_matrix.tsv",
                "locus_id\tsample_id\tstatus\tinclude_in_occupancy\thas_internal_stop_codon\n"
                "l1\ta\tComplete\ttrue\tfalse\n"
                "l1\tb\tComplete\ttrue\tfalse\n"
                "l1\tc\tComplete\ttrue\tfalse\n"
                "l1\td\tComplete\ttrue\tfalse\n",
            )
            write_text(
                tmp_path / "retained_loci.tsv",
                "locus_id\tdecision\toccupancy\toccupancy_threshold\tduplicated_taxa\tstop_codon_taxa\tfailure_reasons\tqc_warnings\tlength_dispersion_observed\n"
                "l1\tretain\t1.0000\t0.8\t0\t0\t\t\tfalse\n",
            )
            write_text(
                tmp_path / "gene_tree_manifest.tsv",
                "locus_id\ttreefile\treport\tselected_model\tsupport_mode\tsupport_values_present\n"
                "l1\tgene_trees/l1.treefile\tgene_trees/l1.iqtree\tLG\tabayes\ttrue\n",
            )
            write_text(tmp_path / "species_tree.tre", "((a,b)1:0.1,(c,d):0.1);\n")
            write_text(
                tmp_path / "gcf.cf.stat",
                "# Concordance factor statistics\n"
                "ID\tgCF\tgCF_N\tgDF1\tgDF1_N\tgDF2\tgDF2_N\tgDFP\tgDFP_N\tgN\tLabel\tLength\n"
                "5\t80\t8\t10\t1\t5\t1\t5\t1\t10\t1\t0.1\n",
            )
            write_text(tmp_path / "gcf.cf.branch", "((a,b)5:0.1,(c,d):0.1);\n")
            write_text(
                tmp_path / "scfl.cf.stat",
                "# Concordance factor statistics\n"
                "ID\tsCF\tsCF_N\tsDF1\tsDF1_N\tsDF2\tsDF2_N\tsN\tLabel\tLength\n"
                "9\t70\t7\t20\t2\t10\t1\t10\t1\t0.1\n",
            )
            write_text(tmp_path / "scfl.cf.branch", "((a,b)9:0.1,(c,d):0.1);\n")
            write_text(
                tmp_path / "wastral_quartets.freqquad.tsv",
                "N1\tt1\t{a}|{b}#{c}|{d}\t1\t8\t10\n"
                "N1\tt2\t{a}|{c}#{b}|{d}\t0\t1\t10\n"
                "N1\tt3\t{a}|{d}#{b}|{c}\t0\t1\t10\n",
            )
            write_text(align_dir / "l1.aln.faa", ">a\nAAAA\n>b\nAAAA\n>c\nAAAA\n>d\nAAAA\n")
            gene_tree_dir = tmp_path / "gene_trees"
            gene_tree_dir.mkdir()
            write_text(gene_tree_dir / "l1.treefile", "((a,b),(c,d));\n")
            write_text(gene_tree_dir / "l1.iqtree", "Best-fit model according to BIC: LG\n")

            output_dir = tmp_path / "report_data"
            build_report_data_bundle(
                repo_root=tmp_path,
                busco_summary_path=tmp_path / "busco_summary.tsv",
                locus_matrix_path=tmp_path / "locus_matrix.tsv",
                retained_loci_path=tmp_path / "retained_loci.tsv",
                gene_tree_manifest_path=tmp_path / "gene_tree_manifest.tsv",
                species_tree_path=tmp_path / "species_tree.tre",
                gcf_stat_path=tmp_path / "gcf.cf.stat",
                gcf_branch_path=tmp_path / "gcf.cf.branch",
                scfl_stat_path=tmp_path / "scfl.cf.stat",
                scfl_branch_path=tmp_path / "scfl.cf.branch",
                quartet_freqquad_path=tmp_path / "wastral_quartets.freqquad.tsv",
                alignment_dir=align_dir,
                output_dir=output_dir,
            )

            branch_metrics = (output_dir / "branch_metrics.tsv").read_text(encoding="utf-8")
            report_tree = (output_dir / "species_tree.report.tre").read_text(encoding="utf-8")
            topology_counts = (output_dir / "topology_counts.tsv").read_text(encoding="utf-8")

        self.assertIn("B1", report_tree)
        self.assertIn("\t80.0000\t", branch_metrics)
        self.assertIn("\t70.0000\t", branch_metrics)
        self.assertIn("topology_id", topology_counts)


class VisualReportWorkflowTests(unittest.TestCase):
    def test_visual_report_rule_dry_run_renders_expected_targets(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "--cores",
                "4",
                "results/report/report.html",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        if "Nothing to be done" not in result.stdout:
            self.assertIn("rule prepare_visual_report_data:", result.stdout)
            self.assertIn("rule render_visual_report:", result.stdout)
            self.assertIn("results/report/data/branch_metrics.tsv", result.stdout)
        else:
            self.assertTrue((REPO_ROOT / "results/report/report.html").is_file())


if __name__ == "__main__":
    unittest.main()
