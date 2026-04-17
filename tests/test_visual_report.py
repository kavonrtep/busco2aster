import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.report_data import build_report_data_bundle
from scripts.tree_utils import (
    canonical_topology_key,
    fill_missing_branch_lengths,
    parse_newick,
    relabel_tree_with_branch_ids,
    render_newick,
)


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

    def test_fill_missing_branch_lengths_assigns_default_to_leaves(self):
        root = parse_newick("((a,b)B1:4.09,c)B2:0.3;")
        fill_missing_branch_lengths(root, default="0")
        rendered = render_newick(root)
        self.assertIn("a:0", rendered)
        self.assertIn("b:0", rendered)
        self.assertIn("c:0", rendered)
        self.assertIn("B1:4.09", rendered)
        self.assertIn("B2:0.3", rendered)

    def test_fill_missing_branch_lengths_leaves_root_unchanged(self):
        root = parse_newick("(a,b);")
        fill_missing_branch_lengths(root, default="0")
        self.assertIsNone(root.length)


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
                "locus_id\ttreefile\ttree_row_index\treport\tselected_model\tsupport_mode\tsupport_values_present\n"
                "l1\tgene_trees/all.raw.tre\t1\tgene_trees/l1.iqtree\tLG\tabayes\ttrue\n",
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
            write_text(gene_tree_dir / "all.raw.tre", "((a,b),(c,d));\n")
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

    def test_render_html_report_writes_self_contained_output(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            data_dir = tmp_path / "data"
            assets_dir = tmp_path / "assets"
            output_path = tmp_path / "results" / "report" / "report.html"
            template_path = tmp_path / "template.html.j2"

            data_dir.mkdir(parents=True)
            assets_dir.mkdir(parents=True)

            write_text(data_dir / "dataset_summary.tsv",
                "sequence_type\tsequence_length_unit\tsample_count\tcandidate_loci\tretained_loci\texcluded_loci\tmean_retained_missing_fraction\talignment_count\tmean_alignment_length_sites\tmedian_alignment_length_sites\tmin_alignment_length_sites\tmax_alignment_length_sites\tgene_tree_count\tcomplete_taxon_gene_tree_count\n"
                "protein\taa\t4\t10\t8\t2\t0.01\t8\t500.0\t450.0\t100\t1000\t8\t6\n")
            write_text(data_dir / "sample_qc.tsv",
                "sample_id\ttaxon_id\tsanitized_taxon_id\tcomplete_percent\tsingle_copy_percent\tmulti_copy_percent\tfragmented_percent\tmissing_percent\tinternal_stop_codon_count\tretained_loci_total\tretained_present_loci\tretained_missing_loci\tretained_missing_fraction\tretained_duplicated_loci\tretained_fragmented_loci\tretained_internal_stop_loci\n"
                "a\tA\ta\t99\t99\t0\t0\t1\t0\t8\t8\t0\t0.0\t0\t0\t0\n")
            write_text(data_dir / "locus_summary.tsv",
                "locus_id\tdecision\tfailure_reasons\n"
                "l1\tretain\t\n")
            write_text(data_dir / "alignment_summary.tsv",
                "locus_id\tsequence_count\talignment_length_sites\tgap_fraction\toccupancy\tlength_dispersion_observed\n"
                "l1\t4\t500\t0.01\t1.0\tfalse\n")
            write_text(data_dir / "branch_metrics.tsv",
                "label\tbranch_key\tsplit_left\tsplit_right\tspecies_tree_support\tspecies_tree_branch_length\tgcf_branch_id\tgcf\tgdf1\tgdf2\tgdfp\tgn\tscfl_branch_id\tscfl\tsdf1\tsdf2\tsn\taster_node_id\taster_local_posterior\taster_q1\taster_q2\taster_q3\taster_f1\taster_f2\taster_f3\taster_total_weight\n"
                "B1\ta,b||c,d\ta,b\tc,d\t0.99\t0.1\t5\t80\t10\t5\t5\t10\t9\t70\t20\t10\t100\tN1\t0.99\t0.8\t0.1\t0.1\t8\t1\t1\t10\n")
            write_text(data_dir / "branch_alternatives.tsv",
                "label\tbranch_key\tsource\talternative_id\tvalue\tsplit\n"
                "B1\ta,b||c,d\tgcf\tmain\t80\t\n")
            write_text(data_dir / "gene_tree_heterogeneity.tsv",
                "locus_id\tselected_model\tsupport_values_present\ttip_count\tis_complete_taxon\ttopology_key\tmatches_species_tree\trf_distance_to_species_tree\n"
                "l1\tLG\ttrue\t4\ttrue\ta,b||c,d\ttrue\t0\n")
            write_text(data_dir / "topology_counts.tsv",
                "topology_rank\ttopology_id\ttopology_key\texample_locus_id\tmatches_species_tree\tcount\tfraction\n"
                "1\tT1\ta,b||c,d\tl1\ttrue\t6\t1.0\n")
            write_text(data_dir / "species_tree.report.tre", "((a,b)B1:0.1,(c,d):0.1);\n")

            write_text(assets_dir / "plotly-basic.min.js", "/* plotly stub */")
            write_text(assets_dir / "underscore.min.js", "/* underscore stub */")
            write_text(assets_dir / "lodash.min.js", "/* lodash stub */")
            write_text(assets_dir / "phylotree.js", "/* phylotree stub */")

            write_text(template_path,
                '<!DOCTYPE html><html><body>'
                '<h1>{{ summary.sample_count }} taxa</h1>'
                '{% for b in branch_metrics %}<p>{{ b.label }}: gcf={{ b.gcf }}</p>{% endfor %}'
                '<script>{{ plotly_js }}</script>'
                '<script>{{ phylotree_js }}</script>'
                '{% if topo_tests_available %}<p>topo tests</p>{% endif %}'
                '</body></html>')

            subprocess.run(
                [
                    "python3", "-m", "scripts.render_html_report",
                    "--template", template_path.as_posix(),
                    "--data-dir", data_dir.as_posix(),
                    "--assets-dir", assets_dir.as_posix(),
                    "--output", output_path.as_posix(),
                ],
                cwd=REPO_ROOT,
                check=True,
            )

            self.assertTrue(output_path.is_file())
            html = output_path.read_text(encoding="utf-8")
            self.assertIn("4 taxa", html)
            self.assertIn("B1: gcf=80", html)
            self.assertIn("plotly stub", html)
            self.assertNotIn("topo tests", html)

    def test_render_html_report_with_real_template_and_test_data(self):
        """Render the actual template against the yeast test data if available."""
        test_data_dir = Path("/mnt/ceph/tmp/busco2aster_test/results/report/data")
        if not test_data_dir.is_dir():
            self.skipTest("Test data not available at /mnt/ceph/tmp/busco2aster_test")

        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "report.html"
            subprocess.run(
                [
                    "python3", "-m", "scripts.render_html_report",
                    "--template", (REPO_ROOT / "reports/template.html.j2").as_posix(),
                    "--data-dir", test_data_dir.as_posix(),
                    "--assets-dir", (REPO_ROOT / "reports/assets").as_posix(),
                    "--output", output_path.as_posix(),
                ],
                cwd=REPO_ROOT,
                check=True,
            )
            self.assertTrue(output_path.is_file())
            html = output_path.read_text(encoding="utf-8")
            self.assertIn("saccharomyces", html)
            self.assertIn("Plotly.newPlot", html)
            self.assertIn("phylotree", html)
            self.assertIn("Topological Uncertainty", html)
            self.assertGreater(len(html), 100_000)


class VisualReportWorkflowTests(unittest.TestCase):
    def test_visual_report_rule_dry_run_renders_expected_targets(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "--rerun-incomplete",
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
            self.assertTrue(
                "rule prepare_visual_report_data:" in result.stdout
                or "Updating job prepare_visual_report_data." in result.stdout
            )
            self.assertIn("rule render_visual_report:", result.stdout)
            self.assertTrue(
                "results/report/data/branch_metrics.tsv" in result.stdout
                or "results/report/data/dataset_summary.tsv" in result.stdout
            )
        else:
            self.assertTrue((REPO_ROOT / "results/report/report.html").is_file())


if __name__ == "__main__":
    unittest.main()
