rule prepare_visual_report_data:
    input:
        busco_summary=BUSCO_SUMMARY_TABLE,
        locus_matrix=LOCUS_TAXON_MATRIX,
        retained_loci=RETAINED_LOCI_TABLE,
        gene_tree_manifest=GENE_TREE_MANIFEST,
        species_tree=DEFAULT_SPECIES_TREE_OUTPUTS["treefile"],
        gcf_stat=GCF_OUTPUTS["stat"],
        gcf_branch=GCF_OUTPUTS["branch"],
        scfl_stat=SCFL_OUTPUTS["stat"],
        scfl_branch=SCFL_OUTPUTS["branch"],
        quartet_freqquad=WASTRAL_QUARTET_OUTPUTS["freqquad"],
        alignments_complete=ALIGNMENTS_COMPLETE,
    output:
        dataset_summary=REPORT_DATASET_SUMMARY,
        sample_qc=REPORT_SAMPLE_QC,
        locus_summary=REPORT_LOCUS_SUMMARY,
        alignment_summary=REPORT_ALIGNMENT_SUMMARY,
        branch_metrics=REPORT_BRANCH_METRICS,
        branch_alternatives=REPORT_BRANCH_ALTERNATIVES,
        gene_tree_heterogeneity=REPORT_GENE_TREE_HETEROGENEITY,
        topology_counts=REPORT_TOPOLOGY_COUNTS,
        species_tree_report=REPORT_SPECIES_TREE,
    params:
        repo_root=REPO_ROOT.as_posix(),
        alignment_dir=ALIGNMENT_DIR,
        output_dir=REPORT_DATA_DIR,
    shell:
        (
            "python3 -m scripts.prepare_visual_report_data "
            "--repo-root {params.repo_root:q} "
            "--busco-summary {input.busco_summary:q} "
            "--locus-matrix {input.locus_matrix:q} "
            "--retained-loci {input.retained_loci:q} "
            "--gene-tree-manifest {input.gene_tree_manifest:q} "
            "--species-tree {input.species_tree:q} "
            "--gcf-stat {input.gcf_stat:q} "
            "--gcf-branch {input.gcf_branch:q} "
            "--scfl-stat {input.scfl_stat:q} "
            "--scfl-branch {input.scfl_branch:q} "
            "--quartet-freqquad {input.quartet_freqquad:q} "
            "--alignment-dir {params.alignment_dir:q} "
            "--output-dir {params.output_dir:q}"
        )


rule render_visual_report:
    input:
        dataset_summary=REPORT_DATASET_SUMMARY,
        sample_qc=REPORT_SAMPLE_QC,
        locus_summary=REPORT_LOCUS_SUMMARY,
        alignment_summary=REPORT_ALIGNMENT_SUMMARY,
        branch_metrics=REPORT_BRANCH_METRICS,
        branch_alternatives=REPORT_BRANCH_ALTERNATIVES,
        gene_tree_heterogeneity=REPORT_GENE_TREE_HETEROGENEITY,
        topology_counts=REPORT_TOPOLOGY_COUNTS,
        species_tree_report=REPORT_SPECIES_TREE,
        qmd="reports/report.qmd",
    output:
        REPORT_HTML,
    params:
        data_dir=REPORT_DATA_DIR,
    conda:
        "../envs/report.yaml"
    shell:
        (
            "python3 -m scripts.render_visual_report "
            "--qmd {input.qmd:q} "
            "--data-dir {params.data_dir:q} "
            "--output {output:q}"
        )
