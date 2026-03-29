rule render_report:
    input:
        busco_summary=BUSCO_SUMMARY_TABLE,
        retained_loci=RETAINED_LOCI_TABLE,
        gene_tree_manifest=GENE_TREE_MANIFEST,
        species_tree=DEFAULT_SPECIES_TREE_OUTPUTS["treefile"],
        species_tree_log=DEFAULT_SPECIES_TREE_OUTPUTS["log"],
        species_tree_complete=DEFAULT_SPECIES_TREE_OUTPUTS["completion"],
        gcf_stat=GCF_OUTPUTS["stat"],
        scfl_stat=SCFL_OUTPUTS["stat"],
        quartet_tree=WASTRAL_QUARTET_OUTPUTS["tree"],
        quartet_freqquad=WASTRAL_QUARTET_OUTPUTS["freqquad"],
    output:
        REPORT_MARKDOWN,
    params:
        backend=SPECIES_TREE_BACKEND,
    shell:
        (
            "python3 -m scripts.render_report "
            "--busco-summary {input.busco_summary:q} "
            "--retained-loci {input.retained_loci:q} "
            "--gene-tree-manifest {input.gene_tree_manifest:q} "
            "--species-tree {input.species_tree:q} "
            "--species-tree-log {input.species_tree_log:q} "
            "--backend {params.backend:q} "
            "--concordance-path {input.gcf_stat:q} "
            "--concordance-path {input.scfl_stat:q} "
            "--concordance-path {input.quartet_tree:q} "
            "--concordance-path {input.quartet_freqquad:q} "
            "--output {output:q}"
        )
