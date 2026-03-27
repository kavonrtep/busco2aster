rule build_locus_matrix:
    input:
        summary=BUSCO_SUMMARY_TABLE,
        records=BUSCO_RECORDS_TABLE,
    output:
        LOCUS_TAXON_MATRIX,
    shell:
        (
            "python3 -m scripts.build_locus_matrix "
            "--summary {input.summary:q} "
            "--records {input.records:q} "
            "--output {output:q}"
        )


rule select_loci:
    input:
        matrix=LOCUS_TAXON_MATRIX,
    output:
        RETAINED_LOCI_TABLE,
    params:
        occupancy_threshold=float(config["occupancy_threshold"]),
    shell:
        (
            "python3 -m scripts.select_loci "
            "--matrix {input.matrix:q} "
            "--output {output:q} "
            "--occupancy-threshold {params.occupancy_threshold}"
        )
