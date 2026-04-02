rule extract_retained_dna_sample:
    input:
        matrix=LOCUS_TAXON_MATRIX,
        retained=RETAINED_LOCI_TABLE,
        assembly=ASSEMBLY_PREPARED_PLAIN_PATTERN,
        assembly_index=ASSEMBLY_PREPARED_PLAIN_INDEX_PATTERN,
    output:
        gff=DNA_AGGREGATED_GFF_PATTERN,
        fasta=DNA_EXTRACTED_FASTA_PATTERN,
        records=DNA_RECORDS_PATTERN,
        command=DNA_COMMAND_PATTERN,
        log=DNA_LOG_PATTERN,
        completion=DNA_COMPLETION_PATTERN,
    params:
        repo_root=REPO_ROOT.as_posix(),
        gffread_executable=str(config.get("gffread_executable", "gffread")),
    conda:
        "../envs/dna_extract.yaml"
    shell:
        (
            "python3 -m scripts.extract_retained_dna "
            "--matrix {input.matrix:q} "
            "--retained {input.retained:q} "
            "--sample-id {wildcards.sample:q} "
            "--genome-fasta {input.assembly:q} "
            "--repo-root {params.repo_root:q} "
            "--gffread-executable {params.gffread_executable:q}"
        )
