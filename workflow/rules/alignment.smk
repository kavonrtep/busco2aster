rule export_locus_fasta:
    input:
        matrix=LOCUS_TAXON_MATRIX,
        retained=RETAINED_LOCI_TABLE,
    output:
        "results/loci/raw_fastas/{locus_id}.faa",
    params:
        repo_root=str(REPO_ROOT),
    shell:
        (
            "python3 -m scripts.export_locus_fasta "
            "--locus-id {wildcards.locus_id:q} "
            "--matrix {input.matrix:q} "
            "--retained {input.retained:q} "
            "--output {output:q} "
            "--repo-root {params.repo_root:q}"
        )


rule align_locus:
    input:
        "results/loci/raw_fastas/{locus_id}.faa",
    output:
        alignment="results/loci/alignments/{locus_id}.aln.faa",
    log:
        "results/loci/logs/mafft/{locus_id}.log",
    threads:
        get_thread_count("alignment")
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {output.alignment})"
        mkdir -p "$(dirname {log})"
        mafft --amino --anysymbol --auto --thread {threads} {input:q} > {output.alignment:q} 2> {log:q}
        """
