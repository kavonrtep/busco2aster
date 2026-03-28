rule export_retained_fastas:
    input:
        matrix=LOCUS_TAXON_MATRIX,
        retained=RETAINED_LOCI_TABLE,
    output:
        raw_fastas=directory(RAW_FASTA_DIR),
        manifest=RAW_FASTA_MANIFEST,
        completion="results/loci/raw_fastas.complete",
    params:
        repo_root=str(REPO_ROOT),
    shell:
        (
            "python3 -m scripts.export_retained_fastas "
            "--matrix {input.matrix:q} "
            "--retained {input.retained:q} "
            "--output-dir {output.raw_fastas:q} "
            "--manifest {output.manifest:q} "
            "--repo-root {params.repo_root:q} "
            "&& touch {output.completion:q}"
        )


rule align_locus:
    input:
        export_complete=ancient(rules.export_retained_fastas.output.completion),
    output:
        alignment="results/loci/alignments/{locus_id}.aln.faa",
    log:
        "results/loci/logs/mafft/{locus_id}.log",
    params:
        raw_fasta=lambda wildcards: f"{RAW_FASTA_DIR}/{wildcards.locus_id}.faa",
    threads:
        get_thread_count("alignment")
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {output.alignment})"
        mkdir -p "$(dirname {log})"
        mafft --amino --anysymbol --auto --thread {threads} {params.raw_fasta:q} > {output.alignment:q} 2> {log:q}
        """
