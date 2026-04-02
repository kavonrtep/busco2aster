rule export_retained_fastas:
    input:
        matrix=LOCUS_TAXON_MATRIX,
        retained=RETAINED_LOCI_TABLE,
        dna_sequences=lambda wildcards: DNA_COMPLETION_TARGETS if SEQUENCE_TYPE == "dna" else [],
    output:
        raw_fastas=directory(RAW_FASTA_DIR),
        manifest=RAW_FASTA_MANIFEST,
        completion="results/loci/raw_fastas.complete",
    params:
        repo_root=str(REPO_ROOT),
        sequence_type=SEQUENCE_TYPE,
    shell:
        (
            "python3 -m scripts.export_retained_fastas "
            "--matrix {input.matrix:q} "
            "--retained {input.retained:q} "
            "--output-dir {output.raw_fastas:q} "
            "--manifest {output.manifest:q} "
            "--repo-root {params.repo_root:q} "
            "--sequence-type {params.sequence_type:q} "
            "&& touch {output.completion:q}"
        )


rule align_locus_batch:
    input:
        manifest=RAW_FASTA_MANIFEST,
        export_complete=ancient(rules.export_retained_fastas.output.completion),
    output:
        completion="results/loci/batches/alignment_batch_{batch_id}.complete",
    params:
        output_dir=ALIGNMENT_DIR,
        log_dir="results/loci/logs/mafft",
        batch_size=int(config.get("alignment_batch_size", 200)),
        sequence_type=SEQUENCE_TYPE,
    threads:
        get_thread_count("alignment")
    conda:
        "../envs/alignment.yaml"
    shell:
        (
            "python3 -m scripts.run_alignment_batch "
            "--manifest {input.manifest:q} "
            "--output-dir {params.output_dir:q} "
            "--log-dir {params.log_dir:q} "
            "--batch-id {wildcards.batch_id:q} "
            "--batch-size {params.batch_size} "
            "--threads-per-alignment {threads} "
            "--sequence-type {params.sequence_type:q} "
            "&& touch {output.completion:q}"
        )
