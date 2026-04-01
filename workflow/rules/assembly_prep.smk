rule prepare_assembly:
    input:
        assembly=lambda wildcards: SAMPLE_TO_ASSEMBLY[wildcards.sample],
    output:
        prepared=ASSEMBLY_PREPARED_PATTERN,
        qc=ASSEMBLY_PREP_QC_PATTERN,
    params:
        wrap_width=ASSEMBLY_WRAP_WIDTH,
    threads:
        get_thread_count("assembly_prep")
    conda:
        "../envs/assembly_prep.yaml"
    shell:
        (
            "python3 -m scripts.prepare_assembly "
            "--input {input.assembly:q} "
            "--output {output.prepared:q} "
            "--qc-output {output.qc:q} "
            "--wrap-width {params.wrap_width} "
            "--threads {threads}"
        )
