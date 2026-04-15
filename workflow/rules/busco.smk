rule capture_busco_metadata:
    output:
        versions=BUSCO_TOOL_VERSIONS,
        datasets=BUSCO_DATASETS,
    conda:
        "../envs/busco.yaml"
    params:
        download_path=BUSCO_DOWNLOAD_PATH,
    shell:
        r"""
        mkdir -p "$(dirname {output.versions})"
        printf "tool\tversion\n" > {output.versions}
        printf "busco\t%s\n" "$(busco --version | tr -d '\r')" >> {output.versions}
        busco --download_path {params.download_path:q} --list-datasets > {output.datasets}
        """


rule verify_busco_lineage:
    input:
        datasets=BUSCO_DATASETS,
    output:
        BUSCO_LINEAGE_VERIFIED,
    params:
        lineage=BUSCO_LINEAGE,
    shell:
        (
            "python3 -m scripts.verify_busco_lineage "
            "--dataset-list {input.datasets} "
            "--lineage {params.lineage} "
            "--output {output}"
        )


rule fetch_busco_lineage:
    """Pre-download the configured BUSCO lineage once so parallel run_busco jobs don't race the same tarball."""
    input:
        verified=ancient(BUSCO_LINEAGE_VERIFIED),
    output:
        dataset_cfg=BUSCO_LINEAGE_DATASET_CFG,
    params:
        lineage=BUSCO_LINEAGE,
        download_path=BUSCO_DOWNLOAD_PATH,
    conda:
        "../envs/busco.yaml"
    shell:
        r"""
        mkdir -p {params.download_path:q}
        busco --download {params.lineage:q} --download_path {params.download_path:q}
        test -f {output.dataset_cfg:q}
        """


rule run_busco:
    input:
        assembly=ASSEMBLY_PREPARED_PATTERN,
        validated=ancient(VALIDATED_MANIFEST),
        lineage=ancient(BUSCO_LINEAGE_VERIFIED),
        lineage_dataset=ancient(BUSCO_LINEAGE_DATASET_CFG),
    output:
        command="results/busco/{sample}/command.sh",
        paths="results/busco/{sample}/paths.tsv",
        sequences=directory("results/busco/{sample}/busco_sequences"),
        short_summary="results/busco/{sample}/short_summary.txt",
        full_table="results/busco/{sample}/full_table.tsv",
        completion="results/busco/{sample}/run.complete",
    params:
        sample_dir="results/busco/{sample}",
        raw_root=lambda wildcards: f"{BUSCO_WORK_ROOT}/{wildcards.sample}/raw",
        lineage=BUSCO_LINEAGE,
        download_path=BUSCO_DOWNLOAD_PATH,
    threads:
        get_thread_count("busco")
    conda:
        "../envs/busco.yaml"
    shell:
        r"""
        mkdir -p {params.sample_dir:q}
        mkdir -p {params.raw_root:q}
        printf '#!/usr/bin/env bash\nset -euo pipefail\n' > {output.command:q}
        printf '%s\n' "busco --in {input.assembly:q} --mode genome --lineage_dataset {params.lineage:q} --cpu {threads} --out {wildcards.sample:q} --out_path {params.raw_root:q} --download_path {params.download_path:q} -f" >> {output.command:q}
        bash {output.command:q}
        python3 -m scripts.standardize_busco_run --sample-id {wildcards.sample:q} --sample-dir {params.sample_dir:q} --raw-root {params.raw_root:q}
        """
