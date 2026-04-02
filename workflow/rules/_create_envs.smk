"""
Dummy rules for pre-creating conda environments during container build.

Each rule references its env YAML via ../envs/, matching the same relative path
used by the real workflow rules so that Snakemake computes the identical env
hash inside and outside the container.
"""


rule create_env_busco:
    output:
        "/tmp/busco2aster_env_busco"
    conda:
        "../envs/busco.yaml"
    shell:
        "touch {output}"


rule create_env_assembly_prep:
    output:
        "/tmp/busco2aster_env_assembly_prep"
    conda:
        "../envs/assembly_prep.yaml"
    shell:
        "seqkit version >/dev/null && touch {output}"


rule create_env_alignment:
    output:
        "/tmp/busco2aster_env_alignment"
    conda:
        "../envs/alignment.yaml"
    shell:
        "touch {output}"


rule create_env_dna_extract:
    output:
        "/tmp/busco2aster_env_dna_extract"
    conda:
        "../envs/dna_extract.yaml"
    shell:
        "gffread --version >/dev/null 2>&1 && touch {output}"


rule create_env_report:
    output:
        "/tmp/busco2aster_env_report"
    conda:
        "../envs/report.yaml"
    shell:
        (
            "quarto --version >/dev/null && "
            "Rscript -e 'library(ggtree); library(phangorn)' >/dev/null && "
            "touch {output}"
        )
