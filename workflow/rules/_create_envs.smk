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


rule create_env_alignment:
    output:
        "/tmp/busco2aster_env_alignment"
    conda:
        "../envs/alignment.yaml"
    shell:
        "touch {output}"
