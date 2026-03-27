rule validate_manifest:
    input:
        manifest=config["samples"],
    output:
        validated=VALIDATED_MANIFEST,
        taxon_map=TAXON_NAME_MAP,
    params:
        repo_root=str(REPO_ROOT),
    shell:
        (
            "python3 -m scripts.validate_manifest "
            "--input {input.manifest} "
            "--validated-output {output.validated} "
            "--taxon-map-output {output.taxon_map} "
            "--repo-root {params.repo_root}"
        )

