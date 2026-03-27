"""
Phase 1 Snakemake entrypoint.

This file currently exposes only the manifest validation stage so the workflow
has a real executable baseline before heavier biological stages are added.
"""

from pathlib import Path

configfile: "config/config.yaml"

REPO_ROOT = Path(workflow.basedir).resolve()
VALIDATED_MANIFEST = "results/metadata/samples.validated.tsv"
TAXON_NAME_MAP = "results/metadata/taxon_name_map.tsv"

include: "workflow/rules/manifest.smk"

localrules: all

rule all:
    input:
        [
            VALIDATED_MANIFEST,
            TAXON_NAME_MAP,
        ]
