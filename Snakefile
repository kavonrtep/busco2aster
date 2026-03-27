"""
Phase 3 Snakemake entrypoint.

This file currently exposes manifest validation, BUSCO tool preflight, and
per-sample BUSCO execution.
"""

import csv
from pathlib import Path

from scripts.busco import busco_output_paths

configfile: "config/config.yaml"


def load_sample_records(manifest_path: str) -> list[dict[str, str]]:
    with Path(manifest_path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def get_thread_count(name: str) -> int:
    thread_config = config.get("threads", {})
    return int(thread_config.get(name, thread_config.get("default", 4)))


REPO_ROOT = Path(workflow.basedir).resolve()
VALIDATED_MANIFEST = "results/metadata/samples.validated.tsv"
TAXON_NAME_MAP = "results/metadata/taxon_name_map.tsv"
BUSCO_TOOL_VERSIONS = "results/metadata/busco_tool_versions.tsv"
BUSCO_DATASETS = "results/metadata/busco_datasets.txt"
BUSCO_LINEAGE_VERIFIED = "results/metadata/busco_lineage_verified.tsv"
SAMPLE_RECORDS = load_sample_records(config["samples"])
SAMPLES = [row["sample_id"] for row in SAMPLE_RECORDS]
SAMPLE_TO_ASSEMBLY = {row["sample_id"]: row["assembly_fasta"] for row in SAMPLE_RECORDS}
BUSCO_COMPLETION_TARGETS = expand("results/busco/{sample}/run.complete", sample=SAMPLES)

include: "workflow/rules/manifest.smk"
include: "workflow/rules/busco.smk"

localrules: all

rule all:
    input:
        [
            VALIDATED_MANIFEST,
            TAXON_NAME_MAP,
            BUSCO_TOOL_VERSIONS,
            BUSCO_DATASETS,
            BUSCO_LINEAGE_VERIFIED,
            *BUSCO_COMPLETION_TARGETS,
        ]
