"""
Phase 6 Snakemake entrypoint.

This file currently exposes manifest validation, BUSCO tool preflight,
per-sample BUSCO execution, BUSCO QC summarization, locus selection,
batched retained-locus FASTA export, and alignment.
"""

import csv
from pathlib import Path

from scripts.alignment import load_retained_locus_ids, locus_output_paths
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
BUSCO_SUMMARY_TABLE = "results/qc/busco_summary.tsv"
BUSCO_RECORDS_TABLE = "results/qc/busco_records.tsv"
LOCUS_TAXON_MATRIX = "results/qc/locus_taxon_matrix.tsv"
RETAINED_LOCI_TABLE = "results/qc/retained_loci.tsv"
RAW_FASTA_DIR = "results/loci/raw_fastas"
RAW_FASTA_MANIFEST = "results/loci/raw_fastas_manifest.tsv"
ALIGNMENT_DIR = "results/loci/alignments"
ALIGNMENTS_COMPLETE = "results/loci/alignments.complete"
SAMPLE_RECORDS = load_sample_records(config["samples"])
SAMPLES = [row["sample_id"] for row in SAMPLE_RECORDS]
SAMPLE_TO_ASSEMBLY = {row["sample_id"]: row["assembly_fasta"] for row in SAMPLE_RECORDS}
BUSCO_OUTPUTS = {sample: busco_output_paths(sample) for sample in SAMPLES}
BUSCO_COMPLETION_TARGETS = [BUSCO_OUTPUTS[sample]["completion"] for sample in SAMPLES]
BUSCO_SUMMARY_INPUTS = [
    BUSCO_OUTPUTS[sample][artifact]
    for sample in SAMPLES
    for artifact in ("completion", "short_summary", "full_table", "paths")
]

include: "workflow/rules/manifest.smk"
include: "workflow/rules/busco.smk"
include: "workflow/rules/busco_summary.smk"
include: "workflow/rules/locus_matrix.smk"
include: "workflow/rules/alignment.smk"


def retained_locus_ids() -> list[str]:
    retained_output = checkpoints.select_loci.get().output[0]
    return load_retained_locus_ids(Path(retained_output))

def retained_alignment_targets(wildcards):
    return [locus_output_paths(locus_id)["alignment"] for locus_id in retained_locus_ids()]

localrules: all

rule alignments_complete:
    input:
        [
            BUSCO_SUMMARY_TABLE,
            BUSCO_RECORDS_TABLE,
            LOCUS_TAXON_MATRIX,
            RETAINED_LOCI_TABLE,
            RAW_FASTA_MANIFEST,
            retained_alignment_targets,
        ]
    output:
        touch(ALIGNMENTS_COMPLETE)

rule all:
    input:
        [
            VALIDATED_MANIFEST,
            TAXON_NAME_MAP,
            BUSCO_TOOL_VERSIONS,
            BUSCO_DATASETS,
            BUSCO_LINEAGE_VERIFIED,
            *BUSCO_COMPLETION_TARGETS,
            BUSCO_SUMMARY_TABLE,
            BUSCO_RECORDS_TABLE,
            LOCUS_TAXON_MATRIX,
            RETAINED_LOCI_TABLE,
            RAW_FASTA_MANIFEST,
            retained_alignment_targets,
        ]
