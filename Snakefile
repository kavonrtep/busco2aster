"""
Phase 8 Snakemake entrypoint.

This file currently exposes manifest validation, BUSCO tool preflight,
per-sample BUSCO execution, BUSCO QC summarization, locus selection,
batched retained-locus FASTA export, alignment, per-locus gene trees,
and ASTER species-tree inference.
"""

import csv
from pathlib import Path

from scripts.alignment import load_retained_locus_ids, locus_output_paths
from scripts.busco import busco_output_paths
from scripts.gene_trees import gene_tree_output_paths
from scripts.species_tree import species_tree_output_paths

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
GENE_TREE_DIR = "results/gene_trees/per_locus"
GENE_TREE_MANIFEST = "results/gene_trees/gene_tree_manifest.tsv"
GENE_TREE_AGGREGATE = "results/gene_trees/gene_trees.raw.tre"
WASTRAL_GENE_TREE_INPUT = "results/gene_trees/gene_trees.wastral.tre"
GENE_TREES_COMPLETE = "results/gene_trees/gene_trees.complete"
SPECIES_TREE_DIR = "results/species_tree"
SPECIES_TREE_BACKEND = str(config.get("species_tree_backend", "wastral"))
DEFAULT_SPECIES_TREE_OUTPUTS = species_tree_output_paths(SPECIES_TREE_BACKEND)
WASTRAL_OUTPUTS = species_tree_output_paths("wastral")
ASTRAL4_OUTPUTS = species_tree_output_paths("astral4")
SPECIES_TREE_COMPLETE = f"{SPECIES_TREE_DIR}/species_tree.complete"
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

def retained_locus_ids() -> list[str]:
    retained_output = checkpoints.select_loci.get().output[0]
    return load_retained_locus_ids(Path(retained_output))

def retained_alignment_targets(wildcards):
    return [locus_output_paths(locus_id)["alignment"] for locus_id in retained_locus_ids()]

def retained_gene_tree_reports(wildcards):
    return [gene_tree_output_paths(locus_id)["report"] for locus_id in retained_locus_ids()]


def retained_gene_tree_treefiles(wildcards):
    return [gene_tree_output_paths(locus_id)["treefile"] for locus_id in retained_locus_ids()]

include: "workflow/rules/manifest.smk"
include: "workflow/rules/busco.smk"
include: "workflow/rules/busco_summary.smk"
include: "workflow/rules/locus_matrix.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/gene_trees.smk"
include: "workflow/rules/species_tree.smk"


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

rule gene_trees_complete:
    input:
        [
            ALIGNMENTS_COMPLETE,
            GENE_TREE_MANIFEST,
            GENE_TREE_AGGREGATE,
            retained_gene_tree_treefiles,
        ]
    output:
        touch(GENE_TREES_COMPLETE)


rule species_tree_complete:
    input:
        [
            GENE_TREES_COMPLETE,
            DEFAULT_SPECIES_TREE_OUTPUTS["treefile"],
            DEFAULT_SPECIES_TREE_OUTPUTS["log"],
            DEFAULT_SPECIES_TREE_OUTPUTS["completion"],
        ]
    output:
        touch(SPECIES_TREE_COMPLETE)

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
            GENE_TREE_MANIFEST,
            GENE_TREE_AGGREGATE,
            GENE_TREES_COMPLETE,
            WASTRAL_GENE_TREE_INPUT,
            DEFAULT_SPECIES_TREE_OUTPUTS["treefile"],
            DEFAULT_SPECIES_TREE_OUTPUTS["log"],
            SPECIES_TREE_COMPLETE,
        ]
