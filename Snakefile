"""
Snakemake entrypoint for the full busco2aster workflow.

This file currently exposes manifest validation, BUSCO tool preflight,
per-sample BUSCO execution, BUSCO QC summarization, locus selection,
batched retained-locus FASTA export, batched alignment, directory-mode
IQ-TREE gene trees, ASTER species-tree inference, IQ-TREE gCF concordance
scoring, and final report generation.
"""

import csv
from pathlib import Path

from scripts.alignment import alignment_batch_output_paths, build_batch_ids, load_retained_locus_ids
from scripts.busco import busco_output_paths
from scripts.concordance import concordance_output_paths, quartet_output_paths
from scripts.gene_trees import gene_tree_directory_output_paths
from scripts.species_tree import species_tree_output_paths

configfile: "config/config.yaml"


def load_sample_records(manifest_path: str) -> list[dict[str, str]]:
    with Path(manifest_path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def get_thread_count(name: str) -> int:
    thread_config = config.get("threads", {})
    return int(thread_config.get(name, thread_config.get("default", 4)))


WORKFLOW_ROOT = Path(workflow.basedir).resolve()
REPO_ROOT = Path(str(config.get("repo_root", WORKFLOW_ROOT))).resolve()


def resolve_repo_path(path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else REPO_ROOT / path


SAMPLES_MANIFEST = resolve_repo_path(str(config["samples"])).as_posix()
BUSCO_DOWNLOAD_PATH = resolve_repo_path(str(config["busco_download_path"])).as_posix()
BUSCO_WORK_ROOT = resolve_repo_path(str(config.get("busco_work_root", "work/busco"))).as_posix()

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
GENE_TREE_MANIFEST = "results/gene_trees/gene_tree_manifest.tsv"
GENE_TREE_AGGREGATE = "results/gene_trees/gene_trees.raw.tre"
GENE_TREE_DIRECTORY_OUTPUTS = gene_tree_directory_output_paths()
WASTRAL_GENE_TREE_INPUT = "results/gene_trees/gene_trees.wastral.tre"
GENE_TREES_COMPLETE = "results/gene_trees/gene_trees.complete"
SPECIES_TREE_DIR = "results/species_tree"
SPECIES_TREE_BACKEND = str(config.get("species_tree_backend", "wastral"))
DEFAULT_SPECIES_TREE_OUTPUTS = species_tree_output_paths(SPECIES_TREE_BACKEND)
WASTRAL_OUTPUTS = species_tree_output_paths("wastral")
ASTRAL4_OUTPUTS = species_tree_output_paths("astral4")
GCF_OUTPUTS = concordance_output_paths("gcf")
SCFL_OUTPUTS = concordance_output_paths("scfl")
WASTRAL_QUARTET_OUTPUTS = quartet_output_paths()
SPECIES_TREE_COMPLETE = f"{SPECIES_TREE_DIR}/species_tree.complete"
REPORT_MARKDOWN = "results/report/report.md"
REPORT_DATA_DIR = "results/report/data"
REPORT_DATASET_SUMMARY = f"{REPORT_DATA_DIR}/dataset_summary.tsv"
REPORT_SAMPLE_QC = f"{REPORT_DATA_DIR}/sample_qc.tsv"
REPORT_LOCUS_SUMMARY = f"{REPORT_DATA_DIR}/locus_summary.tsv"
REPORT_ALIGNMENT_SUMMARY = f"{REPORT_DATA_DIR}/alignment_summary.tsv"
REPORT_BRANCH_METRICS = f"{REPORT_DATA_DIR}/branch_metrics.tsv"
REPORT_BRANCH_ALTERNATIVES = f"{REPORT_DATA_DIR}/branch_alternatives.tsv"
REPORT_GENE_TREE_HETEROGENEITY = f"{REPORT_DATA_DIR}/gene_tree_heterogeneity.tsv"
REPORT_TOPOLOGY_COUNTS = f"{REPORT_DATA_DIR}/topology_counts.tsv"
REPORT_SPECIES_TREE = f"{REPORT_DATA_DIR}/species_tree.report.tre"
REPORT_HTML = "results/report/report.html"
SAMPLE_RECORDS = load_sample_records(SAMPLES_MANIFEST)
SAMPLES = [row["sample_id"] for row in SAMPLE_RECORDS]
SAMPLE_TO_ASSEMBLY = {
    row["sample_id"]: resolve_repo_path(row["assembly_fasta"]).as_posix()
    for row in SAMPLE_RECORDS
}
BUSCO_OUTPUTS = {sample: busco_output_paths(sample) for sample in SAMPLES}
BUSCO_COMPLETION_TARGETS = [BUSCO_OUTPUTS[sample]["completion"] for sample in SAMPLES]
BUSCO_SUMMARY_INPUTS = [
    BUSCO_OUTPUTS[sample][artifact]
    for sample in SAMPLES
    for artifact in ("completion", "short_summary", "full_table", "paths")
]
ALIGNMENT_BATCH_SIZE = int(config.get("alignment_batch_size", 200))

def retained_locus_ids() -> list[str]:
    retained_output = checkpoints.select_loci.get().output[0]
    return load_retained_locus_ids(Path(retained_output))

def alignment_batch_markers(wildcards):
    batch_ids = build_batch_ids(retained_locus_ids(), ALIGNMENT_BATCH_SIZE)
    return [alignment_batch_output_paths(batch_id)["completion"] for batch_id in batch_ids]

include: "workflow/rules/manifest.smk"
include: "workflow/rules/busco.smk"
include: "workflow/rules/busco_summary.smk"
include: "workflow/rules/locus_matrix.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/gene_trees.smk"
include: "workflow/rules/species_tree.smk"
include: "workflow/rules/concordance.smk"
include: "workflow/rules/report.smk"
include: "workflow/rules/visual_report.smk"


localrules: all

rule alignments_complete:
    input:
        summary=BUSCO_SUMMARY_TABLE,
        records=BUSCO_RECORDS_TABLE,
        matrix=LOCUS_TAXON_MATRIX,
        retained=RETAINED_LOCI_TABLE,
        manifest=RAW_FASTA_MANIFEST,
        raw_fastas_complete="results/loci/raw_fastas.complete",
        batch_markers=alignment_batch_markers,
    output:
        ALIGNMENTS_COMPLETE
    params:
        alignment_dir=ALIGNMENT_DIR,
        log_dir="results/loci/logs/mafft",
    shell:
        (
            "python3 -m scripts.finalize_alignments "
            "--manifest {input.manifest:q} "
            "--output-dir {params.alignment_dir:q} "
            "--log-dir {params.log_dir:q} "
            "&& touch {output:q}"
        )

rule gene_trees_complete:
    input:
        [
            ALIGNMENTS_COMPLETE,
            GENE_TREE_MANIFEST,
            GENE_TREE_AGGREGATE,
            GENE_TREE_DIRECTORY_OUTPUTS["report"],
            GENE_TREE_DIRECTORY_OUTPUTS["treefile"],
        ]
    output:
        touch(GENE_TREES_COMPLETE)


rule species_tree_complete:
    input:
        [
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
            ALIGNMENTS_COMPLETE,
            GENE_TREE_MANIFEST,
            GENE_TREE_AGGREGATE,
            GENE_TREES_COMPLETE,
            WASTRAL_GENE_TREE_INPUT,
            DEFAULT_SPECIES_TREE_OUTPUTS["treefile"],
            DEFAULT_SPECIES_TREE_OUTPUTS["log"],
            SPECIES_TREE_COMPLETE,
            GCF_OUTPUTS["stat"],
            SCFL_OUTPUTS["stat"],
            WASTRAL_QUARTET_OUTPUTS["freqquad"],
            REPORT_MARKDOWN,
            REPORT_HTML,
        ]
