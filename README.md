# busco2aster

`busco2aster` is a Snakemake workflow for building an unrooted species tree from genome assemblies using BUSCO-defined ortholog markers. It is designed to keep the run inspectable: every major stage writes explicit intermediate tables, per-locus files, a machine-readable Markdown audit report, and a visual HTML report.

## What the Workflow Does

Given one genome assembly per taxon, the workflow:

1. validates a normalized sample manifest
2. runs BUSCO in genome mode with one shared lineage
3. parses BUSCO outputs into sample-level and locus-level QC tables
4. retains loci that meet the v1 filter policy
5. exports one protein FASTA per retained locus
6. aligns loci with MAFFT
7. infers one gene tree per locus with IQ-TREE 3
8. infers an unrooted species tree with ASTER `wastral`
9. computes concordance metrics on the final species tree
10. renders final Markdown and HTML reports

Current v1 defaults in [config/config.yaml](/home/petr/PycharmProjects/get_phylo/config/config.yaml):

- BUSCO lineage: `solanales_odb12`
- occupancy threshold: `0.8`
- sequence type: amino acid
- gene-tree support mode: `abayes`
- species-tree backend: `wastral`
- default thread budget in config: `4`

## Repository Layout

- [config/samples.tsv](/home/petr/PycharmProjects/get_phylo/config/samples.tsv): normalized input manifest
- [config/config.yaml](/home/petr/PycharmProjects/get_phylo/config/config.yaml): workflow settings
- [workflow/rules/](/home/petr/PycharmProjects/get_phylo/workflow/rules): Snakemake rules by stage
- [scripts/](/home/petr/PycharmProjects/get_phylo/scripts): Python helpers and installer utilities
- [docs/](/home/petr/PycharmProjects/get_phylo/docs): design and implementation notes
- `results/`: generated workflow outputs
- `work/`: downloaded tools and transient workflow state

## Requirements

For native execution you need:

- Python 3
- Snakemake
- Conda
- `curl`
- `make` and a working C/C++ toolchain for ASTER

For containerized execution you need:

- Apptainer or Singularity

Native mode uses rule-specific Conda environments for BUSCO, MAFFT, and the
Quarto/R report step. IQ-TREE 3 and ASTER can be installed locally into
`work/tools/`. The container image bakes Snakemake, IQ-TREE 3, ASTER, Quarto,
and the rule Conda environments into a single `.sif` image.

## Input Format

The workflow reads [config/samples.tsv](/home/petr/PycharmProjects/get_phylo/config/samples.tsv) with three required columns:

```tsv
sample_id	taxon_id	assembly_fasta
```

If you start from the reference sheet in `test_data/genome_set2.csv`, normalize it with:

```bash
python3 -m scripts.normalize_manifest --input test_data/genome_set2.csv --output /tmp/samples.tsv
```

## Quick Start

### Native Mode

Install external tree-building tools:

```bash
python3 -m scripts.install_iqtree3
python3 -m scripts.install_aster
```

Verify the manifest and BUSCO lineage:

```bash
snakemake --cores 1 results/metadata/samples.validated.tsv results/metadata/taxon_name_map.tsv
snakemake --use-conda --conda-frontend conda --cores 1 results/metadata/busco_lineage_verified.tsv
```

Run the workflow in stages:

```bash
snakemake --use-conda --conda-frontend conda --cores 4 results/qc/busco_summary.tsv
snakemake --use-conda --conda-frontend conda --cores 4 results/qc/retained_loci.tsv
snakemake --use-conda --conda-frontend conda --cores 4 results/loci/alignments.complete
snakemake --cores 4 results/gene_trees/gene_trees.complete
snakemake --cores 4 results/species_tree/species_tree.complete
snakemake --cores 4 results/concordance/gcf.complete results/concordance/scfl.complete results/concordance/wastral_quartets.complete
snakemake --use-conda --conda-frontend conda --cores 4 results/report/report.html
```

For a dry-run:

```bash
snakemake -n -p --cores 4 results/report/report.html
```

### Container Mode

Build the image locally:

```bash
apptainer build busco2aster.sif busco2aster.def
```

Run the wrapper from a bind-mounted project directory:

```bash
apptainer run -B "$(pwd)" busco2aster.sif \
  -c config/config.yaml \
  --repo-root "$(pwd)" \
  --directory "$(pwd)" \
  -t 4 \
  --target results/report/report.html \
  --snakemake-args="--dry-run"
```

The same command works with `singularity run`. The wrapper:

- validates that the config, manifest, and assembly paths are visible
- prints suggested `-B` bind mounts if anything is missing
- injects container-internal paths for IQ-TREE 3 and ASTER
- runs Snakemake with `/opt/conda/envs` as the shared Conda prefix

Important runtime rule:

- keep `results/`, `work/`, `.snakemake/`, `.cache/`, and BUSCO downloads on the
  host bind-mounted directory, not inside the image

## Main Outputs

Key outputs are:

- `results/qc/busco_summary.tsv`
- `results/qc/locus_taxon_matrix.tsv`
- `results/qc/retained_loci.tsv`
- `results/loci/alignments/`
- `results/gene_trees/gene_tree_manifest.tsv`
- `results/gene_trees/gene_trees.raw.tre`
- `results/species_tree/species_tree.wastral.tre`
- `results/species_tree/species_tree.wastral.log`
- [results/report/report.md](/home/petr/PycharmProjects/get_phylo/results/report/report.md)
- [results/report/report.html](/home/petr/PycharmProjects/get_phylo/results/report/report.html)

The final species tree is unrooted. The Markdown report is the compact audit summary. The HTML report adds dataset overview plots, a pipeline diagram, a labeled species tree, branch-level concordance panels, conflict summaries, and gene-tree heterogeneity plots.

## Testing and Development

Useful validation commands:

```bash
snakemake -n
snakemake -s workflow/Snakefile_create_envs -n
python3 run_pipeline.py --config config/config.yaml --repo-root . --directory . --target results/metadata/samples.validated.tsv --snakemake-args="--dry-run"
python3 -m unittest discover -s tests -v
git status --short
```

Design notes live in [docs/problem_formulation.md](/home/petr/PycharmProjects/get_phylo/docs/problem_formulation.md), [docs/implementation_consolidated.md](/home/petr/PycharmProjects/get_phylo/docs/implementation_consolidated.md), [docs/implementation_phases.md](/home/petr/PycharmProjects/get_phylo/docs/implementation_phases.md), [docs/containerization_plan.md](/home/petr/PycharmProjects/get_phylo/docs/containerization_plan.md), and [docs/visual_report_plan.md](/home/petr/PycharmProjects/get_phylo/docs/visual_report_plan.md).
