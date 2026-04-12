# busco2aster

`busco2aster` is a Snakemake workflow for building an unrooted species tree from genome assemblies using BUSCO-defined ortholog markers. It is designed to keep the run inspectable: every major stage writes explicit intermediate tables, stable aggregate outputs, a machine-readable Markdown audit report, and a visual HTML report.

## What the Workflow Does

Given one genome assembly per taxon, the workflow:

1. validates a normalized sample manifest
2. normalizes each assembly into a wrapped gzipped cache for BUSCO
3. runs BUSCO in genome mode with one shared lineage
4. parses BUSCO outputs into sample-level and locus-level QC tables
5. retains loci that meet the v1 filter policy
6. exports one protein FASTA per retained locus
7. aligns loci with batched MAFFT jobs
8. infers gene trees with IQ-TREE 3 directory mode
9. infers an unrooted species tree with ASTER `wastral`
10. computes concordance metrics on the final species tree (gCF, sCF, wASTRAL quartet posteriors)
11. identifies contested branches and generates alternative topology hypotheses
12. tests whether alternative topologies can be statistically rejected with the IQ-TREE AU test
13. renders final Markdown and HTML reports

Current v1 defaults in [config/config.yaml](/home/petr/PycharmProjects/get_phylo/config/config.yaml):

- BUSCO lineage: `solanales_odb12`
- occupancy threshold: `0.8`
- sequence type: amino acid
- gene-tree support mode: `abayes`
- species-tree backend: `wastral`
- default thread budget in config: `4`
- topology test threshold: `0.95` (branches with pp1 below this are contested)
- AU test replicates: `10000`

## Repository Layout

- [config/samples.tsv](/home/petr/PycharmProjects/get_phylo/config/samples.tsv): normalized input manifest
- [config/config.yaml](/home/petr/PycharmProjects/get_phylo/config/config.yaml): workflow settings
- [config/config.template.yaml](/home/petr/PycharmProjects/get_phylo/config/config.template.yaml): minimal starting template for new runs
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

Executable path settings are optional. If `iqtree_executable`,
`wastral_executable`, or `astral4_executable` are omitted from the config, the
workflow resolves `iqtree3`, `wastral`, and `astral4` from `PATH`. In container
mode, [run_pipeline.py](/home/petr/PycharmProjects/get_phylo/run_pipeline.py) injects the
container-internal tool paths automatically.

Current optimization-related behavior:

- assembly FASTA files are rewritten with `seqkit` into wrapped cached `.fa.gz` files before BUSCO
- BUSCO raw scratch is written under `work/busco/`
- stable BUSCO sequence files are copied into `results/busco/*/busco_sequences/`
- alignments run in batches to reduce scheduler overhead
- gene trees use one IQ-TREE directory-mode run instead of one workflow job per locus

## Input Format

The workflow reads [config/samples.tsv](/home/petr/PycharmProjects/get_phylo/config/samples.tsv) with three required columns:

```tsv
sample_id	taxon_id	assembly_fasta
```

Start from the minimal config template:

```bash
cp config/config.template.yaml config/my_run.yaml
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
snakemake --cores 4 results/topology_tests/topology_tests.complete
snakemake --use-conda --conda-frontend conda --cores 4 results/report/report.html
```

For a dry-run:

```bash
snakemake -n -p --cores 4 results/report/report.html
```

Optional cleanup after a finished run:

```bash
python3 -m scripts.cleanup_outputs --mode resume_alignments --dry-run
python3 -m scripts.cleanup_outputs --mode resume_gene_trees --dry-run
python3 -m scripts.cleanup_outputs --mode final_report --dry-run
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
- defaults to the workflow `all` target if `--target` is omitted
- uses all visible CPUs by default if `-t/--threads` is omitted
- auto-scales per-rule thread settings from the global core budget unless `thread_policy: fixed`
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
- `results/busco/*/busco_sequences/`
- `results/gene_trees/gene_tree_manifest.tsv`
- `results/gene_trees/gene_trees.raw.tre`
- `results/species_tree/species_tree.wastral.tre`
- `results/species_tree/species_tree.wastral.log`
- `results/topology_tests/branch_quartet_support.tsv`
- `results/topology_tests/contested_branches.tsv`
- `results/topology_tests/candidate_trees.tre`
- `results/topology_tests/candidate_trees_manifest.tsv`
- `results/topology_tests/supermatrix.phy`
- `results/topology_tests/au_test_results.tsv`
- [results/report/report.md](/home/petr/PycharmProjects/get_phylo/results/report/report.md)
- [results/report/report.html](/home/petr/PycharmProjects/get_phylo/results/report/report.html)

The final species tree is unrooted. The Markdown report is the compact audit summary. The HTML report adds dataset overview plots, a pipeline diagram, a labeled species tree, branch-level concordance panels, conflict summaries, gene-tree heterogeneity plots, and a topological uncertainty section (quartet posteriors + AU test results).

## Topology Testing

After the species tree is inferred, the workflow identifies branches where the
topology is uncertain and formally tests whether alternative topologies can be
statistically rejected.

**Stage A — Contested branch identification**: wASTRAL is re-run with `-u 2`
to annotate every internal branch with local posterior probabilities (pp1/pp2/pp3)
under the multispecies coalescent. Branches with `pp1 < contested_branch_threshold`
(default `0.95`) are flagged as contested.

**Stage B — Alternative topology generation**: For each contested branch, an NNI
swap produces one alternative species tree. All 2^k − 1 combinations of swaps are
generated for k contested branches (capped at `max_contested_branches`, default 4).
User-supplied hypothesis trees can be appended via `hypothesis_trees` in the config.

**Stage C — AU topology test**: IQ-TREE evaluates each candidate tree on a
concatenated supermatrix with the approximately unbiased (AU) test. Trees with
p-AU ≥ 0.05 belong to the 95% confidence set and cannot be statistically rejected.

To skip the extension, set `run_topology_tests: false` in the config. To run
only the topology test stage:

```bash
snakemake --cores 4 results/topology_tests/topology_tests.complete
```

Relevant config keys:

| Key | Default | Description |
|-----|---------|-------------|
| `run_topology_tests` | `true` | Enable/disable the extension |
| `contested_branch_threshold` | `0.95` | pp1 below which a branch is contested |
| `max_contested_branches` | `4` | Cap on branches used for NNI combinations |
| `hypothesis_trees` | `null` | Optional Newick file of additional hypothesis trees |
| `au_test_replicates` | `10000` | RELL replicates for the AU test |
| `au_test_model` | `from_gene_trees` | Model for supermatrix partitions; `from_gene_trees` reuses IQ-TREE selected models |

The HTML report includes a **Topological uncertainty** section with a per-branch
quartet posterior table, an AU test results table, and an interpretation guide.

## Testing and Development

Useful validation commands:

```bash
snakemake -n
snakemake -n -p --cores 4 results/loci/alignments.complete
snakemake -n -p --cores 4 results/gene_trees/gene_trees.raw.tre
snakemake -s workflow/Snakefile_create_envs -n
python3 run_pipeline.py --config config/config.yaml --repo-root . --directory . --target results/metadata/samples.validated.tsv --snakemake-args="--dry-run"
python3 -m scripts.cleanup_outputs --mode final_report --dry-run
python3 -m unittest discover -s tests -v
git status --short
```

Design notes live in [docs/problem_formulation.md](/home/petr/PycharmProjects/get_phylo/docs/problem_formulation.md), [docs/implementation_consolidated.md](/home/petr/PycharmProjects/get_phylo/docs/implementation_consolidated.md), [docs/implementation_phases.md](/home/petr/PycharmProjects/get_phylo/docs/implementation_phases.md), [docs/containerization_plan.md](/home/petr/PycharmProjects/get_phylo/docs/containerization_plan.md), [docs/visual_report_plan.md](/home/petr/PycharmProjects/get_phylo/docs/visual_report_plan.md), [docs/optimization_review.md](/home/petr/PycharmProjects/get_phylo/docs/optimization_review.md), [docs/optimization_plan.md](/home/petr/PycharmProjects/get_phylo/docs/optimization_plan.md), and [docs/iqtree_directory_mode_evaluation.md](/home/petr/PycharmProjects/get_phylo/docs/iqtree_directory_mode_evaluation.md).
