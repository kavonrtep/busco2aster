# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

`busco2aster` is a Snakemake phylogenomics workflow: genome assemblies → BUSCO orthologs → gene trees (IQ-TREE 3) → coalescent species tree (ASTER `wastral`). Every stage produces explicit intermediate tables, QC outputs, and a final Markdown + HTML audit report.

Design source of truth: `docs/problem_formulation.md` (biology) and `docs/implementation_consolidated.md` (architecture). Current v1 scope: shared BUSCO lineage, complete single-copy loci only, protein sequences by default.

## Commands

```bash
# Dry-run (always validate first)
snakemake -n -p --cores 4 results/report/report.html

# Install external tools (IQ-TREE 3 and ASTER into work/tools/)
python3 -m scripts.install_iqtree3
python3 -m scripts.install_aster

# Run individual stages
snakemake --use-conda --conda-frontend conda --cores 4 results/qc/busco_summary.tsv
snakemake --use-conda --conda-frontend conda --cores 4 results/qc/retained_loci.tsv
snakemake --use-conda --conda-frontend conda --cores 4 results/loci/alignments.complete
snakemake --cores 4 results/gene_trees/gene_trees.complete
snakemake --cores 4 results/species_tree/species_tree.complete
snakemake --cores 4 results/concordance/gcf.complete results/concordance/scfl.complete results/concordance/wastral_quartets.complete
snakemake --cores 4 results/topology_tests/topology_tests.complete
snakemake --use-conda --conda-frontend conda --cores 4 results/report/report.html

# Unit tests
python3 -m unittest discover -s tests -v

# Container build and run
apptainer build busco2aster.sif busco2aster.def
apptainer run -B "$(pwd)" busco2aster.sif -c config/config.yaml --repo-root "$(pwd)" --directory "$(pwd)" -t 4 --target results/report/report.html

# Normalize a sample manifest from CSV
python3 -m scripts.normalize_manifest --input test_data/genome_set2.csv --output /tmp/samples.tsv

# Post-run cleanup (always use --dry-run first)
python3 -m scripts.cleanup_outputs --mode final_report --dry-run
```

Rules requiring Conda environments: manifest validation, BUSCO, assembly prep, DNA extraction, alignment, and the HTML report step. Gene tree and species tree rules use tools from PATH or `work/tools/`.

## Architecture

### Pipeline Stages (in order)

1. **Manifest validation** (`workflow/rules/manifest.smk`) — normalize `config/samples.tsv`, verify assembly paths, produce `results/metadata/samples.validated.tsv` and `taxon_name_map.tsv`
2. **Assembly prep** (`workflow/rules/assembly_prep.smk`) — rewrite assemblies with `seqkit` into wrapped gzipped FASTA cache under `work/assemblies_prepared/`
3. **BUSCO** (`workflow/rules/busco.smk`) — genome-mode BUSCO with shared lineage; raw scratch in `work/busco/`, stable sequence files copied to `results/busco/*/busco_sequences/`
4. **BUSCO QC** (`workflow/rules/busco_summary.smk`) — parse BUSCO outputs into `results/qc/busco_summary.tsv` and `busco_records.tsv`
5. **Locus selection** (`workflow/rules/locus_matrix.smk`) — occupancy matrix → filter loci ≥ threshold (default 0.8) → `results/qc/retained_loci.tsv`
6. **DNA extraction** (`workflow/rules/dna_extract.smk`, optional) — extract CDS sequences via BUSCO GFF + `gffread`
7. **FASTA export + alignment** (`workflow/rules/alignment.smk`) — one FASTA per retained locus, batched MAFFT alignment (default 200 loci/batch), output to `results/loci/alignments/`
8. **Gene trees** (`workflow/rules/gene_trees.smk`) — single IQ-TREE 3 directory-mode run over all alignments; `abayes` + UFBoot; output `results/gene_trees/gene_trees.raw.tre`
9. **Species tree** (`workflow/rules/species_tree.smk`) — `wastral` (default) or `astral4` fallback; output `results/species_tree/species_tree.wastral.tre`
10. **Concordance** (`workflow/rules/concordance.smk`) — gCF, sCFL, and wASTRAL quartet scores on the species tree
11. **Topology tests** (`workflow/rules/topology_tests.smk`) — re-run wASTRAL with `-u 2` to annotate branch posteriors; NNI alternatives at contested branches (pp1 < threshold); IQ-TREE AU test on concatenated supermatrix; outputs under `results/topology_tests/`; gated by `run_topology_tests` config key
12. **Reports** (`workflow/rules/report.smk`, `visual_report.smk`) — Markdown audit report + Quarto/R HTML report with plots, concordance panels, and topological uncertainty section

### Key Files

| Path | Purpose |
|------|---------|
| `Snakefile` | Top-level entrypoint; imports all rule modules |
| `config/config.yaml` | Runtime settings (lineage, occupancy, threads, tool paths) |
| `config/samples.tsv` | Input manifest: `sample_id`, `taxon_id`, `assembly_fasta` |
| `config/config.template.yaml` | Minimal starting template for new runs |
| `run_pipeline.py` | Container entrypoint wrapper; validates inputs, injects tool paths, auto-scales threads |
| `workflow/rules/` | One `.smk` file per stage; `topology_tests.smk` is the newest |
| `scripts/` | Python helper modules; `quartet_support.py`, `generate_alternatives.py`, `build_supermatrix.py`, `parse_au_test.py`, `topology_tests.py` are the topology-test helpers |
| `workflow/envs/` | Conda environment specs (one per stage group) |
| `reports/template.html.j2` | Jinja2 template for the self-contained HTML visual report |
| `busco2aster.def` | Apptainer/Singularity image definition |
| `.github/workflows/build-sif.yaml` | CI: builds `.sif` on version tags, runs smoke tests, uploads to release |

### Design Decisions to Preserve

- **Complete single-copy loci only** in v1 (ASTRAL-Pro for duplicates is out of scope)
- **wASTRAL** (branch-support-weighted) is the default; `astral4` available for topology-only runs
- **IQ-TREE directory mode** — one run for all loci, not one job per locus
- **Batched alignment** — reduces Snakemake scheduler overhead
- **Explicit QC tables at every stage** — no hidden heuristics
- **Topology tests use existing tools only** — wASTRAL (`-u 2`) and IQ-TREE 3 (AU test); no new conda environments or binaries; `run_topology_tests: false` skips the entire extension
- **`run_pipeline.py` `SINGLE_JOB_RULES`** must include `wastral_u2_annotate` and `au_topology_test` so the container wrapper scales their thread budget to the full core allocation

## Coding Style

- Rule names: `snake_case`, one biological stage per rule (e.g., `run_busco`, `build_locus_matrix`)
- Config keys: `snake_case`
- Logic location: workflow rules in `workflow/rules/`, helper Python in `scripts/`, report templates in `reports/`
- Generated outputs (`results/`, `work/`) are out of version control

## Commit Style

Short imperative subjects scoped to one stage: `Add BUSCO QC summary tables`, `Fix locus matrix occupancy filter`. PR descriptions should state which phase changed, how it was validated (dry-run, unit tests, targeted run), and whether QC criteria or output schemas changed.