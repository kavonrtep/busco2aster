# Consolidated Implementation Plan

## Purpose
This document consolidates the original implementation draft, your inline comments, and the current upstream tool documentation into one working plan for v1 of the phylogenomics workflow.

The workflow goal stays the same: start from one genome assembly per taxon, identify comparable BUSCO loci, infer per-locus gene trees, and estimate an unrooted species tree with strong QC and full traceability.

## Adopted Decisions
The following points are now treated as settled for v1:

- Inputs are only the genome FASTA files listed in [`test_data/genome_set2.csv`](/home/petr/PycharmProjects/get_phylo/test_data/genome_set2.csv).
- Everything under `test_data/*/analysis/` is out of scope for v1 and should be ignored by the workflow.
- Development should move to a new normalized sample sheet, `config/samples.tsv`, with only `sample_id`, `taxon_id`, and `assembly_fasta` for now.
- The original `test_data/genome_set2.csv` stays as a reference input example, not the long-term internal manifest.
- BUSCO lineage verification must be an explicit pipeline step. For the current Solanum test dataset, the intended lineage family is Solanales, but the exact BUSCO identifier must be validated against the installed BUSCO dataset names before execution.
- Length-dispersion and similar sequence pathology metrics should be reported in v1, not used as hard exclusion thresholds yet.
- Conda environments should be small and rule-specific so each rule can be developed and tested independently.

## Tool Choice Decision
### Gene-tree inference
Use IQ-TREE 3 as the default gene-tree engine. The upstream project currently publishes official binaries for Linux, macOS, and Windows, so the implementation should not depend on Conda availability for this tool.

### Species-tree inference
Use the ASTER package as the species-tree toolkit, with `wastral` as the default backend for v1.

This choice is based on current upstream documentation:

- The legacy Java ASTRAL repository explicitly points users to ASTER as the newer C implementation and says, in effect, to use the new code.
- The ASTER README describes Weighted ASTRAL as operating on unrooted gene trees with branch lengths and/or support values.
- The ASTER README FAQ recommends Weighted ASTRAL when the input gene trees include branch lengths and bootstrap or Bayesian support and when transfer or hybridization is not the main concern.
- The ASTER tutorial shows `wastral` as the execution path for this use case and documents `conda install aster` as an installation route.

This is an inference from those sources: because our planned inputs are single-copy BUSCO loci and our gene trees will come from IQ-TREE with branch lengths and support values, `wastral` is a better default fit than legacy ASTRAL-III. `astral4` should still be exposed as a fallback backend if we need a pure topology-only comparison. `astral-pro3` is the future path if we later decide to admit duplicated loci.

## Consequence of Choosing `wastral`
The old plan to collapse low-support branches before species-tree inference was tailored to ASTRAL-III. That is no longer the default path.

For v1:

- the primary species-tree rule should consume raw supported gene trees
- branch support values must be preserved in a format that `wastral -B` can use
- a collapsed-gene-tree branch can still exist as an optional diagnostic or `astral4` fallback path

This is also an inference from the ASTER tutorial: `wastral` uses branch supports directly, while ASTRAL-III and ASTRAL-IV ignore them.

## Proposed Repository Layout
The repository should grow into this structure:

```text
docs/
  problem_formulation.md
  implementation.md
  implementation_consolidated.md
config/
  config.yaml
  samples.tsv
workflow/
  rules/
  envs/
scripts/
tests/
results/            # untracked
work/               # untracked
Snakefile
```

## Configuration
`config/config.yaml` should hold run-wide settings:

```yaml
busco_lineage: null
occupancy_threshold: 0.8
sequence_type: protein
alignment_tool: mafft
trimming_mode: none
gene_tree_tool: iqtree3
species_tree_backend: wastral
threads:
  busco: 8
  iqtree: 4
```

`config/samples.tsv` should have exactly these columns in v1:

```text
sample_id    taxon_id    assembly_fasta
```

Taxon names should be sanitized once, recorded in `results/metadata/taxon_name_map.tsv`, and reused everywhere downstream.

## Workflow Stages
1. `normalize_manifest`
   Convert or replace the current reference CSV with `config/samples.tsv`.
2. `validate_manifest`
   Resolve input paths, check file existence, and generate the taxon-name map.
3. `verify_busco_lineage`
   Check that the configured lineage exists in the current BUSCO installation and fail fast otherwise.
4. `run_busco`
   Run BUSCO genome mode per sample into `results/busco/{sample_id}/`.
5. `summarize_busco`
   Parse BUSCO summaries and build `results/qc/busco_summary.tsv`.
6. `build_locus_matrix`
   Build `results/qc/locus_taxon_matrix.tsv` with status, sequence paths, sequence lengths, translated lengths, stop-codon flags, invalid-character flags, and QC notes.
7. `select_loci`
   Apply only the hard v1 rules: occupancy threshold and complete single-copy status. Record length dispersion and other warnings, but do not filter on them yet.
8. `export_locus_fastas`
   Write one retained protein FASTA per locus under `results/loci/raw_fastas/`.
9. `align_loci`
   Align each locus with MAFFT into `results/loci/alignments/`.
10. `infer_gene_trees`
    Run IQ-TREE 3 per locus, preserving branch lengths and support values needed by `wastral`.
11. `run_species_tree`
    Run `wastral` as the default backend and write `results/species_tree/species_tree.wastral.tre`.
12. `optional_astral4_compare`
    Optional comparison path using topology-only trees and the `astral4` backend.
13. `concordance_and_report`
    Compute downstream summaries, render `results/report/report.md`, and render
    the visual HTML report at `results/report/report.html`.

## Filtering Policy
The v1 filter logic should stay strict but minimal:

- keep loci with occupancy `>= 0.8`
- keep only complete single-copy BUSCO hits
- report duplicated, fragmented, missing, stop-codon, and length-dispersion signals
- defer hard QC cutoffs on dispersion until after inspecting the current test dataset

This keeps the first implementation inspectable and avoids burying biological judgment inside silent thresholds.

## Environment Strategy
Build environments gradually and per rule. Each major rule or tight rule group should get its own small environment file under `workflow/envs/`. The intended workflow is:

1. define a rule
2. define the smallest environment needed for that rule
3. test that rule in isolation
4. move to the next rule

This avoids rebuilding a large all-in-one environment and lets Snakemake support incremental development.

## Outputs
The stable deliverables for v1 should be:

- `results/qc/busco_summary.tsv`
- `results/qc/locus_taxon_matrix.tsv`
- `results/qc/retained_loci.tsv`
- `results/loci/alignments/`
- `results/gene_trees/gene_trees.raw.tre`
- `results/species_tree/species_tree.wastral.tre`
- `results/species_tree/species_tree.wastral.log`
- `results/report/report.md`
- `results/report/report.html`

Optional comparison outputs can include:

- `results/gene_trees/gene_trees.collapsed.tre`
- `results/species_tree/species_tree.astral4.tre`

## Remaining Checks
These are implementation checks, not open design questions:

- verify the exact BUSCO lineage identifier string for the Solanum dataset
- confirm the exact IQ-TREE 3 support-generation command we want to standardize on for `wastral -B`
- verify the exact ASTER command-line options we want to capture in the Snakemake rule wrappers

## Sources
- ASTRAL repository: <https://github.com/smirarab/ASTRAL>
- ASTER repository: <https://github.com/chaoszhang/ASTER>
- ASTER tutorial PDF: <https://tandy.cs.illinois.edu/ASTER-Tutorial-Evolution2025.pdf>
- IQ-TREE project site: <https://iqtree.github.io/>
