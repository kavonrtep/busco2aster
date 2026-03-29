# Visual HTML Report Implementation Plan

## Goal
Add a visual HTML report built around the analysis chain:

`IQ-TREE locus trees -> wASTRAL species tree -> IQ-TREE gCF/sCF on the wASTRAL tree -> Quarto/R report`

The deliverable should be a stable workflow target such as
`results/report/report.html`, with the current Markdown report kept as a
machine-readable audit summary.

## Current Status
- Completed: Phases 1 through 5
- Current target: none
- Verified on the real dataset:
  - report data bundle generated under `results/report/data/`
  - Quarto HTML render completed to `results/report/report.html`
  - full test suite passed after report integration

## Key Constraints Addressed
- The report toolchain is isolated in `workflow/envs/report.yaml`, so Quarto and
  the required R packages do not need to exist in the base host environment.
- `gCF` and `sCFL` branch IDs are not a stable shared key. The implementation
  therefore joins branch metrics by canonical bipartition, not by raw branch ID.

## Recommended Report Contents
### Core sections
- Dataset overview: taxa, BUSCO completeness, retained/excluded loci, occupancy,
  alignment length distribution, and per-taxon missingness.
- Pipeline overview: one compact diagram showing BUSCO -> alignments ->
  IQ-TREE locus trees -> wASTRAL species tree -> IQ-TREE concordance ->
  Quarto report.
- Main species tree: one tree figure annotated with ASTER branch confidence and
  aligned side panels for `gCF`, `sCFL`, and quartet frequency.
- Conflict and alternative resolutions: one branch-level panel per internal edge
  showing `gCF/gDF1/gDF2/gDFP` and ASTER `q1/q2/q3`.
- Gene-tree heterogeneity: distributions of gene-tree support, topologies, and
  distance to the species tree.

### My recommendation
- Keep the main report focused on interpretation, not raw file inventory.
- Show support and concordance side by side, never as one merged score.
- Use summary plots for missingness; avoid a giant loci-by-taxon heatmap in the
  main narrative. A dense matrix can go in an appendix if needed.
- Include a short “takeaways” section up front with the weakest species-tree
  branches and the strongest alternative resolutions.

## Architecture
### Data preparation layer
Add a deterministic preprocessing step that writes report-ready tables under
`results/report/data/`. This step should normalize and join:
- BUSCO/sample QC
- retained-locus and alignment summaries
- gene-tree manifest summaries
- branch-level concordance metrics
- ASTER quartet summaries

Required output tables:
- `dataset_summary.tsv`
- `sample_qc.tsv`
- `locus_summary.tsv`
- `alignment_summary.tsv`
- `branch_metrics.tsv`
- `branch_alternatives.tsv`
- `gene_tree_heterogeneity.tsv`

### Key design rule
Build a canonical branch key from the species tree bipartitions. Then map
`gCF`, `sCFL`, and ASTER quartet outputs onto that key. This should happen in a
tested prep script, not ad hoc inside the Quarto document.

## Proposed Phases
### Phase 1: Reporting environment (Completed)
- add `workflow/envs/report.yaml`
- install `quarto`
- install `ggtree` and `treeio` through Bioconductor
- ensure the container image includes the report toolchain

Tests:
- `quarto --version`
- `Rscript -e "library(ggtree); library(phangorn)"`
- Snakemake dry-run for `results/report/report.html`

### Phase 2: Report data bundle (Completed)
- add `scripts/report_data.py`
- derive alignment sizes from `results/loci/alignments/`
- compute per-sample missingness and per-locus occupancy summaries
- compute canonical branch keys from the species tree
- map `gCF`, `sCFL`, and ASTER quartet outputs onto shared branch rows

Tests:
- unit test: canonical branch keys are stable for reordered trees
- unit test: `gCF` and `sCFL` map onto the same species-tree branches
- integration test: data bundle is generated from real workflow outputs

### Phase 3: Quarto report skeleton (Completed)
- add `reports/report.qmd`
- add a pipeline overview diagram
- add dataset overview plots and tables
- render to `results/report/report.html`

Tests:
- smoke render on a tiny fixture bundle
- workflow test: HTML render succeeds after data-prep outputs exist

### Phase 4: Tree and conflict figures (Completed)
- plot the main species tree with `ggtree`
- add aligned branch metrics: ASTER local PP, `gCF`, `sCFL`, quartet frequency
- add per-branch conflict panels with `gDF1/gDF2/gDFP` and `q1/q2/q3`

Tests:
- unit test: branch metrics table has one row per internal species-tree branch
- visual smoke test: tree plot renders for the real six-taxon dataset

### Phase 5: Gene-tree heterogeneity (Completed)
- summarize gene-tree support distribution
- summarize topology frequencies for complete-taxon gene trees
- compute distance-to-species-tree summaries with `phangorn`
- add an appendix with downloadable detailed tables

Tests:
- unit test: topology counts sum to the gene-tree subset used
- regression test: report still renders when some loci are missing taxa

## Suggested Figure Set
- BUSCO completeness by taxon
- Retained-locus occupancy histogram
- Alignment length distribution
- Per-taxon missingness bar chart
- Main annotated species tree
- Branch-level concordance vs confidence dot/bar panel
- Alternative-resolution stacked bars by branch
- Gene-tree topology frequency plot
- Gene-tree distance-to-species-tree distribution

## Current v1 Decisions
- The HTML report is fully static.
- Gene-tree heterogeneity uses complete-taxon gene trees for topology and
  species-tree distance summaries.
- The Markdown report remains the compact audit artifact beside the HTML report.
- A loci-by-taxa heatmap is deferred; the appendix currently exposes tables and
  summary plots instead.
