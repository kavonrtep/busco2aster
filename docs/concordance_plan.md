# Concordance Metrics Implementation Plan

## Purpose
Add branch-level concordance outputs to `busco2aster` in a way that matches the
current workflow structure: one fixed `wastral` species tree, one aggregate
gene-tree file, and one directory of retained per-locus alignments.

## Recommended Order

### Phase 1: Gene Concordance Factor (gCF)
Status: completed

Use IQ-TREE 3 to score the final `wastral` species tree against the aggregate
gene-tree set.

Inputs:
- `results/species_tree/species_tree.wastral.tre`
- `results/gene_trees/gene_trees.raw.tre`

Outputs:
- `results/concordance/gcf.cf.stat`
- `results/concordance/gcf.cf.tree`
- `results/concordance/gcf.cf.branch`
- `results/concordance/gcf.log`

Why first:
- uses artifacts already present in the workflow
- writes a stable tabular summary (`.cf.stat`) that is easy to test and report
- does not require any new sequence processing

Required tests:
- unit test for gCF output-path helpers
- unit test for IQ-TREE gCF command rendering
- smoke test on a tiny reference tree plus gene-tree fixture
- workflow dry-run for `results/concordance/gcf.cf.stat`

### Phase 2: Site Concordance Factor (sCFL)
Status: next

Use IQ-TREE 3 `--scfl` on the same fixed species tree and the retained
alignment directory.

Inputs:
- `results/species_tree/species_tree.wastral.tre`
- `results/loci/alignments/`

Outputs:
- `results/concordance/scfl.cf.stat`
- `results/concordance/scfl.cf.tree`
- `results/concordance/scfl.cf.branch`
- `results/concordance/scfl.log`

Why second:
- complements gCF with site-level signal
- reuses the same reference-tree workflow pattern
- requires more runtime than gCF, so it is better added after the gCF contract
  is stable

### Phase 3: ASTER Quartet Annotation (Optional)
Status: pending

Score the fixed `wastral` tree with `wastral -C -c ... -u 3` to emit ASTER-side
quartet summaries.

Inputs:
- `results/species_tree/species_tree.wastral.tre`
- `results/gene_trees/gene_trees.wastral.tre`

Outputs:
- annotated species tree
- `results/concordance/freqQuad.csv`

Why optional:
- useful for quartet-specific diagnostics
- output naming is less workflow-friendly than IQ-TREE `.cf.stat`
- better treated as an extra analytical layer after gCF and sCFL are in place

## Reporting Plan
- Phase 1 report update: list `gCF` artifacts and summarize branch-level gCF
  ranges or low-concordance branches.
- Phase 2 report update: add `sCFL` summary alongside `gCF`.
- Phase 3 report update: optionally list ASTER quartet outputs and selected
  `q1/q2/q3` summaries.

## Current Decision
Implement all three phases in order, keeping each phase as a separate commit:
`gCF`, then `sCFL`, then ASTER quartet annotation.
