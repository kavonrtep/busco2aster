# Implementation Plan

## Scope
This document turns the design in [`problem_formulation.md`](/home/petr/PycharmProjects/get_phylo/docs/problem_formulation.md) into a concrete v1 implementation target. The workflow will ingest one genome assembly per taxon, run a shared BUSCO lineage in genome mode, retain complete single-copy loci, infer per-locus gene trees, collapse weak branches, and estimate an unrooted species tree with a pluggable coalescent backend.

The current test dataset is manifest-driven. The workflow should read [`test_data/genome_set2.csv`](/home/petr/PycharmProjects/get_phylo/test_data/genome_set2.csv) and use only the FASTA paths listed there as primary inputs. Extra per-sample files already present under `test_data/*/analysis/` should be ignored by v1.


> Review needed: confirm that `test_data/*/analysis/` is intentionally out of scope for the first implementation and should not be used as input metadata.

COMMENT - ignore analysis dir - the only input are genome files specifies in csv file.



## Proposed Repository Layout
The codebase should grow into this structure:

```text
docs/
  problem_formulation.md
  implementation.md
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

`config/samples.tsv` should be the normalized internal manifest, even if we keep a small adapter for the current CSV example. Required columns should be `sample_id`, `taxon_id`, and `assembly_fasta`. Optional columns can include `group`, `outgroup`, and `notes`.

> Review needed: keep compatibility with the current two-column CSV, or require an explicit normalized TSV from the start.

COMMENT: For purposes of development make new sample file with the required columns , do not include options columns for now, and update the workflow to read the new file. The old CSV can be kept as a reference or archived.


## Configuration
`config/config.yaml` should carry run-wide settings:

```yaml
busco_lineage: null
occupancy_threshold: 0.8
sequence_type: protein
alignment_tool: mafft
trimming_mode: none
collapse_support_threshold: 10
gene_tree_backend: iqtree2
species_tree_backend: astral3
threads:
  busco: 8
  iqtree: 4
```

Taxon names should be sanitized once at import time and stored in a mapping table so FASTA headers, IQ-TREE labels, and species-tree tips stay consistent.

> Review needed: for the current Solanum dataset, `solanales` is the obvious lineage guess, but the exact BUSCO lineage needs to be fixed before implementation.

COMMENT - verification of lineages must be part of pipeline, solanale lineage for test dataset is correct but exact string must be verified against BUSCO


 
 

## Workflow Stages
1. `validate_manifest`
   Read the sample sheet, resolve paths, ensure every assembly exists, and emit `results/metadata/samples.validated.tsv` plus `results/metadata/taxon_name_map.tsv`.
2. `run_busco`
   Run BUSCO genome mode once per sample into `results/busco/{sample_id}/`.
3. `summarize_busco`
   Parse `short_summary*`, `full_table.tsv`, and BUSCO sequence directories into `results/qc/busco_summary.tsv`.
4. `build_locus_matrix`
   Create `results/qc/locus_taxon_matrix.tsv` with per-locus, per-taxon status, sequence path, sequence length, translated length, and sequence QC flags.
5. `filter_loci`
   Apply occupancy and pathology filters and emit `results/qc/retained_loci.tsv`.
6. `export_locus_fastas`
   Write one FASTA per retained locus under `results/loci/raw_fastas/`.
7. `align_loci`
   Align each locus independently into `results/loci/alignments/`.
8. `trim_loci`
   Optionally create trimmed alignments in `results/loci/trimmed/` and record whether raw or trimmed alignments feed tree inference.
9. `infer_gene_trees`
   Run IQ-TREE per locus and collect trees, model files, and logs under `results/gene_trees/per_locus/`.
10. `collapse_gene_trees`
    Contract low-support branches and concatenate outputs into `results/gene_trees/gene_trees.collapsed.tre`.
11. `run_species_tree`
    Run the selected backend and emit `results/species_tree/species_tree.<backend>.tre`, an annotated tree, and the backend log.
12. `concordance_and_report`
    Compute gCF/sCF against the final species tree and render `results/report/report.md`.

## Tool Choices
The initial choices below keep v1 strict and inspectable:

- BUSCO v6 in genome mode.
- Protein sequences as the default downstream analysis layer.
- MAFFT for per-locus alignment.
- Optional trimming with default `none`.
- IQ-TREE 2 for per-locus ML trees with support enabled. - COMMENT  - there is version 3 - https://iqtree.github.io/  - probably not is conda - need to get binaries
- ASTRAL-III as the default species-tree backend behind a simple wrapper interface. - COMMENT - exact version of software must be verified against https://github.com/chaoszhang/ASTER - we may need actually ASTRAL-Pro to allow duplicated genes in future versions

Suggested IQ-TREE defaults are `-m MFP -B 1000 --bnni`, with outputs preserved per locus.

> Review needed: decide whether v1 should stay on IQ-TREE 2 explicitly or allow IQ-TREE 3 if it is easier to package on the target systems.

COMMENT - IQ-TREE 3 

> Review needed: decide whether the first backend should be `astral3` only, or whether we should scaffold `aster` in the interface now and leave it disabled.

COMMENT see above - check againt https://github.com/chaoszhang/ASTER 
 

## QC and Filtering Rules
The filter layer should be explicit and tabular, not hidden inside rule code. Each taxon-locus cell should carry at least:

- BUSCO status: `single`, `duplicated`, `fragmented`, `missing`
- amino-acid length
- internal-stop flag
- invalid-character flag
- length-z-score or another outlier metric
- keep/drop reason

The default locus rule for v1 should be:

- keep loci with occupancy `>= 0.8`
- keep only complete single-copy hits per retained taxon
- drop sequences flagged as pathological
- optionally drop loci with extreme length dispersion

> Review needed: define the exact cutoff for “extreme length dispersion.” My suggestion is to postpone hard filtering and only report dispersion in v1 unless the signal is clearly bad.

COMMENT - no filtering yet, just reporting of length dispersion and other QC flags, we can decide on hard cutoffs after we see the distribution on the test dataset


## Outputs
The workflow should produce these stable deliverables:

- `results/qc/busco_summary.tsv`
- `results/qc/locus_taxon_matrix.tsv`
- `results/qc/retained_loci.tsv`
- `results/loci/alignments/`
- `results/gene_trees/gene_trees.raw.tre`
- `results/gene_trees/gene_trees.collapsed.tre`
- `results/species_tree/species_tree.astral3.tre`
- `results/species_tree/species_tree.astral3.annotated.tre`
- `results/species_tree/species_tree.gcf_scf.tre`
- `results/report/report.md`

The report should summarize taxon completeness, retained-locus counts, occupancy distribution, duplicated-locus rates, support summaries, and discordance summaries.

## Testing Strategy
Tests should be split into three layers:

- parser tests for manifest normalization and BUSCO output parsing
- rule-level tests for locus filtering and branch-collapse behavior
- one smoke workflow run on a tiny fixture dataset once fixtures exist

The full `test_data/` directory is too large for routine CI-scale execution, so it should remain a manual integration dataset.

## Implementation Order
Build the repository in this order:

1. repository skeleton plus `Snakefile`
2. config parsing and manifest validation
3. BUSCO execution and parsing
4. locus matrix and filtering logic
5. alignment and IQ-TREE stages
6. species-tree backend wrapper
7. concordance metrics and final report

> Review needed: confirm whether you want the first coding pass to include environment packaging in `workflow/envs/` immediately, or whether we should start with the workflow skeleton and add locked environments after the rule graph is stable.

COMMENT - gradually build environments for rules, keep enviroments  simple -do not create env will all tools, snakemake enable graduall testing without reruning complete pipeline. So the workflow should be to define rule, env and test the rule, then move to next rule, this way we can build the workflow and environments in parallel without needing to rerun the whole workflow after each change.