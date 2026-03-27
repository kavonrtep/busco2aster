# Detailed Implementation Phases

## Purpose
This document turns the consolidated design into an execution plan. It is organized by implementation phase, with clear scope, outputs, and required tests for each phase.

Use this document as the build order for the workflow. [`implementation_consolidated.md`](/home/petr/PycharmProjects/get_phylo/docs/implementation_consolidated.md) remains the architectural decision record; this file is the delivery plan.

## Phase 0: Repository Bootstrap
### Scope
Create the repository skeleton needed to develop the workflow incrementally.

### Deliverables
- `Snakefile`
- `config/`
- `workflow/rules/`
- `workflow/envs/`
- `scripts/`
- `tests/`
- updated `.gitignore` entries for `results/` and `work/`

### Tasks
- create the top-level directory layout
- define the initial Snakemake entrypoint
- add a minimal `README` or header comments to explain how to run dry-runs
- choose the initial test runner and project layout for Python helpers

### Required Tests
- `snakemake -n` runs successfully against an empty placeholder target set
- repository tree matches the documented layout
- import smoke test for local helper modules if Python package structure is introduced

### Exit Criteria
The repository can execute a dry-run without workflow logic errors and has a stable place for rules, environments, scripts, and tests.

## Phase 1: Configuration and Manifest Normalization
### Scope
Introduce the normalized internal manifest and the run-wide configuration file.

### Deliverables
- `config/config.yaml`
- `config/samples.tsv`
- manifest normalization or conversion helper in `scripts/`
- `results/metadata/samples.validated.tsv`
- `results/metadata/taxon_name_map.tsv`

### Tasks
- define the canonical manifest schema: `sample_id`, `taxon_id`, `assembly_fasta`
- convert the current `test_data/genome_set2.csv` into the normalized TSV
- implement path resolution and taxon-name sanitization
- fail early on missing files, duplicate sample IDs, or invalid columns

### Required Tests
- unit test: valid manifest is parsed correctly
- unit test: missing required columns fail with a clear error
- unit test: duplicate sample IDs fail
- unit test: nonexistent assembly paths fail
- unit test: taxon sanitization is deterministic
- workflow test: `validate_manifest` target produces the expected metadata files

### Exit Criteria
All primary inputs are normalized into a stable internal format and the workflow can validate them before any heavy computation starts.

## Phase 2: BUSCO Environment and Lineage Verification
### Scope
Prepare BUSCO execution and verify that the configured lineage exists before launching per-sample jobs.

### Deliverables
- `workflow/envs/busco.yaml`
- `workflow/rules/busco.smk`
- lineage verification helper in `scripts/`

### Tasks
- define a dedicated BUSCO environment
- add a rule that verifies the BUSCO executable is available
- add a rule that checks the configured lineage string against the installed BUSCO datasets
- record the verified lineage in a small metadata artifact for downstream provenance

### Required Tests
- unit test: lineage verification passes for a mocked valid lineage name
- unit test: lineage verification fails for an invalid lineage name
- workflow test: `verify_busco_lineage` fails before `run_busco` when lineage is invalid
- environment test: BUSCO executable is resolvable inside the rule environment

### Exit Criteria
The workflow can prove that BUSCO is runnable and that the requested lineage identifier is valid before any BUSCO jobs are submitted.

## Phase 3: BUSCO Per-Sample Execution
### Scope
Run BUSCO genome mode once per sample and collect raw BUSCO outputs in a stable directory structure.

### Deliverables
- `results/busco/{sample_id}/`
- BUSCO logs and summaries per sample
- rule parameters for threads and output layout

### Tasks
- implement the per-sample BUSCO rule
- support compressed FASTA inputs from the manifest
- standardize BUSCO output paths so parsers do not depend on BUSCO’s volatile naming conventions
- capture the exact command line used for each sample

### Required Tests
- workflow dry-run: all BUSCO targets are created for the manifest samples
- unit test: expected output paths are derived correctly from a sample ID
- smoke test: BUSCO rule command renders correctly for one sample
- manual integration test: run BUSCO for one small sample and verify expected directories exist

### Exit Criteria
The workflow can launch BUSCO reproducibly for every validated sample and produce parseable outputs.

## Phase 4: BUSCO Parsing and Summary Tables
### Scope
Parse BUSCO outputs into stable tabular QC artifacts.

### Deliverables
- `results/qc/busco_summary.tsv`
- parser helpers in `scripts/`

### Tasks
- parse `short_summary*`
- parse `full_table.tsv`
- collect complete, duplicated, fragmented, and missing counts per sample
- expose paths to BUSCO-derived sequence files needed downstream

### Required Tests
- unit test: parse a BUSCO short summary fixture
- unit test: parse a BUSCO full table fixture
- unit test: malformed BUSCO output fails clearly
- integration test: parsed summary table contains one row per sample
- workflow test: `summarize_busco` depends on completed BUSCO outputs, not raw input files alone

### Exit Criteria
BUSCO results are represented in stable machine-readable tables, independent of BUSCO’s internal file naming.

## Phase 5: Locus Matrix Construction and Selection
### Scope
Build the core locus-by-taxon matrix and apply the strict v1 selection rules.

### Deliverables
- `results/qc/locus_taxon_matrix.tsv`
- `results/qc/retained_loci.tsv`
- selection and QC helper scripts

### Tasks
- assemble one row per locus and taxon cell
- annotate BUSCO status, protein length, stop-codon flags, invalid-character flags, and source paths
- compute occupancy per locus
- apply only the hard v1 filters: occupancy and complete single-copy status
- record QC warnings without excluding on dispersion yet

### Required Tests
- unit test: occupancy is computed correctly
- unit test: duplicated and fragmented BUSCO hits are excluded from retained loci
- unit test: stop-codon and invalid-character flags are propagated into the matrix
- unit test: loci below occupancy threshold are excluded
- integration test: retained loci table is consistent with the locus matrix
- regression test: no hard filtering is applied on length dispersion in v1

### Exit Criteria
The workflow produces a transparent locus matrix and a retained-loci decision table that fully explain why each locus passed or failed.

## Phase 6: Locus FASTA Export and Alignment
### Scope
Export one protein FASTA per retained locus and align each locus independently.

### Deliverables
- `workflow/envs/alignment.yaml`
- `workflow/rules/alignment.smk`
- `results/loci/raw_fastas/`
- `results/loci/alignments/`

### Tasks
- export retained sequences into locus FASTA files
- ensure FASTA headers use sanitized taxon labels
- define the MAFFT rule
- preserve a direct mapping from retained loci to alignment outputs

### Required Tests
- unit test: locus FASTA exporter includes only retained taxa
- unit test: FASTA headers use sanitized taxon IDs
- workflow dry-run: one alignment target is created per retained locus
- smoke test: run alignment on a tiny synthetic FASTA fixture
- integration test: alignment inputs match retained-loci outputs exactly

### Exit Criteria
Each retained locus has a deterministic aligned protein FASTA ready for tree inference.

## Phase 7: Gene-Tree Inference with IQ-TREE 3
### Scope
Infer one ML tree per locus and preserve branch lengths and support values needed by `wastral`.

### Deliverables
- `workflow/envs/iqtree.yaml` or binary wrapper strategy
- `workflow/rules/gene_trees.smk`
- `results/gene_trees/per_locus/`
- `results/gene_trees/gene_trees.raw.tre`

### Tasks
- standardize the IQ-TREE 3 command line
- decide the default support mode used for `wastral -B`
- collect tree files, reports, and selected models per locus
- concatenate per-locus trees into one aggregate input file

### Required Tests
- unit test: per-locus output paths are derived correctly
- unit test: aggregate tree list preserves all retained loci exactly once
- workflow dry-run: one gene-tree target is generated per alignment
- smoke test: run IQ-TREE on a tiny alignment fixture
- integration test: aggregate gene-tree file count matches retained locus count
- manual verification: support values are present in a format accepted by `wastral`

### Exit Criteria
The workflow produces a complete, support-aware gene-tree set suitable for species-tree inference.

## Phase 8: Species-Tree Inference with ASTER
### Scope
Run the default species-tree backend and keep the topology-only comparison path optional.

### Deliverables
- `workflow/envs/aster.yaml`
- `workflow/rules/species_tree.smk`
- `results/species_tree/species_tree.wastral.tre`
- `results/species_tree/species_tree.wastral.log`
- optional `results/species_tree/species_tree.astral4.tre`

### Tasks
- define the `wastral` execution rule
- keep `astral4` as an optional comparison rule, not the default
- standardize logging and provenance capture
- verify that the aggregate gene-tree format is accepted by ASTER

### Required Tests
- unit test: species-tree command rendering uses the configured backend
- workflow dry-run: default backend target resolves to `wastral`
- smoke test: run `wastral` on a tiny gene-tree fixture
- optional smoke test: run `astral4` comparison on the same fixture
- integration test: species-tree rule consumes the aggregate gene-tree output from Phase 7

### Exit Criteria
The workflow can infer a default species tree with ASTER and optionally compare it against a topology-only `astral4` run.

## Phase 9: Reporting and Concordance Summaries
### Scope
Produce the final audit trail and downstream QC report.

### Deliverables
- `results/report/report.md`
- optional concordance outputs if included in v1
- report generation helper scripts

### Tasks
- summarize BUSCO completeness by sample
- summarize candidate and retained locus counts
- report occupancy distribution and duplicated-locus rates
- summarize species-tree outputs and key runtime metadata
- add placeholders for concordance metrics if they are not implemented immediately

### Required Tests
- unit test: report generator renders all required sections
- integration test: report references the expected upstream artifacts
- regression test: missing optional concordance metrics do not break report generation
- manual review: report is readable and explains locus retention decisions

### Exit Criteria
The workflow emits a readable final report that explains what was run, what was retained, and where the main outputs are located.

## Phase 10: End-to-End Validation and Hardening
### Scope
Run the workflow on the real test dataset in stages and close gaps exposed by real execution.

### Deliverables
- validated end-to-end run notes in `docs/` or `results/report/`
- finalized per-rule environments
- updated defaults based on observed runtime behavior

### Tasks
- execute the workflow incrementally on the provided Solanum dataset
- record runtime bottlenecks and failure modes
- confirm that intermediate outputs are sufficient for debugging
- decide whether any QC warnings should become hard filters after inspecting real distributions

### Required Tests
- staged manual run: manifest validation through BUSCO on one sample
- staged manual run: locus matrix generation on all completed BUSCO outputs
- staged manual run: one small subset through species-tree inference
- full dry-run on the real manifest
- full integration run once environments and rules are stable

### Exit Criteria
The workflow is operational on the real dataset, the rule graph is stable, and any remaining issues are reduced to tuning rather than missing functionality.

## Cross-Phase Testing Rules
- Add unit tests as soon as a helper script or parser is introduced.
- Prefer synthetic fixtures for fast tests and keep the full `test_data/` dataset for manual integration.
- Every phase must add at least one workflow-level test, not just pure Python unit tests.
- No later phase should begin until the current phase has a dry-run path and at least one passing targeted test.

## Recommended Build Order
Implement phases in order from 0 through 10. The only acceptable overlap is:

- Phase 2 can start while finishing Phase 1 fixtures.
- Phase 6 environment work can start while finalizing Phase 5 logic.
- Phase 8 wrapper work can start once Phase 7 command shapes are stable.

Anything else will create too much churn for a workflow that still has no executable baseline.
