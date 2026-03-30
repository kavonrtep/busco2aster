# Optimization Implementation Plan

## Goal
Reduce disk usage, reduce file count, and improve throughput without weakening
the current reproducible output contract.

This plan is based on the measured state in
[`docs/optimization_review.md`](/home/petr/PycharmProjects/get_phylo/docs/optimization_review.md):
- about `14G` under `results/`
- about `294,929` files
- storage dominated by BUSCO raw internals
- scheduler and filesystem overhead dominated by thousands of tiny locus jobs

## Current Status
- Completed: Phases 1 through 5 and 7
- Current target: Phase 6
- Phase 1 findings are recorded in
  [`docs/iqtree_directory_mode_evaluation.md`](/home/petr/PycharmProjects/get_phylo/docs/iqtree_directory_mode_evaluation.md)
- Current implemented optimization changes:
  - stable BUSCO sequence artifacts copied into `results/busco/*/busco_sequences/`
  - BUSCO raw scratch redirected to `work/busco/*/raw`
  - batched MAFFT execution
  - IQ-TREE 3 directory mode for gene-tree inference
  - cleanup tiers via `python3 -m scripts.cleanup_outputs`
- Remaining benchmark target:
  - BUSCO thread-scaling evaluation on the real dataset

## Design Principles
- Keep stable user-facing outputs in `results/`.
- Move large tool-internal restart artifacts to `work/` where possible.
- Prefer fewer, larger jobs over thousands of tiny jobs.
- Preserve restart points deliberately instead of keeping every tool byproduct.
- Benchmark before replacing a working stage.

## Upstream Basis for IQ-TREE Evaluation
The first optimization target should be IQ-TREE directory mode.

Relevant upstream documentation:
- IQ-TREE concordance-factor docs state that `iqtree3 -S ALN_DIR --prefix loci`
  can infer locus trees directly from a directory of alignments.
- The same docs state that `-p ALN_DIR` can use a directory of alignments for
  concatenation and `sCF` runs.

Sources:
- https://iqtree.github.io/doc/Concordance-Factor
- https://iqtree.github.io/doc/Command-Reference

## Phase 1: Evaluate IQ-TREE Directory Mode (Completed)
### Scope
Test whether IQ-TREE 3 directory mode can replace the current
one-job-per-locus gene-tree stage.

### Questions to answer
- Does `iqtree3 -S results/loci/alignments --prefix loci` produce all data we
  need for downstream aggregation?
- Does it preserve per-locus model selection and support values in a form
  compatible with `wastral -B`?
- Does it materially reduce runtime and file count?
- Does it still let us recover enough per-locus metadata for QC and reporting?

### Tasks
- create a benchmark fixture using a small subset of real retained loci
- run current per-locus workflow and directory mode on the same subset
- compare:
  - wall time
  - CPU utilization
  - file count
  - tree count
  - selected-model availability
  - support-label format
- inspect exactly which files IQ-TREE emits in directory mode
- decide whether per-locus manifests can be rebuilt from directory-mode outputs

### Required tests
- smoke run: `iqtree3 -S <alignment_dir>` on a small fixture
- regression test: tree count equals locus count
- regression test: all output trees contain the configured support format
- regression test: manifest reconstruction still yields `locus_id`,
  `selected_model`, and `support_values_present`

### Exit criteria
Choose one of:
- adopt IQ-TREE directory mode for Phase 4 implementation
- reject it and proceed with batched worker jobs instead

### Outcome
Adopt IQ-TREE directory mode as the leading implementation candidate. The
benchmark showed strong file-count reduction and acceptable output
compatibility, while the small slowdown in the tool-only benchmark is likely to
be offset by reducing workflow-level scheduler overhead.

## Phase 2: Define Retention Tiers and Cleanup Policy (Completed)
### Scope
Introduce an explicit output-retention model before changing heavy stages.

### Tasks
- define retention tiers such as:
  - `debug`
  - `resume_alignments`
  - `resume_gene_trees`
  - `final_report`
- classify each major subtree as:
  - stable deliverable
  - restart artifact
  - disposable tool scratch
- write a cleanup matrix for:
  - BUSCO raw internals
  - raw locus FASTAs
  - alignment logs
  - per-locus IQ-TREE byproducts
  - concordance logs

### Required tests
- doc review: each stage has a keep/drop policy
- dry-run design review: no proposed cleanup removes a currently required input

### Exit criteria
There is a documented keep/drop contract for every major output subtree.

### Outcome
Implemented retention tiers in `scripts.cleanup_outputs` with the supported
modes `debug`, `resume_alignments`, `resume_gene_trees`, and `final_report`.

## Phase 3: Slim BUSCO Stable Outputs (Completed)
### Scope
Reduce BUSCO disk usage while preserving downstream locus extraction.

### Tasks
- redesign stable BUSCO artifacts so downstream steps no longer depend on the
  full raw BUSCO tree
- preserve only the files actually required downstream:
  - stable `short_summary`
  - stable `full_table`
  - stable BUSCO `.faa`
  - stable BUSCO `.gff`
- stop preserving heavyweight internals such as:
  - `miniprot_output/ref.mpi`
  - miniprot stdout logs
  - `tmp/refseq_db.faa`
  - `hmmer_output`
  - translated intermediates if not needed downstream
- move disposable BUSCO scratch from `results/` to `work/`

### Required tests
- unit test: stable BUSCO path index resolves copied `.faa` / `.gff` files
- integration test: locus FASTA export succeeds after raw BUSCO scratch is
  removed
- regression test: BUSCO summary parsing yields identical QC tables before and
  after slimming

### Exit criteria
The workflow can rebuild retained-locus FASTAs without the original BUSCO raw
directory tree.

### Outcome
BUSCO raw scratch now lives under `work/busco/*/raw`, and stable BUSCO
sequence files are copied into `results/busco/*/busco_sequences/` so downstream
matrix and FASTA-export steps do not depend on the raw BUSCO tree.

## Phase 4: Reduce Gene-Tree Scheduler Overhead (Completed)
### Scope
Implement the gene-tree execution model chosen in Phase 1.

### Option A
Use IQ-TREE directory mode.

### Option B
Keep per-locus inference logic, but batch many loci into one worker job.

### Tasks
- implement the selected execution model behind the existing workflow contract
- preserve current aggregate outputs:
  - `results/gene_trees/gene_tree_manifest.tsv`
  - `results/gene_trees/gene_trees.raw.tre`
- decide which per-locus outputs remain stable and which become optional debug
  artifacts

### Required tests
- integration test: retained locus count equals gene-tree count
- integration test: aggregate tree file remains valid input to `wastral`
- benchmark: compare runtime and file count against current baseline

### Exit criteria
The gene-tree stage preserves biological outputs while materially reducing
runtime overhead and/or file count.

### Outcome
The workflow now uses IQ-TREE 3 directory mode for gene-tree inference and
reconstructs `gene_tree_manifest.tsv` plus `gene_trees.raw.tre` from the shared
directory-mode outputs.

## Phase 5: Batch Alignment Work (Completed)
### Scope
Reduce Snakemake overhead in the alignment stage.

### Tasks
- replace one-job-per-locus alignment execution with batched workers
- keep MAFFT single-threaded inside each locus unless benchmarking shows a
  clear benefit from larger per-job threading
- choose a batch size empirically, likely `100-500` loci per job

### Required tests
- integration test: alignment count equals retained locus count
- regression test: alignment filenames and manifest remain stable
- benchmark: compare wall time and average CPU utilization to current stage

### Exit criteria
The alignment stage spends most compute time inside MAFFT rather than in
scheduler overhead.

### Outcome
Alignment execution now runs in batches keyed by `alignment_batch_size`, with a
final synchronization step that removes stale alignments and logs.

## Phase 6: Rebalance BUSCO Parallelism
### Scope
Determine whether wider BUSCO execution is better than `1 sample x 4 threads`.

### Tasks
- benchmark BUSCO on the real dataset or a representative subset with:
  - `1 x 4`
  - `2 x 2`
  - `3 x 2` if the host has enough CPUs
- measure:
  - wall time per sample
  - total throughput
  - CPU utilization
  - memory pressure

### Required tests
- manual benchmark notes captured in `docs/`
- no regression in BUSCO summary outputs across thread settings

### Exit criteria
The repo has an evidence-based default BUSCO thread strategy.

## Phase 7: Final Cleanup Targets and Documentation (Completed)
### Scope
Expose cleanup and optimized execution as supported workflow behavior.

### Tasks
- add cleanup targets or helper commands for the approved retention tiers
- update README and contributor docs
- record before/after metrics for:
  - total bytes
  - total file count
  - alignment runtime
  - gene-tree runtime

### Required tests
- dry-run for cleanup targets
- manual verification that final-report mode still keeps all report inputs

### Exit criteria
Optimization behavior is documented, testable, and reproducible.

### Outcome
Cleanup modes are implemented and test-covered. Remaining documentation updates
are limited to routine user-facing command examples, not missing design work.

## Recommended Execution Order
1. Phase 1: evaluate IQ-TREE directory mode
2. Phase 2: define retention tiers
3. Phase 3: slim BUSCO outputs
4. Phase 4: optimize gene-tree execution
5. Phase 5: batch alignments
6. Phase 6: benchmark BUSCO parallelism
7. Phase 7: document and expose cleanup targets

## Recommendation
Start with Phase 1 exactly as requested. If IQ-TREE directory mode preserves the
current manifest/reporting contract and cuts scheduler overhead, it is the best
next change. If it does not, the fallback should be batched worker jobs rather
than the current one-job-per-locus model.
