# Storage and Parallelization Review

## Scope
This note documents the current state of the finished `busco2aster` run and
suggests optimizations. It is analysis only. No workflow behavior is changed
here.

## Current Snapshot
- Total output footprint: about `14G`
- Total file count under `results/`: `294,929`
- Main file-count hotspots:
  - `results/busco/*`: about `198k` files
  - `results/gene_trees/per_locus`: `74,347` files
  - `results/loci/raw_fastas`, `results/loci/alignments`, `results/loci/logs`: about `23.5k` files combined

The largest subtrees are:
- `results/busco/*`: about `13.0G` total
- `results/gene_trees/per_locus`: `586M`
- `results/concordance/scfl.log`: `463M`
- `results/loci/logs`: `31M`

## Storage Findings
### 1. BUSCO raw outputs dominate disk use
The main storage problem is not the stable BUSCO outputs in
`results/busco/{sample}/`, but the preserved raw BUSCO run trees under
`results/busco/{sample}/raw/`.

The biggest components across 6 samples are:
- `miniprot_output/ref.mpi`: `9653.4M`
- `logs/miniprot_align_*_out.log`: `2024.6M`
- `tmp/refseq_db.faa`: `455.4M`
- `hmmer_output`: `104.6M`
- `translated_proteins`: `66.6M`
- `busco_sequences`: `154.9M`

This means most BUSCO disk use comes from heavyweight internal byproducts, not
from the actual per-BUSCO `.faa` / `.gff` files used downstream.

### 2. Gene-tree directories dominate file count after BUSCO
Each retained locus gets its own IQ-TREE directory. Most downstream steps only
need:
- `results/gene_trees/gene_tree_manifest.tsv`
- `results/gene_trees/gene_trees.raw.tre`
- concordance summaries

Current per-locus IQ-TREE storage is mostly:
- `.iqtree`: `163.9M`
- `.log`: `148.2M`
- `.model.gz`: `30.3M`
- `.ckp.gz`: `7.9M`
- `.treefile`: `2.0M`

This is not huge in bytes, but it creates `74k` files.

### 3. Several outputs are useful only for restart, not for final reporting
After later stages complete, these are optional:
- `results/loci/raw_fastas/`
- `results/loci/logs/`
- large parts of `results/gene_trees/per_locus/`
- BUSCO raw internals
- `results/concordance/scfl.log`

## Cleanup Opportunities
### Recommended retention tiers
#### Tier A: Full debug state
Keep everything. Best for method development.

#### Tier B: Resume from locus export / alignment
Keep:
- BUSCO summary tables
- BUSCO `full_table.tsv` and `short_summary`
- sequence files actually needed to rebuild locus FASTAs
- all QC tables
- alignments
- aggregate gene-tree outputs

Drop:
- BUSCO `ref.mpi`
- BUSCO `tmp/refseq_db.faa`
- BUSCO miniprot stdout logs
- BUSCO `hmmer_output`
- BUSCO `translated_proteins`

Estimated savings: about `12.3G`.

#### Tier C: Resume from gene-tree / concordance stage
Keep:
- QC tables
- `results/loci/alignments/`
- `results/gene_trees/gene_tree_manifest.tsv`
- `results/gene_trees/gene_trees.raw.tre`
- species-tree and concordance outputs
- reports

Additionally drop:
- `results/loci/raw_fastas/`
- `results/loci/logs/`
- most of `results/gene_trees/per_locus/`

Estimated additional savings: about `650M` and about `98k` files.

#### Tier D: Final report bundle
Keep only:
- metadata and QC tables
- aggregate gene trees
- species tree
- concordance summary tables
- final reports

This would minimize both file count and storage, but removes convenient restart
points.

### Structural recommendation
The best long-term cleanup design is to separate:
- `results/`: stable user-facing outputs
- `work/`: large tool-internal scratch and restart artifacts

Right now BUSCO raw trees live under `results/`, which makes the final results
look much larger than the actual deliverables.

### Important constraint
Current `paths.tsv` records point into the BUSCO raw tree. That means BUSCO raw
cleanup is only safe:
- after `results/loci/raw_fastas.complete`, or
- after introducing a slimmer stable BUSCO artifact layout

## Parallelization Findings
### 1. BUSCO
Current defaults are:
- `threads.busco = 4`
- typical runs use `--cores 4`

That means only one BUSCO sample runs at a time. On a larger machine this is
conservative, but not necessarily optimal. The main optimization question is
whether `1 x 4 threads` is actually better than `2 x 2 threads` or `3 x 2
threads` for BUSCO genome mode on this dataset.

Recommendation:
- benchmark BUSCO with `2` versus `4` threads per sample
- if scaling is weak beyond `2`, run more samples concurrently instead

### 2. Alignment stage
The current design is already better than the earlier version because FASTA
export is batched and MAFFT runs with `1` thread per locus. That is the right
per-job threading choice for alignments with only `5-6` sequences.

The remaining inefficiency is Snakemake job overhead across `7836` tiny jobs.

Recommendation:
- batch many loci per Snakemake job, for example `100-500` loci per worker
- keep MAFFT itself single-threaded inside each batch

The goal is fewer Snakemake jobs, not more threads per alignment.

### 3. Gene-tree stage
This is the strongest parallelization target.

Current behavior:
- `7836` single-thread IQ-TREE jobs
- one directory per locus
- large scheduler and filesystem overhead

Better options to evaluate:
- chunked worker jobs that process many loci sequentially
- IQ-TREE directory mode for multi-locus input, if it preserves the required
  per-locus outputs and support values

Increasing threads per locus is unlikely to help much, because each locus is
small. Reducing job-launch overhead is the more promising direction.

### 4. Concordance
`gCF` and `wastral` are already cheap. `sCFL` is one heavier whole-dataset
command, which is a reasonable design. The main issue there is log volume, not
parallelism.

## Recommended Optimization Order
1. Introduce retention modes and cleanup targets, starting with BUSCO raw
   internals.
2. Move large tool-internal outputs from `results/` to `work/`.
3. Batch the alignment stage.
4. Batch the gene-tree stage or evaluate IQ-TREE directory mode.
5. Benchmark BUSCO thread scaling and choose between wider versus deeper
   parallelism.
6. Reduce or compress oversized logs, especially BUSCO miniprot logs and
   `results/concordance/scfl.log`.

## Practical First Targets
If the goal is maximum gain for minimal engineering effort, the first two tasks
should be:
- BUSCO cleanup policy
- gene-tree batching

Those two changes should remove most of the current storage burden and most of
the avoidable scheduler overhead.
