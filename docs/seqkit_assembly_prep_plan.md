# SeqKit Assembly Preparation Plan

## Goal

Replace the current Python-based FASTA rewriting path with a `seqkit`-driven
assembly preparation stage that is faster, simpler, and better suited to very
large unwrapped genome FASTA files.

## Why SeqKit

`seqkit` is a good fit for this stage because the upstream documentation
confirms all of the required capabilities:

- installable via Conda/Bioconda
- `seqkit seq` supports configurable FASTA output wrapping with `-w`
- `seqkit seq` writes directly to `.gz` outputs with `-o`
- global `-j/--threads` support is available
- gzip I/O is handled internally and is explicitly optimized in SeqKit

This makes it a strong replacement for our current custom streaming writer.

## Recommended Implementation

### Phase 1: Validation Benchmark

Benchmark the current Python path against:

```bash
seqkit seq -j 4 -w 80 -o work/assemblies_prepared/<sample>.fa.gz <input>
```

Acceptance checks:

- `seqkit` output is valid FASTA and BUSCO `stats.sh` accepts it
- record count and total sequence length match the source file
- wall time and peak memory are not worse than the Python implementation

Benchmark on the problematic large inputs first:

- `V_formosa`
- `Cameor`
- `Champagne`

### Phase 2: Workflow Integration

Add a dedicated Conda env `workflow/envs/assembly_prep.yaml` with `seqkit`.

Update the workflow so:

- `prepare_assembly` runs under the `seqkit` env
- output remains `work/assemblies_prepared/{sample}.fa.gz`
- BUSCO continues to consume only prepared assemblies

Keep the current output location and naming stable to avoid downstream DAG
changes.

### Phase 3: QC Preservation

Retain a small Python wrapper around `seqkit` rather than replacing the rule
with raw shell only.

The wrapper should:

- call `seqkit seq`
- call `seqkit stats -T` on input and prepared output
- write a stable QC TSV for the sample

Recommendation:

- keep current fields like input/output paths and byte counts
- replace `max_input_line_length` with `seqkit`-derived sequence stats unless we
  decide that explicit line-length diagnostics are still required

### Phase 4: Container and CI

Update:

- `workflow/Snakefile_create_envs`
- `workflow/rules/_create_envs.smk`
- container build path

Add tests for:

- env creation
- dry-run DAG includes `prepare_assembly`
- prepared BUSCO input path is unchanged
- one real long-line FASTA smoke test rewrites successfully with `seqkit`

## Recommendation

Migrate the rewrite engine to `seqkit`, but keep a thin Python orchestration
layer for QC and stable metadata output. That gives us the performance benefits
of `seqkit` without losing reproducibility or report-friendly bookkeeping.

## Sources

- SeqKit usage: <https://bioinf.shenwei.me/seqkit/usage/>
- SeqKit download/install: <https://bioinf.shenwei.me/seqkit/download/>
- SeqKit project page: <https://bioinf.shenwei.me/seqkit/>
- Bioconda package: <https://anaconda.org/bioconda/seqkit>
