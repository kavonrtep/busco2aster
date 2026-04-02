# DNA GFFRead Implementation Plan

## Goal

Add an opt-in `sequence_type: dna` mode while keeping `protein` as the default. The DNA path will reuse the existing BUSCO Miniprot workflow and derive CDS sequences from BUSCO Miniprot GFF files with `gffread`, rather than switching BUSCO to Metaeuk.

## Rationale

- BUSCO Miniprot already produces per-locus `.gff` and `.faa` files.
- The Miniprot GFFs contain usable `mRNA` and `CDS` features on both `+` and `-` strands.
- This avoids the large Metaeuk temporary I/O burden observed during testing.
- `gffread` is available through conda and can be packaged in both Snakemake envs and the container.

## Design

### 1. Keep BUSCO in Miniprot mode

- Do not change the current default BUSCO execution path.
- DNA mode will still run BUSCO with Miniprot.
- BUSCO summary, locus selection, and protein QC remain unchanged.

### 2. Add a DNA extraction stage after locus selection

- Input:
  - prepared assembly FASTA for each sample
  - BUSCO `.gff` files for retained loci
- Tool:
  - `gffread`
- Output:
  - per-locus CDS FASTA files for retained loci

Recommended execution shape:
- aggregate retained-locus GFF records per sample
- run one `gffread` job per sample, not one per locus
- parse sample-level CDS FASTA back into locus/sample records for downstream export

This avoids thousands of tiny `gffread` jobs.

### 3. Prepare FASTA for `gffread`

Current prepared assemblies are gzip-compressed. `gffread` works best with plain indexed FASTA.

Implementation options:
- preferred: extend assembly preparation to emit plain wrapped FASTA plus `.fai`
- include `samtools` in the extraction env to build `.fai`

BUSCO can keep using the current prepared assembly artifact unless we decide to unify on plain FASTA for both steps.

### 4. DNA-specific downstream behavior

- raw retained locus files: `.fna`
- alignments: nucleotide MAFFT
- IQ-TREE gene trees: `--seqtype DNA`
- sCF: DNA model default, e.g. `GTR+G4`
- gCF unchanged
- reports should label alignment units as `nt`

## Implementation Phases

### Phase 1. Environment and config

- add `gffread` env with `gffread` and `samtools`
- add config switch for `sequence_type: dna`
- keep `protein` default

Tests:
- env smoke test for `gffread --help`
- config dry-run for both `protein` and `dna`

### Phase 2. Sample-level CDS extraction

- build retained-locus sample GFF bundles
- run `gffread` per sample
- split extracted CDS back to locus/sample records

Tests:
- validate one `+` and one `-` strand BUSCO locus from test data
- compare translated extracted CDS against BUSCO `.faa`

### Phase 3. DNA export/alignment/tree path

- wire DNA exports, alignments, IQ-TREE, sCF, and reports

Tests:
- workflow dry-run for final DNA target
- targeted DNA smoke run on a small retained-locus subset

### Phase 4. Full test-data run

- execute the DNA pipeline on the provided test dataset
- confirm final species tree, concordance outputs, and reports are generated

Tests:
- end-to-end run to `results/report/report.html`
- compare retained-locus counts and taxon set against protein mode

## Open Decisions

- whether to keep BUSCO input as `.fa.gz` and generate a second plain indexed FASTA only for `gffread`
- whether sample-level `gffread` outputs should be kept in `work/` only or partially exposed in `results/`
- whether to validate every extracted CDS against BUSCO `.faa` or only in test/debug mode
