# IQ-TREE Directory Mode Evaluation

## Scope
This document records Phase 1 of the optimization plan:
evaluate whether IQ-TREE 3 directory mode is a viable replacement for the
current one-job-per-locus gene-tree stage.

## Environment
- IQ-TREE executable:
  `work/tools/iqtree3/current/bin/iqtree3`
- Installed version:
  `IQ-TREE version 3.1.0 for Linux x86 64-bit built Mar 12 2026`
- Benchmark root:
  `work/benchmarks/iqtree_dir_mode_jsxlgB`

## Benchmark Design
- subset size: `200` retained alignments
- sequence type: amino acid
- model mode: `MFP`
- support mode: `--abayes`
- seed: `20260327`
- CPU budget: `4`

Compared modes:
1. current execution model
   - one IQ-TREE run per locus
   - `4` concurrent workers
   - `-T 1` per locus
2. directory mode
   - one `iqtree3 -S <alignment_dir>` run
   - `-T 4`

Important caveat:
- the current-mode benchmark used direct IQ-TREE invocation, not Snakemake
  scheduling
- therefore it underestimates the real overhead of the current workflow, where
  `7836` locus jobs must also be scheduled and tracked

## Result Summary
### Timing
- current mode wall time: `15.81 s`
- directory mode wall time: `17.743 s`
- relative speed: directory mode was about `11%` slower in this tool-only
  benchmark

### Output volume
- current mode files: `1503`
- directory mode files: `9`
- file-count reduction: about `167x`

- current mode size: `9.44 MB`
- directory mode size: `1.376 MB`
- size reduction: about `6.9x`

### Output compatibility
- current mode trees: `200`
- directory mode trees: `200`
- per-locus model match rate: `199 / 200` = `99.5%`
- current mode trees with aBayes labels detected: `189 / 200`
- directory mode trees with aBayes labels detected: `200 / 200`

The single model mismatch was:
- `15853at4069`: current mode `FLU`, directory mode `FLU+I`

That locus has only `2` parsimony-informative sites in the current report, so
this looks like a weak-signal edge case rather than a structural incompatibility.

## Output Shape
Directory mode writes a compact shared output set:
- `loci.treefile`
- `loci.iqtree`
- `loci.best_model.nex`
- `loci.best_scheme`
- `loci.best_scheme.nex`
- `loci.model.gz`
- `loci.ckp.gz`
- `loci.log`
- `loci.parstree`

Key observations:
- `loci.treefile` contains one tree per input alignment, in alignment order
- `loci.iqtree` contains a single `Best-fit model according to BIC:` line with
  `model:locus` pairs
- `loci.best_model.nex` also records per-locus model assignments

This means the current manifest can likely be reconstructed without per-locus
IQ-TREE directories.

## Interpretation
### What directory mode does well
- drastically reduces file count
- drastically reduces output sprawl
- preserves one tree per locus
- preserves per-locus model information
- produced aBayes-labeled trees for all benchmark loci

### What it does not prove yet
- that full-dataset wall time will improve
- that full-dataset memory use will remain acceptable
- that users will accept losing per-locus IQ-TREE logs as first-class stable
  outputs

## Recommendation
Phase 1 outcome: **accept IQ-TREE directory mode as the leading implementation
candidate**.

Reasoning:
- compatibility with downstream aggregation looks good
- model reconstruction looks feasible
- support labels are preserved
- file-count reduction is substantial
- even though the tool-only benchmark was slightly slower, the real workflow is
  currently penalized by thousands of Snakemake jobs and per-locus directories

## Recommended Next Step
Proceed with the optimization plan under this assumption:
- keep directory mode as the target for gene-tree optimization
- do not switch defaults blindly
- first implement retention tiers and BUSCO slimming, then prototype directory
  mode in the real workflow behind a controlled benchmark path

## Source Basis
- IQ-TREE command help from the installed binary shows:
  - `-s DIR` for a directory of alignments
  - `-S FILE|DIR` for separate tree inference
- upstream IQ-TREE docs:
  - https://iqtree.github.io/doc/Concordance-Factor
  - https://iqtree.github.io/doc/Command-Reference
