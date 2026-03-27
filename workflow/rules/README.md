# Rules

Snakemake rule files live here.

Current Phase 4 rule set:

- `manifest.smk`: normalized manifest validation and metadata generation
- `busco.smk`: BUSCO tool metadata capture, lineage verification, and per-sample BUSCO execution
- `busco_summary.smk`: BUSCO parsing and stable QC table generation

Later phases will add stage-specific rule files such as `alignment.smk` and
`species_tree.smk`.
