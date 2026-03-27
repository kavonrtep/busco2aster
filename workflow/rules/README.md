# Rules

Snakemake rule files live here.

Current Phase 3 rule set:

- `manifest.smk`: normalized manifest validation and metadata generation
- `busco.smk`: BUSCO tool metadata capture, lineage verification, and per-sample BUSCO execution

Later phases will add stage-specific rule files such as `alignment.smk` and
`species_tree.smk`.
