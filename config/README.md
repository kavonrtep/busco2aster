# Config

This directory will hold normalized workflow inputs such as `samples.tsv` and
run-wide settings such as `config.yaml`.

Current contents:

- `config.yaml`: workflow-level settings and file locations
- `config.template.yaml`: minimal starting template for new runs
- `samples.tsv`: normalized internal manifest with `sample_id`, `taxon_id`,
  and `assembly_fasta`

Executable path keys are optional. If `iqtree_executable`,
`wastral_executable`, or `astral4_executable` are omitted, the workflow uses
`iqtree3`, `wastral`, and `astral4` from `PATH`. This is the expected mode for
container runs.
