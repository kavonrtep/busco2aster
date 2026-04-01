# Config

This directory will hold normalized workflow inputs such as `samples.tsv` and
run-wide settings such as `config.yaml`.

Current contents:

- `config.yaml`: workflow-level settings and file locations
- `config.template.yaml`: minimal starting template for new runs
- `samples.tsv`: normalized internal manifest with `sample_id`, `taxon_id`,
  and `assembly_fasta`

Assembly inputs are normalized before BUSCO into wrapped gzipped FASTA files
under `assembly_prep_root` (default: `work/assemblies_prepared`). This stage is
implemented with `seqkit` in a dedicated Conda environment. The corresponding
`assembly_wrap_width` setting defaults to `80`.

When the workflow is launched through `run_pipeline.py`, `thread_policy: auto`
is the default. In that mode, the wrapper uses the requested global core budget
to derive per-rule thread settings automatically. Set `thread_policy: fixed` if
you want the `threads:` block to be used exactly as written.

Executable path keys are optional. If `iqtree_executable`,
`wastral_executable`, or `astral4_executable` are omitted, the workflow uses
`iqtree3`, `wastral`, and `astral4` from `PATH`. This is the expected mode for
container runs.
