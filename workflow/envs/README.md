# Environments

Place small, rule-specific environment definitions here.

The workflow will prefer one focused environment per rule or rule group instead
of a single large shared environment.

During container builds, `workflow/Snakefile_create_envs` and
`workflow/rules/_create_envs.smk` are used to pre-create these environments into
the image-level Conda prefix without invoking the main workflow DAG.
