# Apptainer / Singularity Containerization Plan

## Purpose

Containerize `busco2aster` so the full workflow can run reproducibly with
Apptainer or Singularity, while preserving the current native Snakemake mode for
development.

This plan is based on the existing container pattern used in the local
`orthoTE` repository, especially:

- `/home/petr/PycharmProjects/orthoTE/orthoTE.def`
- `/home/petr/PycharmProjects/orthoTE/run_pipeline.py`
- `/home/petr/PycharmProjects/orthoTE/workflow/Snakefile_create_envs`
- `/home/petr/PycharmProjects/orthoTE/config/config_comparative_container.yaml`
- `/home/petr/PycharmProjects/orthoTE/.github/workflows/build-sif.yaml`

## Reference Pattern Extracted From `orthoTE`

The relevant design choices in `orthoTE` are:

1. One monolithic image definition at the repository root.
2. Pipeline code copied into `/opt/pipeline` inside the image.
3. Snakemake installed in the image base environment.
4. Per-rule conda environments pre-created during image build.
5. A container-aware wrapper script that:
   - validates input paths
   - prints suggested bind mounts
   - runs Snakemake with container-internal paths
   - keeps outputs on the host via bind-mounted writable directories
6. A GitHub Actions workflow that builds a `.sif` image on tags and uploads it
   to the release.

This is a good fit for `busco2aster`. The main difference is that
`busco2aster` currently mixes per-rule conda environments with external tool
install scripts for IQ-TREE 3 and ASTER, so the container plan must standardize
that runtime.

## Recommended Design For `busco2aster`

Adopt the same high-level model as `orthoTE`:

- one Apptainer definition file at repo root: `busco2aster.def`
- one runtime wrapper at repo root: `run_pipeline.py`
- image-internal pipeline root: `/opt/pipeline`
- image-internal shared conda env prefix: `/opt/conda/envs`
- image-internal tool installs:
  - `/opt/tools/iqtree3/current/bin/iqtree3`
  - `/opt/tools/aster/current/bin/wastral`
  - `/opt/tools/aster/current/bin/astral4`

Key runtime rule:

- the container must never write into `/opt/pipeline`
- all `results/`, `work/`, `.snakemake/`, and BUSCO download caches must live in
  the bind-mounted host working directory

## Explicit Design Decisions

### 1. Keep a single image, not per-rule containers

For this project, a single SIF is the simplest operational model. The workflow
already has only two conda env YAMLs and three external executables. A
single-image model avoids mixing Snakemake `container:` directives with the
current rule graph.

### 2. Pre-create conda envs inside the image

Follow the `orthoTE` pattern and add a small helper Snakefile solely for
environment creation. This avoids running `--conda-create-envs-only` against the
main workflow, which already contains checkpoints and dynamic inputs.

Target envs to pre-create now:

- `workflow/envs/busco.yaml`
- `workflow/envs/alignment.yaml`

### 3. Install IQ-TREE 3 and ASTER during image build

Do not rely on `work/tools/...` inside the host workspace when running the
container. The image should install these once during build and the wrapper
should inject container-internal executable paths into Snakemake.

### 4. Do not bake BUSCO lineage datasets into the image in v1

Keep BUSCO downloads external under the host working directory, using the
existing `work/busco_downloads` path. This keeps the image smaller and avoids
locking the image to one dataset cache snapshot.

### 5. Prefer wrapper-injected path overrides over a second full config

Unlike `orthoTE`, `busco2aster` currently stores executable paths directly in
`config/config.yaml`. To avoid config drift, the container wrapper should pass
Snakemake config overrides for:

- `iqtree_executable`
- `wastral_executable`
- `astral4_executable`

That keeps `config/config.yaml` as the single source of truth for native runs.

## Files To Add

- `busco2aster.def`
- `run_pipeline.py`
- `workflow/Snakefile_create_envs`
- `workflow/rules/_create_envs.smk`
- `.github/workflows/build-sif.yaml`
- optional: `docs/container_usage.md`

## Files To Update

- `Snakefile`
- `README.md`
- `AGENTS.md`
- possibly `config/config.yaml` if any container-friendly defaults are needed

## Implementation Phases

### Phase 1: Runtime Contract And Path Audit

Tasks:

- define the container-internal paths for Snakemake, conda envs, IQ-TREE 3, and ASTER
- verify which paths must remain writable on the host
- confirm whether the current `configfile: "config/config.yaml"` behavior requires
  a wrapper-created `config/` symlink in the working directory, as `orthoTE` does
- export `PYTHONPATH=/opt/pipeline` from the container wrapper if needed for
  top-level `from scripts...` imports

Required tests:

- `apptainer exec <image> snakemake --version`
- `apptainer exec <image> bash -lc 'iqtree3 --version'`
- `apptainer exec <image> bash -lc 'wastral -h >/dev/null'`

### Phase 2: Image Definition And Build-Time Provisioning

Tasks:

- create `busco2aster.def`
- copy the pipeline code into `/opt/pipeline`
- install Snakemake in the image base environment
- run `python3 -m scripts.install_iqtree3` during build against `/opt/tools/iqtree3`
- run `python3 -m scripts.install_aster` during build against `/opt/tools/aster`
- add helper workflow to pre-create `busco.yaml` and `alignment.yaml` envs into
  `/opt/conda/envs`
- clean build-only packages and conda caches at the end of `%post`

Required tests:

- local `apptainer build busco2aster.sif busco2aster.def`
- smoke check that the expected binaries exist inside the image

### Phase 3: Wrapper Script

Tasks:

- add a top-level `run_pipeline.py`
- accept at minimum:
  - `-c/--config`
  - `-t/--threads`
  - `--target`
  - `-S/--snakemake-args`
  - optional `--directory`
- validate that input files in the config are readable from the container
- print suggested `-B` bind mounts on failure
- invoke Snakemake with:
  - internal Snakefile path
  - `--use-conda`
  - `--conda-prefix /opt/conda/envs`
  - container-internal executable overrides

Required tests:

- `apptainer run ... busco2aster.sif --help`
- dry-run to `results/metadata/samples.validated.tsv`
- dry-run to `results/report/report.md`

### Phase 4: CI Build And Release

Tasks:

- add `.github/workflows/build-sif.yaml`
- trigger on tags `v*` and optionally `workflow_call`
- install Apptainer with `eWaterCycle/setup-apptainer`
- build `busco2aster.sif`
- run at least one smoke command inside the built image
- upload the `.sif` artifact to the GitHub release

Recommended smoke test in CI:

- `apptainer exec busco2aster.sif snakemake --version`
- `apptainer exec busco2aster.sif bash -lc 'iqtree3 --version && wastral -h >/dev/null'`

### Phase 5: User Documentation

Tasks:

- add README usage examples for:
  - native mode
  - `apptainer run`
  - `singularity run`
- document required bind mounts
- document that BUSCO download cache and all outputs must be on the host
- provide one example command for staged execution and one for the final report

Required tests:

- README commands match the wrapper CLI exactly
- container usage documentation references the current config fields and targets

## Minimal First Deliverable

The smallest acceptable v1 containerization milestone is:

1. buildable `busco2aster.def`
2. working `run_pipeline.py`
3. CI workflow that builds and uploads `busco2aster.sif`
4. successful containerized dry-run to `results/report/report.md`
5. successful containerized execution of at least:
   - manifest validation
   - BUSCO lineage verification

## Risks To Handle Early

- Snakemake path resolution for `configfile:` and Python helper imports may differ
  when the Snakefile lives in `/opt/pipeline` but execution happens in a
  bind-mounted host directory.
- BUSCO runtime may assume writable cache or download directories; these must not
  resolve inside the read-only image.
- ASTER build requirements may increase image build time enough that CI needs
  extra disk cleanup before `apptainer build`.
- The current repository has no GitHub workflow directory yet, so the CI layout
  must be introduced from scratch rather than adapted.

## Recommended Execution Order

Implement in this order:

1. helper env-creation workflow
2. `busco2aster.def`
3. `run_pipeline.py`
4. local image build and smoke tests
5. GitHub Actions build-and-release workflow
6. README container usage docs

Do not start with CI before the local image can run the workflow wrapper
successfully.
