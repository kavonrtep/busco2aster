# Repository Guidelines

## Project Scope & Structure
This repository is for a reproducible phylogenomics workflow, not for the auxiliary tooling in [`hermit/`](/home/petr/PycharmProjects/get_phylo/hermit). Treat [`docs/problem_formulation.md`](/home/petr/PycharmProjects/get_phylo/docs/problem_formulation.md) as the biological design source of truth, [`docs/implementation_consolidated.md`](/home/petr/PycharmProjects/get_phylo/docs/implementation_consolidated.md) as the architectural plan, [`docs/implementation_phases.md`](/home/petr/PycharmProjects/get_phylo/docs/implementation_phases.md) as the execution plan, [`docs/implementation.md`](/home/petr/PycharmProjects/get_phylo/docs/implementation.md) as comment history, and [`test_data/genome_set2.csv`](/home/petr/PycharmProjects/get_phylo/test_data/genome_set2.csv) as the current sample manifest example. Do not modify `hermit/` unless explicitly asked.

The implementation target is a strict v1 pipeline: shared BUSCO lineage, complete single-copy loci, per-locus IQ-TREE 3 gene trees with support, and an unrooted coalescent-aware species tree from the ASTER toolkit using `wastral` by default. When code is added, keep workflow logic in `workflow/` or a top-level `Snakefile`, helper Python in `scripts/`, environment definitions in `envs/`, and generated outputs out of version control.

## Build, Test, and Development Commands
There is no runnable pipeline in the repository yet, so do not document placeholder commands as if they already work. Current useful commands are:

```bash
git status
sed -n '1,220p' docs/problem_formulation.md
sed -n '1,260p' docs/implementation_consolidated.md
sed -n '1,320p' docs/implementation_phases.md
column -t -s $'\t' test_data/genome_set2.csv
```

Once the workflow exists, prefer Snakemake for orchestration and support a dry-run first, for example `snakemake -n` before any full execution.

## Coding Style & Naming Conventions
Prefer Snakemake plus small Python helpers over monolithic shell scripts. Use `snake_case` for rule names, scripts, and config keys. Keep each rule responsible for one biological stage such as `run_busco`, `build_locus_matrix`, `infer_gene_trees`, or `run_species_tree`. Write explicit filenames and QC tables; avoid hidden heuristics.

## Testing Guidelines
Large genome files should not be committed. Use lightweight fixtures and metadata in [`test_data/`](/home/petr/PycharmProjects/get_phylo/test_data) for parser, filtering, and reporting tests. For workflow changes, require at minimum a Snakemake dry-run, targeted unit tests for filtering logic, and one small end-to-end smoke test when fixtures exist.

## Commit & Pull Request Guidelines
The root repository was initialized on March 27, 2026 and has no commit history yet. Start with short imperative subjects such as `Add Snakemake skeleton` or `Implement BUSCO locus matrix parser`. Pull requests should state which pipeline stage changed, how it was validated, and whether outputs or QC criteria changed.
