# Rules

Snakemake rule files live here.

Current implemented rule set:

- `manifest.smk`: normalized manifest validation and metadata generation
- `busco.smk`: BUSCO tool metadata capture, lineage verification, and per-sample BUSCO execution
- `busco_summary.smk`: BUSCO parsing and stable QC table generation
- `locus_matrix.smk`: long-form locus matrix construction and retained-locus selection
- `alignment.smk`: batched retained-locus FASTA export plus per-locus MAFFT alignment with one thread per locus
- `gene_trees.smk`: per-locus IQ-TREE 3 inference plus aggregate gene-tree exports
- `species_tree.smk`: ASTER integration, including `wastral` input normalization and the default species-tree rule
- `concordance.smk`: IQ-TREE `gCF`/`sCFL` scoring and ASTER quartet scoring on the final species tree
- `report.smk`: final Markdown report generation from QC tables and species-tree outputs
- `visual_report.smk`: report-data bundling plus Quarto HTML rendering
- `_create_envs.smk`: dummy rules used only during container build to pre-create Conda environments with the same relative env YAML paths as the real rules
