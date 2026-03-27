# Rules

Snakemake rule files live here.

Current Phase 8 rule set:

- `manifest.smk`: normalized manifest validation and metadata generation
- `busco.smk`: BUSCO tool metadata capture, lineage verification, and per-sample BUSCO execution
- `busco_summary.smk`: BUSCO parsing and stable QC table generation
- `locus_matrix.smk`: long-form locus matrix construction and retained-locus selection
- `alignment.smk`: batched retained-locus FASTA export plus per-locus MAFFT alignment with one thread per locus
- `gene_trees.smk`: per-locus IQ-TREE 3 inference plus aggregate gene-tree exports
- `species_tree.smk`: ASTER integration, including `wastral` input normalization and the default species-tree rule
