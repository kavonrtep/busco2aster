# get_phylo

This repository will host a reproducible phylogenomics workflow built with
Snakemake.

Current status: Phase 8 species-tree inference is implemented. The repository now has a normalized
sample sheet, BUSCO tool preflight and per-sample execution, stable BUSCO QC
tables, locus-level retention decisions, batched retained-protein FASTA export,
per-locus MAFFT alignments, IQ-TREE 3 gene-tree inference, and ASTER-based
species-tree rules with `wastral` as the default backend.

Useful commands:

```bash
snakemake -n
snakemake --cores 1 results/metadata/samples.validated.tsv results/metadata/taxon_name_map.tsv
snakemake --use-conda --cores 1 results/metadata/busco_lineage_verified.tsv
snakemake -n -p --cores 4 results/busco/solanum_chilense/run.complete
snakemake --cores 4 results/qc/busco_summary.tsv
snakemake --cores 4 results/qc/retained_loci.tsv
snakemake --cores 4 results/loci/raw_fastas_manifest.tsv
snakemake --cores 4 results/loci/alignments/35at4069.aln.faa
snakemake --cores 4 results/loci/alignments.complete
python3 -m scripts.install_iqtree3
python3 -m scripts.install_aster
snakemake --cores 4 results/gene_trees/gene_trees.complete
snakemake --cores 4 results/species_tree/species_tree.complete
python3 -m scripts.normalize_manifest --input test_data/genome_set2.csv --output /tmp/samples.tsv
python3 -m unittest discover -s tests -v
```

Planning documents live in [`docs/`](/home/petr/PycharmProjects/get_phylo/docs).
