# get_phylo

This repository will host a reproducible phylogenomics workflow built with
Snakemake.

Current status: Phase 6 retained-locus FASTA export and alignment. The repository now has a normalized
sample sheet, a workflow config, BUSCO tool preflight, per-sample BUSCO
execution, stable BUSCO QC tables, locus-level retention decisions, retained
protein FASTAs, and per-locus MAFFT alignments.

Useful commands:

```bash
snakemake -n
snakemake --cores 1 results/metadata/samples.validated.tsv results/metadata/taxon_name_map.tsv
snakemake --use-conda --cores 1 results/metadata/busco_lineage_verified.tsv
snakemake -n -p --cores 4 results/busco/solanum_chilense/run.complete
snakemake --cores 4 results/qc/busco_summary.tsv
snakemake --cores 4 results/qc/retained_loci.tsv
snakemake --cores 4 results/loci/alignments/35at4069.aln.faa
python3 -m scripts.normalize_manifest --input test_data/genome_set2.csv --output /tmp/samples.tsv
python3 -m unittest discover -s tests -v
```

Planning documents live in [`docs/`](/home/petr/PycharmProjects/get_phylo/docs).
