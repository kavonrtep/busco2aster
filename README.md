# get_phylo

This repository will host a reproducible phylogenomics workflow built with
Snakemake.

Current status: Phase 5 locus matrix construction and locus selection. The repository now has a normalized
sample sheet, a workflow config, BUSCO tool preflight, per-sample BUSCO
execution, stable BUSCO QC tables, and locus-level retention decisions.

Useful commands:

```bash
snakemake -n
snakemake --cores 1 results/metadata/samples.validated.tsv results/metadata/taxon_name_map.tsv
snakemake --use-conda --cores 1 results/metadata/busco_lineage_verified.tsv
snakemake -n -p --cores 4 results/busco/solanum_chilense/run.complete
snakemake --cores 4 results/qc/busco_summary.tsv
snakemake --cores 4 results/qc/retained_loci.tsv
python3 -m scripts.normalize_manifest --input test_data/genome_set2.csv --output /tmp/samples.tsv
python3 -m unittest discover -s tests -v
```

Planning documents live in [`docs/`](/home/petr/PycharmProjects/get_phylo/docs).
