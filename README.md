# get_phylo

This repository will host a reproducible phylogenomics workflow built with
Snakemake.

Current status: Phase 1 manifest normalization and validation. The repository
now has a normalized sample sheet, a workflow config, and a working
`validate_manifest` rule that writes metadata outputs.

Useful commands:

```bash
snakemake -n
snakemake --cores 1
python3 -m scripts.normalize_manifest --input test_data/genome_set2.csv --output /tmp/samples.tsv
python3 -m unittest discover -s tests -v
```

Planning documents live in [`docs/`](/home/petr/PycharmProjects/get_phylo/docs).
