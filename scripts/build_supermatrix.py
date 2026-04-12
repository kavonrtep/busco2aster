"""Build a concatenated supermatrix and NEXUS partition file for the AU topology test."""

from __future__ import annotations

import csv
from pathlib import Path


def _read_fasta(path: Path) -> dict[str, str]:
    """Return an ordered dict of {taxon_label: sequence} from a FASTA file."""
    sequences: dict[str, str] = {}
    current_label: str | None = None
    parts: list[str] = []

    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if current_label is not None:
                sequences[current_label] = "".join(parts)
                parts = []
            current_label = line[1:].split()[0]
        else:
            parts.append(line.strip())

    if current_label is not None:
        sequences[current_label] = "".join(parts)

    return sequences


def _load_gene_tree_manifest(path: Path) -> dict[str, str]:
    """Return {locus_id: selected_model} from the gene-tree manifest TSV."""
    model_map: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            model_map[row["locus_id"]] = row["selected_model"]
    return model_map


def _load_retained_locus_ids(path: Path) -> list[str]:
    """Return ordered list of retained locus IDs from retained_loci.tsv."""
    locus_ids: list[str] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            locus_ids.append(row["locus_id"])
    return locus_ids


def build_supermatrix(
    *,
    alignment_dir: Path,
    retained_loci_path: Path,
    gene_tree_manifest_path: Path,
    output_phy: Path,
    output_nex: Path,
    au_test_model: str = "from_gene_trees",
    alignment_suffix: str = ".aln",
) -> None:
    """Concatenate per-locus alignments into a PHYLIP supermatrix and NEXUS partition file.

    Missing taxa in any locus are gap-padded to the full locus length.
    """
    locus_ids = _load_retained_locus_ids(retained_loci_path)
    model_map = _load_gene_tree_manifest(gene_tree_manifest_path)

    # Collect all taxon labels across all loci (preserving first-seen order).
    all_taxa: list[str] = []
    taxa_set: set[str] = set()
    locus_alignments: list[dict[str, str]] = []
    locus_lengths: list[int] = []

    for locus_id in locus_ids:
        aln_path = alignment_dir / f"{locus_id}{alignment_suffix}"
        if not aln_path.is_file():
            raise FileNotFoundError(f"Alignment not found for locus {locus_id!r}: {aln_path}")
        seqs = _read_fasta(aln_path)
        locus_alignments.append(seqs)
        lengths = {len(s) for s in seqs.values()}
        if len(lengths) != 1:
            raise ValueError(
                f"Inconsistent sequence lengths in alignment {aln_path}: {sorted(lengths)}"
            )
        locus_lengths.append(lengths.pop())
        for taxon in seqs:
            if taxon not in taxa_set:
                all_taxa.append(taxon)
                taxa_set.add(taxon)

    n_taxa = len(all_taxa)
    total_length = sum(locus_lengths)

    # Build concatenated sequences with gap-padding for missing taxa.
    concat_sequences: dict[str, list[str]] = {taxon: [] for taxon in all_taxa}
    for seqs, length in zip(locus_alignments, locus_lengths):
        gap_seq = "-" * length
        for taxon in all_taxa:
            concat_sequences[taxon].append(seqs.get(taxon, gap_seq))

    # Write relaxed PHYLIP.
    output_phy.parent.mkdir(parents=True, exist_ok=True)
    with output_phy.open("w", encoding="utf-8") as handle:
        handle.write(f" {n_taxa} {total_length}\n")
        for taxon in all_taxa:
            seq = "".join(concat_sequences[taxon])
            handle.write(f"{taxon}  {seq}\n")

    # Write NEXUS partition file.
    output_nex.parent.mkdir(parents=True, exist_ok=True)
    with output_nex.open("w", encoding="utf-8") as handle:
        handle.write("#nexus\n\nbegin sets;\n")
        position = 1
        charset_lines: list[str] = []
        model_lines: list[str] = []
        for locus_id, length in zip(locus_ids, locus_lengths):
            end = position + length - 1
            charset_lines.append(f"  charset {locus_id} = {position}-{end};")
            if au_test_model == "from_gene_trees":
                model = model_map.get(locus_id, "MFP")
            else:
                model = au_test_model
            model_lines.append(f"  {model}, {locus_id}")
            position += length
        handle.write("\n".join(charset_lines))
        handle.write("\nend;\n\nbegin models;\n")
        handle.write("\n".join(model_lines))
        handle.write("\nend;\n")
