"""Prepare report-ready tables for the visual HTML report."""

from __future__ import annotations

import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from .concordance import read_cf_stat_rows, read_freqquad_rows
from .manifest import write_tsv
from .sequence_mode import alignment_suffix, sequence_length_unit
from .tree_utils import (
    branch_label_map,
    canonical_branch_key,
    canonical_topology_key,
    fill_missing_branch_lengths,
    informative_branch_key,
    internal_branch_key_map,
    iter_nodes,
    leaf_labels,
    parse_newick,
    read_newick,
    relabel_tree_with_branch_ids,
    write_newick,
)


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts)))
            header = line[1:]
            parts = []
            continue
        parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def _parse_float(value: str) -> float | None:
    text = value.strip()
    if text in {"", "NA", "N/A", "nan", "NaN"}:
        return None
    return float(text)


def _parse_int(value: str) -> int | None:
    text = value.strip()
    if text in {"", "NA", "N/A"}:
        return None
    return int(float(text))


def _resolve_path(repo_root: Path, path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else repo_root / path


def _format_float(value: float | None, digits: int = 6) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def _extract_sample_qc(
    *,
    busco_summary_rows: list[dict[str, str]],
    locus_matrix_rows: list[dict[str, str]],
    retained_locus_ids: set[str],
) -> list[dict[str, str]]:
    sample_stats: dict[str, Counter] = defaultdict(Counter)
    retained_total = len(retained_locus_ids)

    for row in locus_matrix_rows:
        if row["locus_id"] not in retained_locus_ids:
            continue
        sample_id = row["sample_id"]
        sample_stats[sample_id]["retained_loci_total"] = retained_total
        if row["include_in_occupancy"].strip().lower() == "true":
            sample_stats[sample_id]["retained_present_loci"] += 1
        else:
            sample_stats[sample_id]["retained_missing_loci"] += 1
        if row["status"] == "Duplicated":
            sample_stats[sample_id]["retained_duplicated_loci"] += 1
        if row["status"] == "Fragmented":
            sample_stats[sample_id]["retained_fragmented_loci"] += 1
        if row["has_internal_stop_codon"].strip().lower() == "true":
            sample_stats[sample_id]["retained_internal_stop_loci"] += 1

    rows: list[dict[str, str]] = []
    for row in busco_summary_rows:
        sample_id = row["sample_id"]
        stats = sample_stats[sample_id]
        missing_fraction = (
            stats["retained_missing_loci"] / retained_total if retained_total else 0.0
        )
        rows.append(
            {
                "sample_id": sample_id,
                "taxon_id": row["taxon_id"],
                "sanitized_taxon_id": row["sanitized_taxon_id"],
                "complete_percent": row["complete_percent"],
                "single_copy_percent": row["single_copy_percent"],
                "multi_copy_percent": row["multi_copy_percent"],
                "fragmented_percent": row["fragmented_percent"],
                "missing_percent": row["missing_percent"],
                "internal_stop_codon_count": row["internal_stop_codon_count"],
                "retained_loci_total": str(retained_total),
                "retained_present_loci": str(stats["retained_present_loci"]),
                "retained_missing_loci": str(stats["retained_missing_loci"]),
                "retained_missing_fraction": _format_float(missing_fraction, digits=4),
                "retained_duplicated_loci": str(stats["retained_duplicated_loci"]),
                "retained_fragmented_loci": str(stats["retained_fragmented_loci"]),
                "retained_internal_stop_loci": str(stats["retained_internal_stop_loci"]),
            }
        )
    return rows


def _extract_alignment_summary(
    *,
    repo_root: Path,
    retained_rows: list[dict[str, str]],
    alignment_dir: Path,
    sequence_type: str,
) -> list[dict[str, str]]:
    summary_rows: list[dict[str, str]] = []
    for row in retained_rows:
        if row["decision"] != "retain":
            continue
        alignment_path = _resolve_path(
            repo_root,
            (alignment_dir / f"{row['locus_id']}.{alignment_suffix(sequence_type)}").as_posix(),
        )
        records = read_fasta_records(alignment_path)
        lengths = {len(sequence) for _, sequence in records}
        if len(lengths) != 1:
            raise ValueError(f"Alignment lengths are inconsistent in {alignment_path}.")
        alignment_length = lengths.pop() if lengths else 0
        gap_count = sum(sequence.count("-") for _, sequence in records)
        total_characters = alignment_length * len(records)
        gap_fraction = gap_count / total_characters if total_characters else 0.0
        summary_rows.append(
            {
                "locus_id": row["locus_id"],
                "sequence_count": str(len(records)),
                "alignment_length_sites": str(alignment_length),
                "gap_fraction": _format_float(gap_fraction, digits=6),
                "occupancy": row["occupancy"],
                "length_dispersion_observed": row["length_dispersion_observed"],
            }
        )
    return summary_rows


def _read_branch_stats(path: Path) -> dict[str, dict[str, str]]:
    rows = {}
    for row in read_cf_stat_rows(path):
        rows[row["ID"]] = row
    return rows


def _map_branch_ids(path: Path) -> dict[str, str]:
    root = read_newick(path)
    total_taxa = leaf_labels(root)
    mapped: dict[str, str] = {}
    for node in iter_nodes(root):
        branch_key = informative_branch_key(node, total_taxa)
        if branch_key is None or not node.label:
            continue
        mapped[branch_key] = node.label
    return mapped


def _split_key_parts(branch_key: str) -> tuple[str, str]:
    left, right = branch_key.split("||", 1)
    return left, right


def _quartet_branch_key(split_text: str) -> str:
    left_text, right_text = split_text.split("#", 1)
    left_taxa = set()
    right_taxa = set()
    for part in left_text.split("|"):
        left_taxa.update(taxon.strip() for taxon in part.strip("{}").split(",") if taxon.strip())
    for part in right_text.split("|"):
        right_taxa.update(taxon.strip() for taxon in part.strip("{}").split(",") if taxon.strip())
    return canonical_branch_key(left_taxa, right_taxa)


def _extract_branch_tables(
    *,
    species_tree_path: Path,
    output_tree_path: Path,
    gcf_stat_path: Path,
    gcf_branch_path: Path,
    scfl_stat_path: Path,
    scfl_branch_path: Path,
    quartet_freqquad_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    species_root = read_newick(species_tree_path)
    labeled_root = relabel_tree_with_branch_ids(species_root)
    fill_missing_branch_lengths(labeled_root, default="0")
    write_newick(output_tree_path, labeled_root)

    species_branch_nodes = internal_branch_key_map(species_root)
    branch_labels = branch_label_map(species_root)
    gcf_stats = _read_branch_stats(gcf_stat_path)
    scfl_stats = _read_branch_stats(scfl_stat_path)
    gcf_branch_ids = _map_branch_ids(gcf_branch_path)
    scfl_branch_ids = _map_branch_ids(scfl_branch_path)

    quartet_groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    quartet_branch_keys: dict[str, str] = {}
    for row in read_freqquad_rows(quartet_freqquad_path):
        quartet_groups[str(row["node_id"])].append(row)
        if str(row["topology_id"]) == "t1":
            quartet_branch_keys[str(row["node_id"])] = _quartet_branch_key(str(row["split"]))

    quartet_by_branch: dict[str, dict[str, object]] = {}
    for node_id, rows in quartet_groups.items():
        branch_key = quartet_branch_keys.get(node_id)
        if branch_key is None:
            continue
        row_by_topology = {str(row["topology_id"]): row for row in rows}
        total_weight = float(next(iter(rows))["total_weight"])
        quartet_by_branch[branch_key] = {
            "node_id": node_id,
            "local_posterior": float(row_by_topology["t1"]["local_posterior"]),
            "q1": float(row_by_topology["t1"]["weighted_support"]) / total_weight if total_weight > 0 else None,
            "q2": float(row_by_topology["t2"]["weighted_support"]) / total_weight if total_weight > 0 else None,
            "q3": float(row_by_topology["t3"]["weighted_support"]) / total_weight if total_weight > 0 else None,
            "f1": float(row_by_topology["t1"]["weighted_support"]),
            "f2": float(row_by_topology["t2"]["weighted_support"]),
            "f3": float(row_by_topology["t3"]["weighted_support"]),
            "total_weight": total_weight,
            "rows": rows,
        }

    branch_metric_rows: list[dict[str, str]] = []
    branch_alternative_rows: list[dict[str, str]] = []

    for branch_key in sorted(species_branch_nodes):
        node = species_branch_nodes[branch_key]
        branch_label = branch_labels[branch_key]
        split_left, split_right = _split_key_parts(branch_key)
        support_value = _parse_float(node.label or "")
        gcf_branch_id = gcf_branch_ids.get(branch_key)
        scfl_branch_id = scfl_branch_ids.get(branch_key)
        gcf_row = gcf_stats.get(gcf_branch_id or "")
        scfl_row = scfl_stats.get(scfl_branch_id or "")
        quartet_row = quartet_by_branch.get(branch_key, {})

        branch_metric_rows.append(
            {
                "label": branch_label,
                "branch_key": branch_key,
                "split_left": split_left,
                "split_right": split_right,
                "species_tree_support": _format_float(support_value, digits=6),
                "species_tree_branch_length": node.length or "",
                "gcf_branch_id": gcf_branch_id or "",
                "gcf": _format_float(_parse_float(gcf_row["gCF"]) if gcf_row else None, digits=4),
                "gdf1": _format_float(_parse_float(gcf_row["gDF1"]) if gcf_row else None, digits=4),
                "gdf2": _format_float(_parse_float(gcf_row["gDF2"]) if gcf_row else None, digits=4),
                "gdfp": _format_float(_parse_float(gcf_row["gDFP"]) if gcf_row else None, digits=4),
                "gn": str(_parse_int(gcf_row["gN"]) if gcf_row else "" or ""),
                "scfl_branch_id": scfl_branch_id or "",
                "scfl": _format_float(_parse_float(scfl_row["sCF"]) if scfl_row else None, digits=4),
                "sdf1": _format_float(_parse_float(scfl_row["sDF1"]) if scfl_row else None, digits=4),
                "sdf2": _format_float(_parse_float(scfl_row["sDF2"]) if scfl_row else None, digits=4),
                "sn": _format_float(_parse_float(scfl_row["sN"]) if scfl_row else None, digits=4),
                "aster_node_id": str(quartet_row.get("node_id", "")),
                "aster_local_posterior": _format_float(quartet_row.get("local_posterior"), digits=6),
                "aster_q1": _format_float(quartet_row.get("q1"), digits=6),
                "aster_q2": _format_float(quartet_row.get("q2"), digits=6),
                "aster_q3": _format_float(quartet_row.get("q3"), digits=6),
                "aster_f1": _format_float(quartet_row.get("f1"), digits=4),
                "aster_f2": _format_float(quartet_row.get("f2"), digits=4),
                "aster_f3": _format_float(quartet_row.get("f3"), digits=4),
                "aster_total_weight": _format_float(quartet_row.get("total_weight"), digits=4),
            }
        )

        if gcf_row:
            for key, value_key in (
                ("main", "gCF"),
                ("df1", "gDF1"),
                ("df2", "gDF2"),
                ("dfp", "gDFP"),
            ):
                value = _parse_float(gcf_row[value_key])
                branch_alternative_rows.append(
                    {
                        "label": branch_label,
                        "branch_key": branch_key,
                        "source": "gcf",
                        "alternative_id": key,
                        "value": _format_float(value, digits=4),
                        "split": "",
                    }
                )

        if scfl_row:
            for key, value_key in (
                ("main", "sCF"),
                ("df1", "sDF1"),
                ("df2", "sDF2"),
            ):
                value = _parse_float(scfl_row[value_key])
                branch_alternative_rows.append(
                    {
                        "label": branch_label,
                        "branch_key": branch_key,
                        "source": "scfl",
                        "alternative_id": key,
                        "value": _format_float(value, digits=4),
                        "split": "",
                    }
                )

        if quartet_row:
            for row in quartet_row["rows"]:
                total_weight = float(row["total_weight"])
                value = float(row["weighted_support"]) / total_weight if total_weight > 0 else None
                branch_alternative_rows.append(
                    {
                        "label": branch_label,
                        "branch_key": branch_key,
                        "source": "aster",
                        "alternative_id": str(row["topology_id"]),
                        "value": _format_float(value, digits=6),
                        "split": str(row["split"]),
                    }
                )

    return branch_metric_rows, branch_alternative_rows


def _extract_gene_tree_heterogeneity(
    *,
    repo_root: Path,
    gene_tree_manifest_rows: list[dict[str, str]],
    species_tree_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    species_root = read_newick(species_tree_path)
    species_taxa = leaf_labels(species_root)
    species_topology = canonical_topology_key(species_root)
    species_splits = set(species_topology.split(";")) if species_topology else set()

    heterogeneity_rows: list[dict[str, str]] = []
    topology_examples: dict[str, dict[str, str]] = {}
    topology_counts: Counter = Counter()
    treefile_cache: dict[Path, list[str]] = {}

    for row in gene_tree_manifest_rows:
        tree_path = _resolve_path(repo_root, row["treefile"])
        tree_row_index = (row.get("tree_row_index") or "").strip()
        if tree_row_index:
            tree_lines = treefile_cache.get(tree_path)
            if tree_lines is None:
                tree_lines = [
                    line.strip()
                    for line in tree_path.read_text(encoding="utf-8").splitlines()
                    if line.strip()
                ]
                treefile_cache[tree_path] = tree_lines
            line_index = int(tree_row_index) - 1
            if line_index < 0 or line_index >= len(tree_lines):
                raise ValueError(
                    f"tree_row_index {tree_row_index!r} is out of range for {tree_path}."
                )
            gene_root = parse_newick(tree_lines[line_index])
        else:
            gene_root = read_newick(tree_path)
        gene_taxa = leaf_labels(gene_root)
        is_complete_taxon = gene_taxa == species_taxa
        topology_key = canonical_topology_key(gene_root)
        rf_distance = ""
        matches_species_tree = "false"
        if is_complete_taxon:
            gene_splits = set(topology_key.split(";")) if topology_key else set()
            rf_distance = str(len(species_splits.symmetric_difference(gene_splits)))
            matches_species_tree = str(gene_splits == species_splits).lower()
            topology_counts[topology_key] += 1
            topology_examples.setdefault(
                topology_key,
                {
                    "example_locus_id": row["locus_id"],
                    "matches_species_tree": matches_species_tree,
                },
            )

        heterogeneity_rows.append(
            {
                "locus_id": row["locus_id"],
                "selected_model": row["selected_model"],
                "support_values_present": row["support_values_present"],
                "tip_count": str(len(gene_taxa)),
                "is_complete_taxon": str(is_complete_taxon).lower(),
                "topology_key": topology_key if is_complete_taxon else "",
                "matches_species_tree": matches_species_tree,
                "rf_distance_to_species_tree": rf_distance,
            }
        )

    topology_count_rows: list[dict[str, str]] = []
    total_complete = sum(topology_counts.values())
    for rank, (topology_key, count) in enumerate(topology_counts.most_common(), start=1):
        example = topology_examples[topology_key]
        fraction = count / total_complete if total_complete else 0.0
        topology_count_rows.append(
            {
                "topology_rank": str(rank),
                "topology_id": f"T{rank}",
                "topology_key": topology_key,
                "example_locus_id": example["example_locus_id"],
                "matches_species_tree": example["matches_species_tree"],
                "count": str(count),
                "fraction": _format_float(fraction, digits=6),
            }
        )

    return heterogeneity_rows, topology_count_rows


def _extract_dataset_summary(
    *,
    busco_summary_rows: list[dict[str, str]],
    retained_rows: list[dict[str, str]],
    sample_qc_rows: list[dict[str, str]],
    alignment_summary_rows: list[dict[str, str]],
    gene_tree_heterogeneity_rows: list[dict[str, str]],
    sequence_type: str,
) -> list[dict[str, str]]:
    retained_count = sum(1 for row in retained_rows if row["decision"] == "retain")
    excluded_count = sum(1 for row in retained_rows if row["decision"] != "retain")
    alignment_lengths = [int(row["alignment_length_sites"]) for row in alignment_summary_rows]
    mean_alignment = statistics.fmean(alignment_lengths) if alignment_lengths else 0.0
    median_alignment = statistics.median(alignment_lengths) if alignment_lengths else 0.0
    complete_gene_trees = sum(1 for row in gene_tree_heterogeneity_rows if row["is_complete_taxon"] == "true")
    sample_missing_fractions = [float(row["retained_missing_fraction"]) for row in sample_qc_rows]
    mean_missing_fraction = statistics.fmean(sample_missing_fractions) if sample_missing_fractions else 0.0

    return [
        {
            "sequence_type": sequence_type,
            "sequence_length_unit": sequence_length_unit(sequence_type),
            "sample_count": str(len(busco_summary_rows)),
            "candidate_loci": str(len(retained_rows)),
            "retained_loci": str(retained_count),
            "excluded_loci": str(excluded_count),
            "mean_retained_missing_fraction": _format_float(mean_missing_fraction, digits=6),
            "alignment_count": str(len(alignment_summary_rows)),
            "mean_alignment_length_sites": _format_float(mean_alignment, digits=2),
            "median_alignment_length_sites": _format_float(median_alignment, digits=2),
            "min_alignment_length_sites": str(min(alignment_lengths) if alignment_lengths else 0),
            "max_alignment_length_sites": str(max(alignment_lengths) if alignment_lengths else 0),
            "gene_tree_count": str(len(gene_tree_heterogeneity_rows)),
            "complete_taxon_gene_tree_count": str(complete_gene_trees),
        }
    ]


def build_report_data_bundle(
    *,
    repo_root: Path,
    busco_summary_path: Path,
    locus_matrix_path: Path,
    retained_loci_path: Path,
    gene_tree_manifest_path: Path,
    species_tree_path: Path,
    gcf_stat_path: Path,
    gcf_branch_path: Path,
    scfl_stat_path: Path,
    scfl_branch_path: Path,
    quartet_freqquad_path: Path,
    alignment_dir: Path,
    output_dir: Path,
    sequence_type: str = "protein",
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    busco_summary_rows = read_tsv_rows(busco_summary_path)
    locus_matrix_rows = read_tsv_rows(locus_matrix_path)
    retained_rows = read_tsv_rows(retained_loci_path)
    gene_tree_manifest_rows = read_tsv_rows(gene_tree_manifest_path)
    retained_locus_ids = {row["locus_id"] for row in retained_rows if row["decision"] == "retain"}

    sample_qc_rows = _extract_sample_qc(
        busco_summary_rows=busco_summary_rows,
        locus_matrix_rows=locus_matrix_rows,
        retained_locus_ids=retained_locus_ids,
    )
    alignment_summary_rows = _extract_alignment_summary(
        repo_root=repo_root,
        retained_rows=retained_rows,
        alignment_dir=alignment_dir,
        sequence_type=sequence_type,
    )
    branch_metric_rows, branch_alternative_rows = _extract_branch_tables(
        species_tree_path=species_tree_path,
        output_tree_path=output_dir / "species_tree.report.tre",
        gcf_stat_path=gcf_stat_path,
        gcf_branch_path=gcf_branch_path,
        scfl_stat_path=scfl_stat_path,
        scfl_branch_path=scfl_branch_path,
        quartet_freqquad_path=quartet_freqquad_path,
    )
    heterogeneity_rows, topology_count_rows = _extract_gene_tree_heterogeneity(
        repo_root=repo_root,
        gene_tree_manifest_rows=gene_tree_manifest_rows,
        species_tree_path=species_tree_path,
    )
    dataset_summary_rows = _extract_dataset_summary(
        busco_summary_rows=busco_summary_rows,
        retained_rows=retained_rows,
        sample_qc_rows=sample_qc_rows,
        alignment_summary_rows=alignment_summary_rows,
        gene_tree_heterogeneity_rows=heterogeneity_rows,
        sequence_type=sequence_type,
    )

    write_tsv(
        output_dir / "dataset_summary.tsv",
        dataset_summary_rows,
        [
            "sequence_type",
            "sequence_length_unit",
            "sample_count",
            "candidate_loci",
            "retained_loci",
            "excluded_loci",
            "mean_retained_missing_fraction",
            "alignment_count",
            "mean_alignment_length_sites",
            "median_alignment_length_sites",
            "min_alignment_length_sites",
            "max_alignment_length_sites",
            "gene_tree_count",
            "complete_taxon_gene_tree_count",
        ],
    )
    write_tsv(
        output_dir / "sample_qc.tsv",
        sample_qc_rows,
        [
            "sample_id",
            "taxon_id",
            "sanitized_taxon_id",
            "complete_percent",
            "single_copy_percent",
            "multi_copy_percent",
            "fragmented_percent",
            "missing_percent",
            "internal_stop_codon_count",
            "retained_loci_total",
            "retained_present_loci",
            "retained_missing_loci",
            "retained_missing_fraction",
            "retained_duplicated_loci",
            "retained_fragmented_loci",
            "retained_internal_stop_loci",
        ],
    )
    write_tsv(output_dir / "locus_summary.tsv", retained_rows, list(retained_rows[0].keys()))
    write_tsv(
        output_dir / "alignment_summary.tsv",
        alignment_summary_rows,
        [
            "locus_id",
            "sequence_count",
            "alignment_length_sites",
            "gap_fraction",
            "occupancy",
            "length_dispersion_observed",
        ],
    )
    write_tsv(
        output_dir / "branch_metrics.tsv",
        branch_metric_rows,
        [
            "label",
            "branch_key",
            "split_left",
            "split_right",
            "species_tree_support",
            "species_tree_branch_length",
            "gcf_branch_id",
            "gcf",
            "gdf1",
            "gdf2",
            "gdfp",
            "gn",
            "scfl_branch_id",
            "scfl",
            "sdf1",
            "sdf2",
            "sn",
            "aster_node_id",
            "aster_local_posterior",
            "aster_q1",
            "aster_q2",
            "aster_q3",
            "aster_f1",
            "aster_f2",
            "aster_f3",
            "aster_total_weight",
        ],
    )
    write_tsv(
        output_dir / "branch_alternatives.tsv",
        branch_alternative_rows,
        [
            "label",
            "branch_key",
            "source",
            "alternative_id",
            "value",
            "split",
        ],
    )
    write_tsv(
        output_dir / "gene_tree_heterogeneity.tsv",
        heterogeneity_rows,
        [
            "locus_id",
            "selected_model",
            "support_values_present",
            "tip_count",
            "is_complete_taxon",
            "topology_key",
            "matches_species_tree",
            "rf_distance_to_species_tree",
        ],
    )
    write_tsv(
        output_dir / "topology_counts.tsv",
        topology_count_rows,
        [
            "topology_rank",
            "topology_id",
            "topology_key",
            "example_locus_id",
            "matches_species_tree",
            "count",
            "fraction",
        ],
    )
