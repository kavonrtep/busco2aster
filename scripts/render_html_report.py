"""Render the self-contained HTML visual report from pre-computed TSV data."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import jinja2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True, help="Path to Jinja2 HTML template")
    parser.add_argument("--data-dir", required=True, help="Directory with report TSV files")
    parser.add_argument("--assets-dir", required=True, help="Directory with vendored JS files")
    parser.add_argument("--output", required=True, help="Output HTML file path")
    return parser.parse_args()


def _coerce_value(v: str) -> str | float | int | bool | None:
    """Best-effort coercion of TSV string values to Python types."""
    if v == "":
        return None
    if v in ("true", "True"):
        return True
    if v in ("false", "False"):
        return False
    try:
        if "." in v or "e" in v or "E" in v:
            return float(v)
        return int(v)
    except ValueError:
        return v


def load_tsv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [{k: _coerce_value(v) for k, v in row.items()} for row in reader]


def load_optional_tsv(path: Path) -> list[dict]:
    if path.is_file():
        return load_tsv(path)
    return []


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def build_context(data_dir: Path, assets_dir: Path) -> dict:
    """Build the Jinja2 template context from data files and JS assets."""
    dataset_summary = load_tsv(data_dir / "dataset_summary.tsv")
    summary = dataset_summary[0] if dataset_summary else {}

    sample_qc = load_tsv(data_dir / "sample_qc.tsv")
    locus_summary = load_tsv(data_dir / "locus_summary.tsv")
    alignment_summary = load_tsv(data_dir / "alignment_summary.tsv")
    branch_metrics = load_tsv(data_dir / "branch_metrics.tsv")
    branch_alternatives = load_tsv(data_dir / "branch_alternatives.tsv")
    gene_tree_heterogeneity = load_tsv(data_dir / "gene_tree_heterogeneity.tsv")
    topology_counts = load_tsv(data_dir / "topology_counts.tsv")

    tree_path = data_dir / "species_tree.report.tre"
    species_tree_newick = tree_path.read_text(encoding="utf-8").strip() if tree_path.is_file() else ""

    topo_tests_available = (data_dir / "au_test_results.tsv").is_file()
    au_test_results = load_optional_tsv(data_dir / "au_test_results.tsv")
    branch_quartet_support = load_optional_tsv(data_dir / "branch_quartet_support.tsv")
    contested_branches = load_optional_tsv(data_dir / "contested_branches.tsv")

    return {
        "summary": summary,
        "sample_qc_json": json.dumps(sample_qc),
        "locus_summary_json": json.dumps(locus_summary),
        "alignment_summary_json": json.dumps(alignment_summary),
        "branch_metrics": branch_metrics,
        "branch_metrics_json": json.dumps(branch_metrics),
        "branch_alternatives_json": json.dumps(branch_alternatives),
        "gene_tree_heterogeneity_json": json.dumps(gene_tree_heterogeneity),
        "topology_counts_json": json.dumps(topology_counts),
        "species_tree_newick": species_tree_newick,
        "topo_tests_available": topo_tests_available,
        "au_test_results": au_test_results,
        "branch_quartet_support": branch_quartet_support,
        "contested_branches": contested_branches,
        "plotly_js": load_text(assets_dir / "plotly-basic.min.js"),
        "underscore_js": load_text(assets_dir / "underscore.min.js"),
        "lodash_js": load_text(assets_dir / "lodash.min.js"),
        "phylotree_js": load_text(assets_dir / "phylotree.js"),
    }


def render(template_path: Path, context: dict, output_path: Path) -> None:
    template_text = template_path.read_text(encoding="utf-8")
    env = jinja2.Environment(autoescape=False, undefined=jinja2.StrictUndefined)
    template = env.from_string(template_text)
    html = template.render(**context)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")


def main() -> int:
    args = parse_args()
    data_dir = Path(args.data_dir).resolve()
    assets_dir = Path(args.assets_dir).resolve()
    template_path = Path(args.template).resolve()
    output_path = Path(args.output).resolve()

    context = build_context(data_dir, assets_dir)
    render(template_path, context, output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
