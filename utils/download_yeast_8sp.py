#!/usr/bin/env python3
"""
Download an 8-species yeast test dataset for busco2aster.

Species selection rationale
----------------------------
All eight assemblies are small (9.8–12.7 Mb), have chromosome-level or
high-quality scaffold assemblies, and span multiple orders of Ascomycota:

  Saccharomycetales ──┬── Saccharomyces cerevisiae  S288C        GCF_000146045.2 ~12.1 Mb
                      ├── Saccharomyces paradoxus   NRRL Y-17217 GCF_002079055.1 ~12.0 Mb
                      ├── Zygosaccharomyces rouxii  CBS 732      GCF_000026365.1  ~9.8 Mb
                      ├── Lachancea kluyveri        NRRL Y-12651 GCA_000149225.2 ~11.5 Mb
                      ├── Lachancea thermotolerans  CBS 6340     GCF_000142805.1 ~10.4 Mb
                      ├── Kluyveromyces marxianus   DMKU3-1042   GCA_028768585.1 ~11.0 Mb
                      └── Nakaseomyces glabratus    ATCC 2001    GCF_010111755.1 ~12.7 Mb
  Schizosaccharomycetales (outgroup)
                      └── Schizosaccharomyces pombe 972h-        GCF_000002945.2 ~12.6 Mb

BUSCO lineage: eukaryota_odb12 (129 single-copy markers, OrthoDB v12).
Smallest available lineage — ~18x fewer BUSCOs than saccharomycetes_odb12 (2319).
All 8 species carry the core eukaryotic gene set, giving enough loci for wASTRAL
even at 0.8 occupancy. Upgrade to fungi_odb12 (1122) or saccharomycetes_odb12
(2319) for more rigorous runs.

Total download: ~95 Mb compressed genomes.

Usage
-----
    python3 utils/download_yeast_8sp.py
    python3 utils/download_yeast_8sp.py --output-dir /scratch/yeast8
    python3 utils/download_yeast_8sp.py --dry-run
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

# ── Species manifest ──────────────────────────────────────────────────────────

SPECIES: list[dict[str, str]] = [
    {
        "sample_id": "scerevisiae",
        "taxon_id": "Saccharomyces_cerevisiae",
        "accession": "GCF_000146045.2",
        "strain": "S288C",
        "note": "Saccharomycetales; canonical reference",
    },
    {
        "sample_id": "sparadoxus",
        "taxon_id": "Saccharomyces_paradoxus",
        "accession": "GCF_002079055.1",
        "strain": "NRRL Y-17217",
        "note": "Saccharomycetales; sister to S. cerevisiae",
    },
    {
        "sample_id": "zrouxii",
        "taxon_id": "Zygosaccharomyces_rouxii",
        "accession": "GCF_000026365.1",
        "strain": "CBS 732",
        "note": "Saccharomycetales; osmotolerant",
    },
    {
        "sample_id": "lkluyveri",
        "taxon_id": "Lachancea_kluyveri",
        "accession": "GCA_000149225.2",
        "strain": "NRRL Y-12651",
        "note": "Saccharomycetales; Lachancea clade",
    },
    {
        "sample_id": "lthermotolerans",
        "taxon_id": "Lachancea_thermotolerans",
        "accession": "GCF_000142805.1",
        "strain": "CBS 6340",
        "note": "Saccharomycetales; Lachancea clade",
    },
    {
        "sample_id": "kmarxianus",
        "taxon_id": "Kluyveromyces_marxianus",
        "accession": "GCA_028768585.1",
        "strain": "DMKU3-1042",
        "note": "Saccharomycetales; dairy yeast",
    },
    {
        "sample_id": "nglabrata",
        "taxon_id": "Nakaseomyces_glabratus",
        "accession": "GCF_010111755.1",
        "strain": "ATCC 2001",
        "note": "Saccharomycetales; reference genome",
    },
    {
        "sample_id": "spombe",
        "taxon_id": "Schizosaccharomyces_pombe",
        "accession": "GCF_000002945.2",
        "strain": "972h-",
        "note": "Schizosaccharomycetales; outgroup",
    },
]

NCBI_DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
NCBI_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
REQUEST_DELAY = 0.4  # seconds between NCBI API calls


# ── FTP URL construction ──────────────────────────────────────────────────────

def ncbi_ftp_genome_url(accession: str, assembly_name: str) -> str:
    """Construct NCBI FTP HTTPS URL for the genomic FASTA of an assembly."""
    prefix, num_ver = accession.split("_", 1)
    num = num_ver.split(".")[0]  # e.g. "000146045"
    path_chunks = "/".join([num[i: i + 3] for i in range(0, 9, 3)])
    dir_name = f"{accession}_{assembly_name}"
    filename = f"{dir_name}_genomic.fna.gz"
    return f"{NCBI_FTP_BASE}/{prefix}/{path_chunks}/{dir_name}/{filename}"


def fetch_assembly_name(accession: str) -> str:
    """Return the NCBI assembly_name for an accession (e.g. 'ASM294v3')."""
    url = f"{NCBI_DATASETS_API}/genome/accession/{accession}/dataset_report"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.load(resp)
        return data["reports"][0]["assembly_info"]["assembly_name"]
    except Exception as exc:
        raise RuntimeError(
            f"Could not fetch assembly_name for {accession}: {exc}"
        ) from exc


# ── Download helper ───────────────────────────────────────────────────────────

def _reporthook(count: int, block_size: int, total_size: int) -> None:
    if total_size <= 0:
        print(f"\r  downloaded {count * block_size // 1024} KB ...", end="", flush=True)
    else:
        pct = min(100, count * block_size * 100 // total_size)
        mb_done = count * block_size / 1_048_576
        mb_total = total_size / 1_048_576
        print(f"\r  {pct:3d}%  {mb_done:.1f}/{mb_total:.1f} MB", end="", flush=True)


def download_file(url: str, dest: Path, *, dry_run: bool = False) -> None:
    """Download url to dest, showing progress. Skip if dest already exists."""
    if dest.is_file():
        print(f"  already present: {dest.name}")
        return
    if dry_run:
        print(f"  [dry-run] would download: {url}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".part")
    print(f"  → {dest.name}")
    try:
        urllib.request.urlretrieve(url, tmp, reporthook=_reporthook)
        print()  # newline after progress bar
        tmp.rename(dest)
    except Exception:
        if tmp.is_file():
            tmp.unlink()
        raise


# ── Manifest writers ──────────────────────────────────────────────────────────

def write_samples_tsv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("sample_id\ttaxon_id\tassembly_fasta\n")
        for row in rows:
            fh.write(f"{row['sample_id']}\t{row['taxon_id']}\t{row['assembly_fasta']}\n")
    print(f"Wrote samples manifest: {path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--output-dir",
        default=str(repo_root / "test_datasets" / "yeast_8sp"),
        help="Directory to write downloaded genomes and samples.tsv. "
             "Created if absent. Default: <repo>/test_datasets/yeast_8sp",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print URLs that would be downloaded; do not fetch.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    genomes_dir = output_dir / "genomes"
    dry_run = args.dry_run

    print(f"Output directory : {output_dir}")
    print(f"Genomes directory: {genomes_dir}")
    if dry_run:
        print("Mode             : DRY-RUN (no files will be written)\n")
    else:
        genomes_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, str]] = []
    errors: list[str] = []

    for sp in SPECIES:
        accession = sp["accession"]
        sample_id = sp["sample_id"]
        print(f"\n[{sample_id}] {sp['taxon_id'].replace('_', ' ')} ({accession}, {sp['strain']})")

        try:
            print("  Fetching assembly name from NCBI ...", end="", flush=True)
            assembly_name = fetch_assembly_name(accession)
            print(f" {assembly_name}")
        except RuntimeError as exc:
            print(f" FAILED: {exc}")
            errors.append(str(exc))
            continue

        time.sleep(REQUEST_DELAY)

        fasta_url = ncbi_ftp_genome_url(accession, assembly_name)
        dest = genomes_dir / f"{sample_id}.fna.gz"

        try:
            download_file(fasta_url, dest, dry_run=dry_run)
        except Exception as exc:
            print(f"  FAILED: {exc}")
            errors.append(f"{sample_id}: {exc}")
            continue

        manifest_rows.append(
            {
                "sample_id": sample_id,
                "taxon_id": sp["taxon_id"],
                "assembly_fasta": str(dest),
            }
        )

    print()
    if errors:
        print(f"WARNING: {len(errors)} download(s) failed:")
        for err in errors:
            print(f"  - {err}")
        print()

    if manifest_rows and not dry_run:
        samples_tsv = output_dir / "samples.tsv"
        write_samples_tsv(manifest_rows, samples_tsv)

    if not dry_run:
        repo_root = Path(__file__).resolve().parents[1]
        config_src = repo_root / "utils" / "yeast_8sp.yaml"
        print()
        print("Next steps:")
        print(f"  1. Verify the manifest:  snakemake --cores 1 --configfile {config_src} results/metadata/samples.validated.tsv")
        print(f"  2. Run the pipeline:     python3 run_pipeline.py -c {config_src} --repo-root . --directory . -t 8")
        print()
        print("Or with the container:")
        print(f"  apptainer run -B $(pwd) -B {output_dir} busco2aster.sif -c {config_src} --repo-root $(pwd) --directory $(pwd) -t 8")

    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
