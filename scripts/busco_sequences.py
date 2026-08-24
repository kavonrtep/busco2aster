"""Layout-aware access to BUSCO per-locus sequences.

BUSCO writes one ``.faa`` and one ``.gff`` per locus per category, which costs
two inodes per locus per sample (360,088 files for a 23-genome run). The
workflow instead concatenates each category into a single pair of files per
sample:

``results/busco/<sample>/busco_sequences/single_copy.faa``
    One FASTA record per sequence, the header rewritten to start with the BUSCO
    id so records stay addressable: ``>{busco_id} {original header}``.

``results/busco/<sample>/busco_sequences/single_copy.gff``
    Each locus's GFF copied verbatim behind a ``##busco_id <id>`` marker line.
    The marker is a comment, so parsers that skip ``#`` lines are unaffected.

Runs produced before this layout recorded a *directory* per category in
``paths.tsv`` and are still readable: every accessor here takes the path the
tables recorded and picks the layout from it, so an existing ``results/busco/``
tree keeps working without re-running BUSCO.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

GFF_MARKER_PREFIX = "##busco_id "
CONCATENATED_LAYOUT = "concatenated"
PER_LOCUS_LAYOUT = "per_locus"

# Every multi-locus file the workflow writes, by stem. `standardize_busco_run`
# writes the three BUSCO categories, `extract_retained_sample_dna` writes
# `retained_loci.fna`. Anything else a table points at is a per-locus file from
# a run made before this layout, and holds exactly one locus.
CONCATENATED_FILE_STEMS = frozenset(
    {"single_copy", "multi_copy", "fragmented", "retained_loci"}
)

CATEGORY_BY_STATUS = {
    "Complete": "single_copy",
    "Duplicated": "multi_copy",
    "Fragmented": "fragmented",
}
LEGACY_DIR_BY_CATEGORY = {
    "single_copy": "single_copy_busco_sequences",
    "multi_copy": "multi_copy_busco_sequences",
    "fragmented": "fragmented_busco_sequences",
}
FAA_ARTIFACT_BY_CATEGORY = {
    category: f"{category}_busco_sequences" for category in LEGACY_DIR_BY_CATEGORY
}
GFF_ARTIFACT_BY_CATEGORY = {
    category: f"{category}_busco_gff" for category in LEGACY_DIR_BY_CATEGORY
}


class BuscoSequenceError(ValueError):
    """Raised when a BUSCO sequence artifact is missing or unreadable."""


def resolve_repo_path(repo_root: Path, path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else repo_root / path


def parse_fasta_records_text(text: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence_lines: list[str] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                records.append((header, "".join(sequence_lines)))
            header = line[1:]
            sequence_lines = []
            continue
        sequence_lines.append(line)

    if header:
        records.append((header, "".join(sequence_lines)))
    return records


def is_concatenated_artifact(path_text: str) -> bool:
    """True when `paths.tsv` recorded a concatenated file rather than a directory."""
    return path_text.endswith((".faa", ".gff"))


def artifact_layout(path_text: str) -> str:
    return CONCATENATED_LAYOUT if is_concatenated_artifact(path_text) else PER_LOCUS_LAYOUT


def is_concatenated_file(path: Path) -> bool:
    return path.stem in CONCATENATED_FILE_STEMS


def _cache_key(path: Path) -> tuple[str, int, int]:
    stat = path.stat()
    return (path.as_posix(), stat.st_size, stat.st_mtime_ns)


@lru_cache(maxsize=None)
def _fasta_offsets(key: tuple[str, int, int]) -> dict[str, list[tuple[int, int]]]:
    """Map each BUSCO id to the byte ranges of its records.

    Only offsets are cached, never sequences, so the index costs about a
    megabyte per sample regardless of how large the concatenated file is.
    """
    offsets: dict[str, list[tuple[int, int]]] = {}
    current_id: str | None = None
    start = 0
    position = 0

    with Path(key[0]).open("rb") as handle:
        for raw_line in handle:
            if raw_line.startswith(b">"):
                if current_id is not None:
                    offsets.setdefault(current_id, []).append((start, position))
                fields = raw_line[1:].split()
                current_id = fields[0].decode("utf-8") if fields else ""
                start = position
            position += len(raw_line)

    if current_id is not None:
        offsets.setdefault(current_id, []).append((start, position))
    return offsets


@lru_cache(maxsize=None)
def _gff_offsets(key: tuple[str, int, int]) -> dict[str, tuple[int, int]]:
    offsets: dict[str, tuple[int, int]] = {}
    marker = GFF_MARKER_PREFIX.encode("utf-8")
    current_id: str | None = None
    start = 0
    position = 0

    with Path(key[0]).open("rb") as handle:
        for raw_line in handle:
            if raw_line.startswith(marker):
                if current_id is not None:
                    offsets[current_id] = (start, position)
                current_id = raw_line[len(marker) :].strip().decode("utf-8")
                start = position + len(raw_line)
            position += len(raw_line)

    if current_id is not None:
        offsets[current_id] = (start, position)
    return offsets


def _read_slice(path: Path, start: int, end: int) -> str:
    with path.open("rb") as handle:
        handle.seek(start)
        return handle.read(end - start).decode("utf-8")


def read_fasta_records_for_id(
    fasta_path_text: str,
    busco_id: str,
    repo_root: Path,
) -> list[tuple[str, str]]:
    """Return every FASTA record stored under `busco_id`.

    Duplicated loci carry more than one record, which is what the sequence-count
    QC in the locus matrix checks, so this always returns a list.
    """
    path = resolve_repo_path(repo_root, fasta_path_text)
    if not path.is_file():
        raise BuscoSequenceError(f"Expected FASTA does not exist: {path}")

    if not is_concatenated_file(path):
        return parse_fasta_records_text(path.read_text(encoding="utf-8"))

    ranges = _fasta_offsets(_cache_key(path)).get(busco_id, [])
    records: list[tuple[str, str]] = []
    for start, end in ranges:
        records.extend(parse_fasta_records_text(_read_slice(path, start, end)))
    return records


def read_gff_text(gff_path_text: str, busco_id: str, repo_root: Path) -> str:
    """Return the GFF lines recorded for `busco_id`, verbatim."""
    path = resolve_repo_path(repo_root, gff_path_text)
    if not path.is_file():
        raise BuscoSequenceError(f"Expected BUSCO GFF does not exist: {path}")

    if not is_concatenated_file(path):
        return path.read_text(encoding="utf-8")

    block = _gff_offsets(_cache_key(path)).get(busco_id)
    if block is None:
        raise BuscoSequenceError(f"Concatenated GFF {path} has no block for locus {busco_id!r}.")
    return _read_slice(path, *block)


def has_fasta_records_for_id(fasta_path_text: str, busco_id: str, repo_root: Path) -> bool:
    path = resolve_repo_path(repo_root, fasta_path_text)
    if not path.is_file():
        return False
    if not is_concatenated_file(path):
        return True
    return busco_id in _fasta_offsets(_cache_key(path))


def has_gff_block(gff_path_text: str, busco_id: str, repo_root: Path) -> bool:
    path = resolve_repo_path(repo_root, gff_path_text)
    if not path.is_file():
        return False
    if not is_concatenated_file(path):
        return True
    return busco_id in _gff_offsets(_cache_key(path))


def concatenated_category_paths(sequence_root: Path, category: str) -> tuple[Path, Path]:
    return sequence_root / f"{category}.faa", sequence_root / f"{category}.gff"


def _write_concatenated_fasta(sources: list[Path], target: Path) -> int:
    record_count = 0
    with target.open("w", encoding="utf-8") as handle:
        for source in sources:
            busco_id = source.stem
            for raw_line in source.read_text(encoding="utf-8").splitlines():
                if raw_line.startswith(">"):
                    record_count += 1
                    original = raw_line[1:].strip()
                    first_token = original.split()[0] if original else ""
                    if first_token == busco_id:
                        handle.write(f">{original}\n")
                    else:
                        handle.write(f">{busco_id} {original}".rstrip() + "\n")
                    continue
                handle.write(raw_line + "\n")
    return record_count


def _write_concatenated_gff(sources: list[Path], target: Path) -> int:
    with target.open("w", encoding="utf-8") as handle:
        for source in sources:
            handle.write(f"{GFF_MARKER_PREFIX}{source.stem}\n")
            text = source.read_text(encoding="utf-8")
            handle.write(text if text.endswith("\n") or not text else text + "\n")
    return len(sources)


def concatenate_busco_category(source_dir: Path, sequence_root: Path, category: str) -> dict[str, str]:
    """Fold one BUSCO category directory into a `.faa`/`.gff` pair.

    Returns the written paths, which is what `paths.tsv` records; the per-locus
    filenames are preserved as the BUSCO id on every record.
    """
    faa_sources = sorted(source_dir.glob("*.faa"))
    gff_sources = sorted(source_dir.glob("*.gff"))
    faa_target, gff_target = concatenated_category_paths(sequence_root, category)
    sequence_root.mkdir(parents=True, exist_ok=True)

    record_count = _write_concatenated_fasta(faa_sources, faa_target)
    block_count = _write_concatenated_gff(gff_sources, gff_target)

    return {
        "faa_path": faa_target.as_posix(),
        "gff_path": gff_target.as_posix(),
        "locus_count": str(len(faa_sources)),
        "record_count": str(record_count),
        "gff_block_count": str(block_count),
    }
