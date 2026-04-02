"""Helpers for switching the workflow between protein and DNA modes."""

from __future__ import annotations

from pathlib import Path


SUPPORTED_SEQUENCE_TYPES = {"protein", "dna"}


def normalize_sequence_type(value: str | None) -> str:
    normalized = (value or "protein").strip().lower()
    if normalized not in SUPPORTED_SEQUENCE_TYPES:
        supported = ", ".join(sorted(SUPPORTED_SEQUENCE_TYPES))
        raise ValueError(f"Unsupported sequence_type {value!r}; expected one of: {supported}")
    return normalized


def raw_fasta_suffix(sequence_type: str) -> str:
    return "fna" if normalize_sequence_type(sequence_type) == "dna" else "faa"


def alignment_suffix(sequence_type: str) -> str:
    return f"aln.{raw_fasta_suffix(sequence_type)}"


def iqtree_seqtype(sequence_type: str) -> str:
    return "DNA" if normalize_sequence_type(sequence_type) == "dna" else "AA"


def scfl_default_model(sequence_type: str) -> str:
    return "GTR+G" if normalize_sequence_type(sequence_type) == "dna" else "LG+G4"


def resolve_scfl_model(sequence_type: str, configured_model: str | None) -> str:
    normalized_type = normalize_sequence_type(sequence_type)
    default_model = scfl_default_model(normalized_type)
    model_text = (configured_model or "").strip()
    if not model_text or model_text.lower() == "auto":
        return default_model
    # Backward compatibility: older configs baked in the protein default,
    # which is invalid for DNA-mode sCF runs.
    if normalized_type == "dna" and model_text == "LG+G4":
        return default_model
    return model_text


def mafft_mode_arguments(sequence_type: str) -> list[str]:
    return [] if normalize_sequence_type(sequence_type) == "dna" else ["--amino"]


def sequence_length_unit(sequence_type: str) -> str:
    return "nt" if normalize_sequence_type(sequence_type) == "dna" else "aa"


def raw_fasta_path(output_dir: Path, locus_id: str, sequence_type: str) -> Path:
    return output_dir / f"{locus_id}.{raw_fasta_suffix(sequence_type)}"


def alignment_path(output_dir: Path, locus_id: str, sequence_type: str) -> Path:
    return output_dir / f"{locus_id}.{alignment_suffix(sequence_type)}"


def locus_id_from_alignment_filename(path_text: str) -> str:
    name = Path(path_text).name
    for suffix in (".aln.faa", ".aln.fna"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path_text).stem
