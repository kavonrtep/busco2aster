"""Normalize assembly FASTA files into a wrapped, gzipped cache using SeqKit."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


ASSEMBLY_PREP_QC_FIELDNAMES = [
    "input_path",
    "output_path",
    "seqkit_executable",
    "seqkit_version",
    "threads",
    "wrap_width",
    "input_bytes",
    "output_bytes",
    "record_count",
    "sequence_character_count",
    "input_format",
    "input_type",
    "input_num_seqs",
    "input_sum_len",
    "input_min_len",
    "input_avg_len",
    "input_max_len",
    "output_format",
    "output_type",
    "output_num_seqs",
    "output_sum_len",
    "output_min_len",
    "output_avg_len",
    "output_max_len",
    "count_match",
    "length_match",
]


class AssemblyPreparationError(ValueError):
    """Raised when assembly preparation fails."""


def _write_qc_tsv(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=ASSEMBLY_PREP_QC_FIELDNAMES,
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerow(row)


def _run_command(command: list[str], context: str) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(command, capture_output=True, text=True, check=True)
    except FileNotFoundError as exc:
        raise AssemblyPreparationError(
            f"{context} failed because {command[0]!r} is not available on PATH."
        ) from exc
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        stdout = (exc.stdout or "").strip()
        detail = stderr or stdout or f"exit code {exc.returncode}"
        raise AssemblyPreparationError(f"{context} failed: {detail}") from exc


def build_seqkit_seq_command(
    input_path: Path,
    output_path: Path,
    wrap_width: int,
    threads: int,
    seqkit_executable: str = "seqkit",
) -> list[str]:
    return [
        seqkit_executable,
        "seq",
        "-j",
        str(threads),
        "-w",
        str(wrap_width),
        "-o",
        output_path.as_posix(),
        input_path.as_posix(),
    ]


def build_seqkit_stats_command(
    path: Path,
    threads: int,
    seqkit_executable: str = "seqkit",
) -> list[str]:
    return [
        seqkit_executable,
        "stats",
        "-T",
        "-j",
        str(threads),
        path.as_posix(),
    ]


def parse_seqkit_stats_tsv(text: str) -> dict[str, str]:
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if len(rows) != 1:
        raise AssemblyPreparationError(
            f"Expected one row from `seqkit stats -T`, but received {len(rows)} rows."
        )
    row = rows[0]
    required = ("file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len")
    missing = [column for column in required if column not in row or row[column] == ""]
    if missing:
        raise AssemblyPreparationError(
            "SeqKit stats output is missing required columns: " + ", ".join(missing)
        )
    return row


def query_seqkit_stats(
    path: Path,
    threads: int,
    seqkit_executable: str = "seqkit",
) -> dict[str, str]:
    result = _run_command(
        build_seqkit_stats_command(path, threads, seqkit_executable),
        context=f"SeqKit stats for {path}",
    )
    return parse_seqkit_stats_tsv(result.stdout)


def query_seqkit_version(seqkit_executable: str = "seqkit") -> str:
    result = _run_command([seqkit_executable, "version"], context="SeqKit version check")
    version = result.stdout.strip() or result.stderr.strip()
    if not version:
        raise AssemblyPreparationError("SeqKit version output was empty.")
    return version


def prepare_assembly(
    input_path: Path,
    output_path: Path,
    qc_output_path: Path,
    wrap_width: int,
    threads: int = 1,
    seqkit_executable: str = "seqkit",
) -> dict[str, str]:
    if wrap_width <= 0:
        raise AssemblyPreparationError("Wrap width must be a positive integer.")
    if threads <= 0:
        raise AssemblyPreparationError("Thread count must be a positive integer.")
    if not input_path.is_file():
        raise AssemblyPreparationError(f"Assembly input does not exist: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    seqkit_version = query_seqkit_version(seqkit_executable)
    input_stats = query_seqkit_stats(input_path, threads, seqkit_executable)

    _run_command(
        build_seqkit_seq_command(
            input_path=input_path,
            output_path=output_path,
            wrap_width=wrap_width,
            threads=threads,
            seqkit_executable=seqkit_executable,
        ),
        context=f"SeqKit rewrite for {input_path}",
    )

    output_stats = query_seqkit_stats(output_path, threads, seqkit_executable)
    count_match = input_stats["num_seqs"] == output_stats["num_seqs"]
    length_match = input_stats["sum_len"] == output_stats["sum_len"]
    if not count_match or not length_match:
        raise AssemblyPreparationError(
            "Prepared assembly statistics do not match the input assembly: "
            f"num_seqs {input_stats['num_seqs']} -> {output_stats['num_seqs']}, "
            f"sum_len {input_stats['sum_len']} -> {output_stats['sum_len']}"
        )

    qc_row = {
        "input_path": input_path.as_posix(),
        "output_path": output_path.as_posix(),
        "seqkit_executable": seqkit_executable,
        "seqkit_version": seqkit_version,
        "threads": str(threads),
        "wrap_width": str(wrap_width),
        "input_bytes": str(input_path.stat().st_size),
        "output_bytes": str(output_path.stat().st_size),
        "record_count": input_stats["num_seqs"],
        "sequence_character_count": input_stats["sum_len"],
        "input_format": input_stats["format"],
        "input_type": input_stats["type"],
        "input_num_seqs": input_stats["num_seqs"],
        "input_sum_len": input_stats["sum_len"],
        "input_min_len": input_stats["min_len"],
        "input_avg_len": input_stats["avg_len"],
        "input_max_len": input_stats["max_len"],
        "output_format": output_stats["format"],
        "output_type": output_stats["type"],
        "output_num_seqs": output_stats["num_seqs"],
        "output_sum_len": output_stats["sum_len"],
        "output_min_len": output_stats["min_len"],
        "output_avg_len": output_stats["avg_len"],
        "output_max_len": output_stats["max_len"],
        "count_match": "true" if count_match else "false",
        "length_match": "true" if length_match else "false",
    }
    _write_qc_tsv(qc_output_path, qc_row)
    return qc_row


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Original assembly FASTA or FASTA.gz file.")
    parser.add_argument("--output", required=True, help="Prepared gzipped FASTA output path.")
    parser.add_argument("--qc-output", required=True, help="Per-sample QC TSV output path.")
    parser.add_argument(
        "--wrap-width",
        type=int,
        default=80,
        help="Sequence line width for the prepared FASTA. Default: 80.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of SeqKit threads to use. Default: 1.",
    )
    parser.add_argument(
        "--seqkit-executable",
        default="seqkit",
        help="SeqKit executable to use. Default: seqkit.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    prepare_assembly(
        input_path=Path(args.input),
        output_path=Path(args.output),
        qc_output_path=Path(args.qc_output),
        wrap_width=args.wrap_width,
        threads=args.threads,
        seqkit_executable=args.seqkit_executable,
    )


if __name__ == "__main__":
    main()
