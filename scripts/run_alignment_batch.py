"""Run one batch of MAFFT alignments from the retained-locus FASTA manifest."""

from __future__ import annotations

import argparse
from pathlib import Path

from .alignment import run_alignment_batch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-dir", required=True)
    parser.add_argument("--batch-id", required=True)
    parser.add_argument("--batch-size", required=True, type=int)
    parser.add_argument("--threads-per-alignment", required=True, type=int)
    parser.add_argument("--sequence-type", default="protein")
    parser.add_argument("--mafft-executable", default="mafft")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run_alignment_batch(
        manifest_path=Path(args.manifest),
        output_dir=Path(args.output_dir),
        log_dir=Path(args.log_dir),
        batch_id=args.batch_id,
        batch_size=args.batch_size,
        threads_per_alignment=args.threads_per_alignment,
        sequence_type=args.sequence_type,
        mafft_executable=args.mafft_executable,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
