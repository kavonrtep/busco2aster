import csv
import gzip
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path

from scripts.prepare_assembly import (
    AssemblyPreparationError,
    build_seqkit_seq_command,
    build_seqkit_stats_command,
    parse_seqkit_stats_tsv,
    prepare_assembly,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


class AssemblyPrepUnitTests(unittest.TestCase):
    def make_fake_seqkit(self, directory: Path) -> Path:
        fake_path = directory / "seqkit"
        fake_path.write_text(
            """#!/usr/bin/env python3
import csv
import gzip
import sys
from pathlib import Path


def open_text(path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', newline='')
    return open(path, 'rt', encoding='utf-8', newline='')


def open_write(path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'wt', encoding='utf-8', newline='')
    return open(path, 'wt', encoding='utf-8', newline='')


def load_records(path):
    records = []
    header = None
    chunks = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip('\\n').rstrip('\\r')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(chunks)))
                header = line[1:].strip()
                chunks = []
            else:
                if header is None:
                    raise SystemExit('invalid FASTA')
                chunks.append(''.join(line.split()))
    if header is not None:
        records.append((header, ''.join(chunks)))
    return records


argv = sys.argv[1:]
if argv == ['version']:
    print('seqkit v9.9.9-test')
    raise SystemExit(0)

if argv and argv[0] == 'stats':
    path = Path(argv[-1])
    records = load_records(path)
    if not records:
        raise SystemExit('invalid FASTA')
    lengths = [len(seq) for _, seq in records]
    writer = csv.writer(sys.stdout, delimiter='\\t', lineterminator='\\n')
    writer.writerow(['file', 'format', 'type', 'num_seqs', 'sum_len', 'min_len', 'avg_len', 'max_len'])
    writer.writerow([
        path.as_posix(),
        'FASTA',
        'DNA',
        str(len(records)),
        str(sum(lengths)),
        str(min(lengths)),
        f"{sum(lengths) / len(lengths):.1f}",
        str(max(lengths)),
    ])
    raise SystemExit(0)

if argv and argv[0] == 'seq':
    out_idx = argv.index('-o')
    output_path = Path(argv[out_idx + 1])
    width = int(argv[argv.index('-w') + 1])
    input_path = Path(argv[-1])
    records = load_records(input_path)
    if not records:
        raise SystemExit('invalid FASTA')
    with open_write(output_path) as out_handle:
        for header, sequence in records:
            out_handle.write(f'>{header}\\n')
            for start in range(0, len(sequence), width):
                out_handle.write(sequence[start:start+width] + '\\n')
    raise SystemExit(0)

raise SystemExit('unsupported fake seqkit call: ' + ' '.join(sys.argv))
""",
            encoding="utf-8",
        )
        fake_path.chmod(fake_path.stat().st_mode | stat.S_IXUSR)
        return fake_path

    def test_build_seqkit_seq_command_renders_expected_flags(self):
        command = build_seqkit_seq_command(
            input_path=Path("input.fa.gz"),
            output_path=Path("output.fa.gz"),
            wrap_width=80,
            threads=4,
            seqkit_executable="seqkit",
        )
        self.assertEqual(
            command,
            ["seqkit", "seq", "-j", "4", "-w", "80", "-o", "output.fa.gz", "input.fa.gz"],
        )

    def test_build_seqkit_stats_command_renders_expected_flags(self):
        command = build_seqkit_stats_command(
            path=Path("input.fa.gz"),
            threads=3,
            seqkit_executable="seqkit",
        )
        self.assertEqual(
            command,
            ["seqkit", "stats", "-T", "-j", "3", "input.fa.gz"],
        )

    def test_parse_seqkit_stats_tsv_reads_expected_columns(self):
        parsed = parse_seqkit_stats_tsv(
            "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n"
            "input.fa.gz\tFASTA\tDNA\t2\t140\t20\t70.0\t120\n"
        )
        self.assertEqual(parsed["file"], "input.fa.gz")
        self.assertEqual(parsed["num_seqs"], "2")
        self.assertEqual(parsed["sum_len"], "140")

    def test_prepare_assembly_wraps_long_sequence_lines_and_writes_qc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input.fa"
            output_path = root / "prepared.fa.gz"
            qc_path = root / "prepared.prep.tsv"
            fake_seqkit = self.make_fake_seqkit(root)

            input_path.write_text(
                ">seq1\n" + ("ACGT" * 30) + "\n>seq2\n" + ("NN" * 10) + "\n",
                encoding="utf-8",
            )

            qc_row = prepare_assembly(
                input_path,
                output_path,
                qc_path,
                wrap_width=25,
                threads=2,
                seqkit_executable=fake_seqkit.as_posix(),
            )

            with gzip.open(output_path, "rt", encoding="utf-8") as handle:
                output_lines = handle.read().splitlines()
            with qc_path.open(newline="", encoding="utf-8") as handle:
                qc_rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(output_lines[0], ">seq1")
            self.assertEqual([len(line) for line in output_lines[1:6]], [25, 25, 25, 25, 20])
            self.assertEqual("".join(output_lines[1:6]), "ACGT" * 30)
            self.assertEqual(output_lines[6], ">seq2")
            self.assertEqual(output_lines[7], "NNNNNNNNNNNNNNNNNNNN")
            self.assertEqual(qc_row["record_count"], "2")
            self.assertEqual(qc_row["sequence_character_count"], "140")
            self.assertEqual(qc_row["wrap_width"], "25")
            self.assertEqual(qc_row["threads"], "2")
            self.assertEqual(qc_row["seqkit_version"], "seqkit v9.9.9-test")
            self.assertEqual(qc_row["input_max_len"], "120")
            self.assertEqual(qc_row["output_max_len"], "120")
            self.assertEqual(qc_row["count_match"], "true")
            self.assertEqual(qc_row["length_match"], "true")
            self.assertEqual(len(qc_rows), 1)
            self.assertTrue(qc_path.is_file())

    def test_prepare_assembly_supports_gzipped_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input.fa.gz"
            output_path = root / "prepared.fa.gz"
            qc_path = root / "prepared.prep.tsv"
            fake_seqkit = self.make_fake_seqkit(root)

            with gzip.open(input_path, "wt", encoding="utf-8") as handle:
                handle.write(">seq1\nACGT\n>seq2\nTGCA\n")

            qc_row = prepare_assembly(
                input_path,
                output_path,
                qc_output_path=qc_path,
                wrap_width=80,
                threads=1,
                seqkit_executable=fake_seqkit.as_posix(),
            )

            with gzip.open(output_path, "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), ">seq1\nACGT\n>seq2\nTGCA\n")

            self.assertEqual(qc_row["record_count"], "2")
            self.assertEqual(qc_row["count_match"], "true")
            self.assertEqual(qc_row["length_match"], "true")

    def test_prepare_assembly_rejects_non_fasta_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input.txt"
            output_path = root / "prepared.fa.gz"
            qc_path = root / "prepared.prep.tsv"
            fake_seqkit = self.make_fake_seqkit(root)
            input_path.write_text("not a fasta\n", encoding="utf-8")

            with self.assertRaises(AssemblyPreparationError):
                prepare_assembly(
                    input_path,
                    output_path,
                    qc_path,
                    wrap_width=80,
                    threads=1,
                    seqkit_executable=fake_seqkit.as_posix(),
                )


class AssemblyPrepWorkflowTests(unittest.TestCase):
    def test_busco_dry_run_routes_through_prepare_assembly(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "-R",
                "run_busco",
                "--cores",
                "4",
                "results/busco/solanum_chilense/run.complete",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )

        self.assertIn("rule prepare_assembly:", result.stdout)
        self.assertIn("work/assemblies_prepared/solanum_chilense.fa.gz", result.stdout)
        self.assertIn("--in work/assemblies_prepared/solanum_chilense.fa.gz", result.stdout)


if __name__ == "__main__":
    unittest.main()
