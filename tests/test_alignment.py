import csv
import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.alignment import (
    alignment_batch_output_paths,
    build_batch_ids,
    export_locus_fasta,
    export_retained_fastas,
    load_retained_locus_ids,
    locus_output_paths,
    sync_alignment_outputs,
)
from scripts.locus_matrix import parse_fasta_records
from scripts.manifest import write_tsv


REPO_ROOT = Path(__file__).resolve().parents[1]


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")


class AlignmentUnitTests(unittest.TestCase):
    def build_fixture(self):
        tempdir = TemporaryDirectory()
        repo_root = Path(tempdir.name)

        source_a = repo_root / "source" / "taxon_a.faa"
        source_b = repo_root / "source" / "taxon_b.faa"
        source_c = repo_root / "source" / "taxon_c.faa"
        write_fasta(source_a, [("orig_a", "MKTAA*")])
        write_fasta(source_b, [("orig_b", "MKT*AA*")])
        write_fasta(source_c, [("orig_c", "MKTTT*")])

        matrix_path = repo_root / "matrix.tsv"
        retained_path = repo_root / "retained.tsv"
        output_path = repo_root / "out" / "locus_keep.faa"

        write_tsv(
            matrix_path,
            [
                {
                    "locus_id": "locus_keep",
                    "sample_id": "sample_a",
                    "taxon_id": "Taxon A",
                    "sanitized_taxon_id": "taxon_a",
                    "assembly_fasta": "assemblies/taxon_a.fa.gz",
                    "status": "Complete",
                    "record_count": "1",
                    "sequence_count": "1",
                    "sequence_count_mismatch": "false",
                    "full_table_sequence_ids": "chr1",
                    "gene_starts": "1",
                    "gene_ends": "10",
                    "strands": "+",
                    "scores": "100",
                    "match_lengths_aa": "6",
                    "faa_path": source_a.relative_to(repo_root).as_posix(),
                    "gff_path": "source/taxon_a.gff",
                    "protein_lengths_aa": "5",
                    "representative_protein_length_aa": "5",
                    "has_terminal_stop_codon": "true",
                    "terminal_stop_sequence_count": "1",
                    "has_internal_stop_codon": "false",
                    "has_invalid_amino_acid": "false",
                    "invalid_amino_acids": "",
                    "is_complete_single_copy": "true",
                    "include_in_occupancy": "true",
                    "qc_notes": "",
                    "locus_total_taxa": "3",
                    "locus_complete_single_copy_taxa": "2",
                    "locus_duplicated_taxa": "0",
                    "locus_fragmented_taxa": "0",
                    "locus_missing_taxa": "1",
                    "locus_stop_codon_taxa": "1",
                    "locus_invalid_amino_acid_taxa": "0",
                    "locus_occupancy": "0.6667",
                },
                {
                    "locus_id": "locus_keep",
                    "sample_id": "sample_b",
                    "taxon_id": "Taxon B",
                    "sanitized_taxon_id": "taxon_b",
                    "assembly_fasta": "assemblies/taxon_b.fa.gz",
                    "status": "Complete",
                    "record_count": "1",
                    "sequence_count": "1",
                    "sequence_count_mismatch": "false",
                    "full_table_sequence_ids": "chr2",
                    "gene_starts": "20",
                    "gene_ends": "30",
                    "strands": "+",
                    "scores": "110",
                    "match_lengths_aa": "7",
                    "faa_path": source_b.relative_to(repo_root).as_posix(),
                    "gff_path": "source/taxon_b.gff",
                    "protein_lengths_aa": "6",
                    "representative_protein_length_aa": "6",
                    "has_terminal_stop_codon": "true",
                    "terminal_stop_sequence_count": "1",
                    "has_internal_stop_codon": "true",
                    "has_invalid_amino_acid": "false",
                    "invalid_amino_acids": "",
                    "is_complete_single_copy": "true",
                    "include_in_occupancy": "true",
                    "qc_notes": "internal_stop_codon",
                    "locus_total_taxa": "3",
                    "locus_complete_single_copy_taxa": "2",
                    "locus_duplicated_taxa": "0",
                    "locus_fragmented_taxa": "0",
                    "locus_missing_taxa": "1",
                    "locus_stop_codon_taxa": "1",
                    "locus_invalid_amino_acid_taxa": "0",
                    "locus_occupancy": "0.6667",
                },
                {
                    "locus_id": "locus_keep",
                    "sample_id": "sample_c",
                    "taxon_id": "Taxon C",
                    "sanitized_taxon_id": "taxon_c",
                    "assembly_fasta": "assemblies/taxon_c.fa.gz",
                    "status": "Missing",
                    "record_count": "0",
                    "sequence_count": "0",
                    "sequence_count_mismatch": "false",
                    "full_table_sequence_ids": "",
                    "gene_starts": "",
                    "gene_ends": "",
                    "strands": "",
                    "scores": "",
                    "match_lengths_aa": "",
                    "faa_path": source_c.relative_to(repo_root).as_posix(),
                    "gff_path": "source/taxon_c.gff",
                    "protein_lengths_aa": "",
                    "representative_protein_length_aa": "",
                    "has_terminal_stop_codon": "false",
                    "terminal_stop_sequence_count": "0",
                    "has_internal_stop_codon": "false",
                    "has_invalid_amino_acid": "false",
                    "invalid_amino_acids": "",
                    "is_complete_single_copy": "false",
                    "include_in_occupancy": "false",
                    "qc_notes": "missing_status",
                    "locus_total_taxa": "3",
                    "locus_complete_single_copy_taxa": "2",
                    "locus_duplicated_taxa": "0",
                    "locus_fragmented_taxa": "0",
                    "locus_missing_taxa": "1",
                    "locus_stop_codon_taxa": "1",
                    "locus_invalid_amino_acid_taxa": "0",
                    "locus_occupancy": "0.6667",
                },
            ],
            [
                "locus_id",
                "sample_id",
                "taxon_id",
                "sanitized_taxon_id",
                "assembly_fasta",
                "status",
                "record_count",
                "sequence_count",
                "sequence_count_mismatch",
                "full_table_sequence_ids",
                "gene_starts",
                "gene_ends",
                "strands",
                "scores",
                "match_lengths_aa",
                "faa_path",
                "gff_path",
                "protein_lengths_aa",
                "representative_protein_length_aa",
                "has_terminal_stop_codon",
                "terminal_stop_sequence_count",
                "has_internal_stop_codon",
                "has_invalid_amino_acid",
                "invalid_amino_acids",
                "is_complete_single_copy",
                "include_in_occupancy",
                "qc_notes",
                "locus_total_taxa",
                "locus_complete_single_copy_taxa",
                "locus_duplicated_taxa",
                "locus_fragmented_taxa",
                "locus_missing_taxa",
                "locus_stop_codon_taxa",
                "locus_invalid_amino_acid_taxa",
                "locus_occupancy",
            ],
        )

        write_tsv(
            retained_path,
            [
                {
                    "locus_id": "locus_keep",
                    "total_taxa": "3",
                    "complete_single_copy_taxa": "2",
                    "duplicated_taxa": "0",
                    "fragmented_taxa": "0",
                    "missing_taxa": "1",
                    "stop_codon_taxa": "1",
                    "invalid_amino_acid_taxa": "0",
                    "occupancy": "0.6667",
                    "occupancy_threshold": "0.6000",
                    "passes_occupancy": "true",
                    "decision": "retain",
                    "failure_reasons": "",
                    "qc_warnings": "internal_stop_codon_present",
                    "retained_sample_ids": "sample_a,sample_b",
                    "retained_sanitized_taxon_ids": "taxon_a,taxon_b",
                    "retained_protein_lengths_aa": "5,6",
                    "min_protein_length_aa": "5",
                    "median_protein_length_aa": "5.5000",
                    "max_protein_length_aa": "6",
                    "length_span_aa": "1",
                    "length_ratio": "1.2000",
                    "length_dispersion_observed": "true",
                }
            ],
            [
                "locus_id",
                "total_taxa",
                "complete_single_copy_taxa",
                "duplicated_taxa",
                "fragmented_taxa",
                "missing_taxa",
                "stop_codon_taxa",
                "invalid_amino_acid_taxa",
                "occupancy",
                "occupancy_threshold",
                "passes_occupancy",
                "decision",
                "failure_reasons",
                "qc_warnings",
                "retained_sample_ids",
                "retained_sanitized_taxon_ids",
                "retained_protein_lengths_aa",
                "min_protein_length_aa",
                "median_protein_length_aa",
                "max_protein_length_aa",
                "length_span_aa",
                "length_ratio",
                "length_dispersion_observed",
            ],
        )
        return tempdir, repo_root, matrix_path, retained_path, output_path

    def test_export_locus_fasta_includes_only_retained_taxa(self):
        tempdir, repo_root, matrix_path, retained_path, output_path = self.build_fixture()
        self.addCleanup(tempdir.cleanup)

        export_locus_fasta("locus_keep", matrix_path, retained_path, output_path, repo_root)
        records = parse_fasta_records(output_path)

        self.assertEqual([header for header, _ in records], ["taxon_a", "taxon_b"])
        self.assertEqual(len(records), 2)

    def test_exported_headers_use_sanitized_taxon_ids(self):
        tempdir, repo_root, matrix_path, retained_path, output_path = self.build_fixture()
        self.addCleanup(tempdir.cleanup)

        export_locus_fasta("locus_keep", matrix_path, retained_path, output_path, repo_root)
        records = parse_fasta_records(output_path)

        self.assertEqual(records[0][0], "taxon_a")
        self.assertEqual(records[1][0], "taxon_b")
        self.assertEqual(records[1][1], "MKT*AA*")

    def test_export_retained_fastas_writes_manifest_and_cleans_stale_files(self):
        tempdir, repo_root, matrix_path, retained_path, output_path = self.build_fixture()
        self.addCleanup(tempdir.cleanup)

        output_dir = output_path.parent
        manifest_path = repo_root / "out" / "raw_fastas_manifest.tsv"
        stale_path = output_dir / "stale.faa"
        stale_path.parent.mkdir(parents=True, exist_ok=True)
        stale_path.write_text(">stale\nM\n", encoding="utf-8")

        export_retained_fastas(matrix_path, retained_path, output_dir, manifest_path, repo_root)

        records = parse_fasta_records(output_dir / "locus_keep.faa")
        self.assertEqual([header for header, _ in records], ["taxon_a", "taxon_b"])
        self.assertFalse(stale_path.exists())

        with manifest_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            rows,
            [
                {
                    "locus_id": "locus_keep",
                    "raw_fasta": (output_dir / "locus_keep.faa").as_posix(),
                    "taxon_count": "2",
                }
            ],
        )

    def test_locus_output_paths_are_deterministic(self):
        paths = locus_output_paths("35at4069")
        self.assertEqual(paths["raw_fasta"], "results/loci/raw_fastas/35at4069.faa")
        self.assertEqual(paths["alignment"], "results/loci/alignments/35at4069.aln.faa")

    def test_alignment_batch_helpers_are_deterministic(self):
        self.assertEqual(build_batch_ids(["1at1", "2at1", "3at1", "4at1", "5at1"], 2), ["0000", "0001", "0002"])
        self.assertEqual(
            alignment_batch_output_paths("0003")["completion"],
            "results/loci/batches/alignment_batch_0003.complete",
        )

    def test_sync_alignment_outputs_removes_stale_alignments_and_logs(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            manifest_path = tmp_path / "raw_fastas_manifest.tsv"
            alignment_dir = tmp_path / "alignments"
            log_dir = tmp_path / "logs"
            alignment_dir.mkdir()
            log_dir.mkdir()

            write_tsv(
                manifest_path,
                [{"locus_id": "locus_keep", "raw_fasta": "raw/locus_keep.faa", "taxon_count": "2"}],
                ["locus_id", "raw_fasta", "taxon_count"],
            )
            (alignment_dir / "locus_keep.aln.faa").write_text(">a\nM\n", encoding="utf-8")
            (alignment_dir / "stale.aln.faa").write_text(">a\nM\n", encoding="utf-8")
            (log_dir / "locus_keep.log").write_text("ok\n", encoding="utf-8")
            (log_dir / "stale.log").write_text("stale\n", encoding="utf-8")

            sync_alignment_outputs(manifest_path, alignment_dir, log_dir)

            self.assertTrue((alignment_dir / "locus_keep.aln.faa").is_file())
            self.assertFalse((alignment_dir / "stale.aln.faa").exists())
            self.assertTrue((log_dir / "locus_keep.log").is_file())
            self.assertFalse((log_dir / "stale.log").exists())


class AlignmentSmokeTests(unittest.TestCase):
    def test_mafft_smoke_alignment_preserves_internal_stop_with_anysymbol(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            input_path = tmp_path / "input.faa"
            output_path = tmp_path / "output.faa"
            input_path.write_text(">a\nMA*AA*\n>b\nMAAAA*\n", encoding="utf-8")

            result = subprocess.run(
                [
                    "mafft",
                    "--amino",
                    "--anysymbol",
                    "--auto",
                    str(input_path),
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            output_path.write_text(result.stdout, encoding="utf-8")
            records = parse_fasta_records(output_path)

        self.assertEqual(records[0][1], "MA*AA*")
        self.assertEqual(records[1][1], "MAAAA*")


class AlignmentWorkflowTests(unittest.TestCase):
    def test_alignment_rule_dry_run_renders_batched_execution(self):
        result = subprocess.run(
            [
                "snakemake",
                "-n",
                "-p",
                "-R",
                "export_retained_fastas",
                "align_locus_batch",
                "--cores",
                "4",
                "results/loci/batches/alignment_batch_0000.complete",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("rule export_retained_fastas:", result.stdout)
        self.assertIn("rule align_locus_batch:", result.stdout)
        self.assertIn("python3 -m scripts.run_alignment_batch", result.stdout)

    def test_exported_fasta_matches_retained_loci_output_for_real_locus(self):
        raw_fasta_path = REPO_ROOT / "results/loci/raw_fastas/35at4069.faa"
        manifest_path = REPO_ROOT / "results/loci/raw_fastas_manifest.tsv"
        if not raw_fasta_path.is_file() or not manifest_path.is_file():
            self.skipTest("Real retained-locus FASTA outputs are not present.")

        records = parse_fasta_records(raw_fasta_path)
        headers = [header for header, _ in records]

        with (REPO_ROOT / "results/qc/retained_loci.tsv").open(newline="", encoding="utf-8") as handle:
            retained_rows = list(csv.DictReader(handle, delimiter="\t"))
        retained_map = {row["locus_id"]: row for row in retained_rows}

        self.assertEqual(headers, retained_map["35at4069"]["retained_sanitized_taxon_ids"].split(","))
        self.assertEqual(headers, sorted(headers))

        with manifest_path.open(newline="", encoding="utf-8") as handle:
            manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
        manifest_map = {row["locus_id"]: row for row in manifest_rows}
        self.assertEqual(
            manifest_map["35at4069"]["raw_fasta"],
            "results/loci/raw_fastas/35at4069.faa",
        )

    def test_retained_locus_id_loader_matches_real_table(self):
        with (REPO_ROOT / "results/qc/retained_loci.tsv").open(newline="", encoding="utf-8") as handle:
            retained_rows = list(csv.DictReader(handle, delimiter="\t"))
        retained_ids = [row["locus_id"] for row in retained_rows if row["decision"] == "retain"]
        self.assertEqual(
            load_retained_locus_ids(REPO_ROOT / "results/qc/retained_loci.tsv"),
            sorted(retained_ids, key=lambda value: (int(value.split("at", 1)[0]), value.split("at", 1)[1])),
        )


if __name__ == "__main__":
    unittest.main()
