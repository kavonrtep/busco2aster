import csv
import os
import subprocess
import textwrap
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.dna_extract import (
    build_gffread_command,
    extract_retained_sample_dna,
    per_locus_dna_sequence_path,
    sample_dna_output_paths,
)
from scripts.locus_matrix import parse_fasta_records
from scripts.manifest import write_tsv
from scripts.sequence_mode import resolve_scfl_model


REPO_ROOT = Path(__file__).resolve().parents[1]


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


class DnaExtractUnitTests(unittest.TestCase):
    def test_resolve_scfl_model_maps_legacy_protein_default_to_dna_default(self):
        self.assertEqual(resolve_scfl_model("dna", "LG+G4"), "GTR+G")
        self.assertEqual(resolve_scfl_model("dna", "auto"), "GTR+G")
        self.assertEqual(resolve_scfl_model("dna", "HKY+G"), "HKY+G")

    def test_build_gffread_command_renders_expected_flags(self):
        command = build_gffread_command(
            gffread_executable="gffread",
            annotation_path=Path("results/dna_sequences/sample/retained_loci.gff3"),
            genome_fasta_path=Path("work/assemblies_prepared/sample.fa"),
            output_fasta_path=Path("results/dna_sequences/sample/retained_loci.fna"),
        )
        self.assertEqual(command[0], "gffread")
        self.assertIn("-g", command)
        self.assertIn("-x", command)
        self.assertIn("results/dna_sequences/sample/retained_loci.gff3", command)
        self.assertIn("work/assemblies_prepared/sample.fa", command)
        self.assertIn("results/dna_sequences/sample/retained_loci.fna", command)

    def test_extract_retained_sample_dna_writes_per_locus_fastas(self):
        with TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            matrix_path = repo_root / "matrix.tsv"
            retained_path = repo_root / "retained.tsv"
            genome_path = repo_root / "assemblies" / "sample1.fa"
            gff_dir = repo_root / "gff"
            fake_gffread = repo_root / "fake_gffread.py"

            write_text(
                genome_path,
                ">chr1\nATGGCCGAA\n>chr2\nTCATTTCAT\n",
            )
            write_text(
                gff_dir / "locus_plus.gff",
                "\n".join(
                    [
                        "chr1\tminiprot\tmRNA\t1\t9\t.\t+\t.\tID=MP1;Target=locus_plus_1 1 3",
                        "chr1\tminiprot\tCDS\t1\t3\t.\t+\t0\tParent=MP1;Target=locus_plus_1 1 1",
                        "chr1\tminiprot\tCDS\t4\t9\t.\t+\t0\tParent=MP1;Target=locus_plus_1 2 3",
                        "",
                    ]
                ),
            )
            write_text(
                gff_dir / "locus_minus.gff",
                "\n".join(
                    [
                        "chr2\tminiprot\tmRNA\t1\t9\t.\t-\t.\tID=MP2;Target=locus_minus_1 1 3",
                        "chr2\tminiprot\tCDS\t1\t3\t.\t-\t0\tParent=MP2;Target=locus_minus_1 3 3",
                        "chr2\tminiprot\tCDS\t4\t9\t.\t-\t0\tParent=MP2;Target=locus_minus_1 1 2",
                        "",
                    ]
                ),
            )

            write_text(
                fake_gffread,
                textwrap.dedent(
                    """\
                    #!/usr/bin/env python3
                    import sys
                    from pathlib import Path

                    def read_fasta(path):
                        records = {}
                        header = None
                        parts = []
                        for raw in Path(path).read_text(encoding="utf-8").splitlines():
                            line = raw.strip()
                            if not line:
                                continue
                            if line.startswith(">"):
                                if header is not None:
                                    records[header] = "".join(parts)
                                header = line[1:]
                                parts = []
                            else:
                                parts.append(line)
                        if header is not None:
                            records[header] = "".join(parts)
                        return records

                    def revcomp(seq):
                        table = str.maketrans("ACGTacgt", "TGCAtgca")
                        return seq.translate(table)[::-1]

                    args = sys.argv[1:]
                    gff_path = Path(args[0])
                    genome_path = Path(args[args.index("-g") + 1])
                    output_path = Path(args[args.index("-x") + 1])
                    genome = read_fasta(genome_path)
                    rows = {}
                    strands = {}
                    with gff_path.open(encoding="utf-8") as handle:
                        for raw in handle:
                            line = raw.strip()
                            if not line or line.startswith("#"):
                                continue
                            cols = line.split("\\t")
                            attrs = {}
                            for item in cols[8].split(";"):
                                if not item:
                                    continue
                                key, value = item.split("=", 1)
                                attrs[key] = value
                            if cols[2] == "mRNA":
                                strands[attrs["ID"]] = cols[6]
                            elif cols[2] == "CDS":
                                rows.setdefault(attrs["Parent"], []).append((cols[0], int(cols[3]), int(cols[4])))
                    with output_path.open("w", encoding="utf-8") as out:
                        for transcript_id, parts in rows.items():
                            strand = strands[transcript_id]
                            ordered = sorted(parts, key=lambda item: item[1])
                            seq = "".join(genome[chrom][start - 1:end] for chrom, start, end in ordered)
                            if strand == "-":
                                seq = revcomp(seq)
                            out.write(f">{transcript_id}\\n{seq}\\n")
                    """
                ),
            )
            fake_gffread.chmod(0o755)

            write_tsv(
                matrix_path,
                [
                    {
                        "locus_id": "locus_plus",
                        "sample_id": "sample1",
                        "taxon_id": "Taxon one",
                        "sanitized_taxon_id": "taxon_one",
                        "assembly_fasta": genome_path.as_posix(),
                        "status": "Complete",
                        "record_count": "1",
                        "sequence_count": "1",
                        "sequence_count_mismatch": "false",
                        "full_table_sequence_ids": "chr1",
                        "gene_starts": "1",
                        "gene_ends": "9",
                        "strands": "+",
                        "scores": "100",
                        "match_lengths_aa": "3",
                        "faa_path": "unused/locus_plus.faa",
                        "gff_path": (gff_dir / "locus_plus.gff").relative_to(repo_root).as_posix(),
                        "protein_lengths_aa": "3",
                        "representative_protein_length_aa": "3",
                        "has_terminal_stop_codon": "false",
                        "terminal_stop_sequence_count": "0",
                        "has_internal_stop_codon": "false",
                        "has_invalid_amino_acid": "false",
                        "invalid_amino_acids": "",
                        "is_complete_single_copy": "true",
                        "include_in_occupancy": "true",
                        "qc_notes": "",
                        "locus_total_taxa": "1",
                        "locus_complete_single_copy_taxa": "1",
                        "locus_duplicated_taxa": "0",
                        "locus_fragmented_taxa": "0",
                        "locus_missing_taxa": "0",
                        "locus_stop_codon_taxa": "0",
                        "locus_invalid_amino_acid_taxa": "0",
                        "locus_occupancy": "1.0000",
                    },
                    {
                        "locus_id": "locus_minus",
                        "sample_id": "sample1",
                        "taxon_id": "Taxon one",
                        "sanitized_taxon_id": "taxon_one",
                        "assembly_fasta": genome_path.as_posix(),
                        "status": "Complete",
                        "record_count": "1",
                        "sequence_count": "1",
                        "sequence_count_mismatch": "false",
                        "full_table_sequence_ids": "chr2",
                        "gene_starts": "1",
                        "gene_ends": "9",
                        "strands": "-",
                        "scores": "100",
                        "match_lengths_aa": "3",
                        "faa_path": "unused/locus_minus.faa",
                        "gff_path": (gff_dir / "locus_minus.gff").relative_to(repo_root).as_posix(),
                        "protein_lengths_aa": "3",
                        "representative_protein_length_aa": "3",
                        "has_terminal_stop_codon": "false",
                        "terminal_stop_sequence_count": "0",
                        "has_internal_stop_codon": "false",
                        "has_invalid_amino_acid": "false",
                        "invalid_amino_acids": "",
                        "is_complete_single_copy": "true",
                        "include_in_occupancy": "true",
                        "qc_notes": "",
                        "locus_total_taxa": "1",
                        "locus_complete_single_copy_taxa": "1",
                        "locus_duplicated_taxa": "0",
                        "locus_fragmented_taxa": "0",
                        "locus_missing_taxa": "0",
                        "locus_stop_codon_taxa": "0",
                        "locus_invalid_amino_acid_taxa": "0",
                        "locus_occupancy": "1.0000",
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
                        "locus_id": "locus_plus",
                        "decision": "retain",
                    },
                    {
                        "locus_id": "locus_minus",
                        "decision": "retain",
                    },
                ],
                ["locus_id", "decision"],
            )

            cwd = Path.cwd()
            try:
                os.chdir(repo_root)
                rows = extract_retained_sample_dna(
                    matrix_path=matrix_path,
                    retained_path=retained_path,
                    sample_id="sample1",
                    genome_fasta_path=genome_path,
                    repo_root=repo_root,
                    gffread_executable=fake_gffread.as_posix(),
                )
            finally:
                os.chdir(cwd)

            self.assertEqual([row["locus_id"] for row in rows], ["locus_minus", "locus_plus"])
            plus_records = parse_fasta_records(repo_root / per_locus_dna_sequence_path("sample1", "locus_plus"))
            minus_records = parse_fasta_records(repo_root / per_locus_dna_sequence_path("sample1", "locus_minus"))
            self.assertEqual(plus_records, [("taxon_one", "ATGGCCGAA")])
            self.assertEqual(minus_records, [("taxon_one", "ATGAAATGA")])
            records_path = repo_root / sample_dna_output_paths("sample1")["records"]
            with records_path.open(newline="", encoding="utf-8") as handle:
                record_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(record_rows), 2)


class DnaModeWorkflowTests(unittest.TestCase):
    def test_dna_mode_dry_run_wires_gffread_and_dna_tree_settings(self):
        with TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "dna_config.yaml"
            config_text = (REPO_ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
            config_path.write_text(
                config_text.replace("sequence_type: protein", "sequence_type: dna"),
                encoding="utf-8",
            )
            result = subprocess.run(
                [
                    "snakemake",
                    "-n",
                    "-p",
                    "-R",
                    "extract_retained_dna_sample",
                    "infer_gene_trees",
                    "infer_site_concordance",
                    "--rerun-incomplete",
                    "--configfile",
                    config_path.as_posix(),
                    "--cores",
                    "4",
                    "results/report/report.html",
                ],
                cwd=REPO_ROOT,
                capture_output=True,
                text=True,
                check=True,
            )

        self.assertTrue(
            "rule prepare_assembly_plain:" in result.stdout
            or "work/assemblies_prepared/solanum_chilense.fa.fai" in result.stdout
        )
        self.assertIn("rule extract_retained_dna_sample:", result.stdout)
        self.assertIn("rule infer_gene_trees:", result.stdout)
        self.assertIn("rule infer_site_concordance:", result.stdout)
        self.assertIn("retained_loci.fna", result.stdout)
        self.assertIn("results/dna_sequences/solanum_chilense/run.complete", result.stdout)


if __name__ == "__main__":
    unittest.main()
