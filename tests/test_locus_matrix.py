import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.locus_matrix import build_locus_taxon_matrix_rows, build_retained_loci_rows


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")


def write_gff(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("##gff-version 3\n", encoding="utf-8")


def make_summary_rows() -> list[dict[str, str]]:
    return [
        {
            "sample_id": sample_id,
            "taxon_id": sample_id.replace("_", " ").title(),
            "sanitized_taxon_id": sample_id,
            "assembly_fasta": f"assemblies/{sample_id}.fa.gz",
        }
        for sample_id in ("sample_a", "sample_b", "sample_c", "sample_d")
    ]


def make_record(
    sample_id: str,
    locus_id: str,
    status: str,
    faa_path: str = "",
    gff_path: str = "",
    sequence_id: str = "",
    gene_start: str = "",
    gene_end: str = "",
    strand: str = "",
    score: str = "",
    length: str = "",
) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "taxon_id": sample_id.replace("_", " ").title(),
        "sanitized_taxon_id": sample_id,
        "assembly_fasta": f"assemblies/{sample_id}.fa.gz",
        "busco_id": locus_id,
        "status": status,
        "sequence_id": sequence_id,
        "gene_start": gene_start,
        "gene_end": gene_end,
        "strand": strand,
        "score": score,
        "length": length,
        "orthodb_url": f"https://example.org/{locus_id}",
        "description": f"{locus_id} description",
        "sequence_category": "",
        "sequence_dir": "",
        "faa_path": faa_path,
        "gff_path": gff_path,
    }


class LocusMatrixTests(unittest.TestCase):
    def build_fixture(self):
        tempdir = TemporaryDirectory()
        repo_root = Path(tempdir.name)
        summary_rows = make_summary_rows()
        record_rows: list[dict[str, str]] = []

        def add_complete(locus_id: str, sample_id: str, sequence: str, index: int) -> None:
            faa_path = f"fixtures/{locus_id}/{sample_id}.faa"
            gff_path = f"fixtures/{locus_id}/{sample_id}.gff"
            write_fasta(repo_root / faa_path, [(f"{locus_id}_{sample_id}", sequence)])
            write_gff(repo_root / gff_path)
            record_rows.append(
                make_record(
                    sample_id=sample_id,
                    locus_id=locus_id,
                    status="Complete",
                    faa_path=faa_path,
                    gff_path=gff_path,
                    sequence_id=f"chr{index}",
                    gene_start=str(index * 10),
                    gene_end=str(index * 10 + 9),
                    strand="+",
                    score=str(100 + index),
                    length=str(len(sequence.replace('*', ''))),
                )
            )

        for idx, (sample_id, sequence) in enumerate(
            {
                "sample_a": "M" * 10,
                "sample_b": "M" * 12,
                "sample_c": "M" * 11,
                "sample_d": "M" * 13,
            }.items(),
            start=1,
        ):
            add_complete("locus_keep", sample_id, sequence, idx)

        for idx, (sample_id, sequence) in enumerate(
            {
                "sample_a": "A" * 8,
                "sample_b": "A" * 8,
                "sample_c": "A" * 9,
            }.items(),
            start=1,
        ):
            add_complete("locus_missing", sample_id, sequence, idx)
        record_rows.append(make_record(sample_id="sample_d", locus_id="locus_missing", status="Missing"))

        add_complete("locus_dup_frag", "sample_a", "C" * 7, 1)
        dup_faa = "fixtures/locus_dup_frag/sample_b.faa"
        dup_gff = "fixtures/locus_dup_frag/sample_b.gff"
        write_fasta(
            repo_root / dup_faa,
            [("dup1", "CCCCCCC"), ("dup2", "CCCCCCC")],
        )
        write_gff(repo_root / dup_gff)
        record_rows.extend(
            [
                make_record(
                    sample_id="sample_b",
                    locus_id="locus_dup_frag",
                    status="Duplicated",
                    faa_path=dup_faa,
                    gff_path=dup_gff,
                    sequence_id="chr2",
                    gene_start="20",
                    gene_end="29",
                    strand="+",
                    score="210",
                    length="7",
                ),
                make_record(
                    sample_id="sample_b",
                    locus_id="locus_dup_frag",
                    status="Duplicated",
                    faa_path=dup_faa,
                    gff_path=dup_gff,
                    sequence_id="chr3",
                    gene_start="30",
                    gene_end="39",
                    strand="-",
                    score="211",
                    length="7",
                ),
            ]
        )
        frag_faa = "fixtures/locus_dup_frag/sample_c.faa"
        frag_gff = "fixtures/locus_dup_frag/sample_c.gff"
        write_fasta(repo_root / frag_faa, [("frag1", "CCCCCCC")])
        write_gff(repo_root / frag_gff)
        record_rows.append(
            make_record(
                sample_id="sample_c",
                locus_id="locus_dup_frag",
                status="Fragmented",
                faa_path=frag_faa,
                gff_path=frag_gff,
                sequence_id="chr4",
                gene_start="40",
                gene_end="49",
                strand="+",
                score="120",
                length="7",
            )
        )
        add_complete("locus_dup_frag", "sample_d", "C" * 7, 4)

        add_complete("locus_qc", "sample_a", "MA*AA*", 1)
        add_complete("locus_qc", "sample_b", "MA?AA*", 2)
        add_complete("locus_qc", "sample_c", "MAAAA*", 3)
        add_complete("locus_qc", "sample_d", "MAAAA*", 4)

        matrix_rows = build_locus_taxon_matrix_rows(summary_rows, record_rows, repo_root)
        retained_rows = build_retained_loci_rows(matrix_rows, 0.8)
        return tempdir, matrix_rows, retained_rows

    def test_occupancy_is_computed_correctly(self):
        tempdir, matrix_rows, _ = self.build_fixture()
        self.addCleanup(tempdir.cleanup)
        locus_rows = [row for row in matrix_rows if row["locus_id"] == "locus_missing"]
        self.assertEqual({row["locus_occupancy"] for row in locus_rows}, {"0.7500"})
        self.assertEqual({row["locus_complete_single_copy_taxa"] for row in locus_rows}, {3})

    def test_duplicated_and_fragmented_hits_are_excluded_from_retained_loci(self):
        tempdir, _, retained_rows = self.build_fixture()
        self.addCleanup(tempdir.cleanup)
        retained = {row["locus_id"]: row for row in retained_rows}
        locus = retained["locus_dup_frag"]
        self.assertEqual(locus["complete_single_copy_taxa"], 2)
        self.assertEqual(locus["duplicated_taxa"], 1)
        self.assertEqual(locus["fragmented_taxa"], 1)
        self.assertEqual(locus["retained_sample_ids"], "sample_a,sample_d")
        self.assertEqual(locus["decision"], "exclude")

    def test_stop_codon_and_invalid_character_flags_propagate(self):
        tempdir, matrix_rows, _ = self.build_fixture()
        self.addCleanup(tempdir.cleanup)
        cell_map = {(row["locus_id"], row["sample_id"]): row for row in matrix_rows}
        stop_row = cell_map[("locus_qc", "sample_a")]
        invalid_row = cell_map[("locus_qc", "sample_b")]

        self.assertEqual(stop_row["has_internal_stop_codon"], "true")
        self.assertIn("internal_stop_codon", stop_row["qc_notes"])
        self.assertEqual(invalid_row["has_invalid_amino_acid"], "true")
        self.assertEqual(invalid_row["invalid_amino_acids"], "?")

    def test_loci_below_threshold_are_excluded(self):
        tempdir, _, retained_rows = self.build_fixture()
        self.addCleanup(tempdir.cleanup)
        retained = {row["locus_id"]: row for row in retained_rows}
        self.assertEqual(retained["locus_missing"]["occupancy"], "0.7500")
        self.assertEqual(retained["locus_missing"]["passes_occupancy"], "false")
        self.assertEqual(retained["locus_missing"]["decision"], "exclude")
        self.assertIn("occupancy_below_threshold", retained["locus_missing"]["failure_reasons"])

    def test_retained_loci_table_is_consistent_with_matrix(self):
        tempdir, matrix_rows, retained_rows = self.build_fixture()
        self.addCleanup(tempdir.cleanup)
        retained_map = {row["locus_id"]: row for row in retained_rows}

        for locus_id in {row["locus_id"] for row in matrix_rows}:
            locus_rows = [row for row in matrix_rows if row["locus_id"] == locus_id]
            retained_row = retained_map[locus_id]
            complete_count = sum(row["include_in_occupancy"] == "true" for row in locus_rows)
            duplicated_count = sum(row["status"] == "Duplicated" for row in locus_rows)
            fragmented_count = sum(row["status"] == "Fragmented" for row in locus_rows)
            missing_count = sum(row["status"] in {"Missing", "NotReported"} for row in locus_rows)

            self.assertEqual(int(retained_row["complete_single_copy_taxa"]), complete_count)
            self.assertEqual(int(retained_row["duplicated_taxa"]), duplicated_count)
            self.assertEqual(int(retained_row["fragmented_taxa"]), fragmented_count)
            self.assertEqual(int(retained_row["missing_taxa"]), missing_count)

    def test_length_dispersion_is_reported_but_not_filtered(self):
        tempdir, _, retained_rows = self.build_fixture()
        self.addCleanup(tempdir.cleanup)
        retained = {row["locus_id"]: row for row in retained_rows}
        locus = retained["locus_keep"]
        self.assertEqual(locus["length_dispersion_observed"], "true")
        self.assertIn("length_dispersion_observed", locus["qc_warnings"])
        self.assertEqual(locus["decision"], "retain")


if __name__ == "__main__":
    unittest.main()
