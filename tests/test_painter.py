from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from busco_alg_painter.painter import (
    build_location_rows,
    build_sequence_layout,
    choose_profile,
    paint_buscos,
    parse_busco_table,
    place_query_buscos,
)


def write_busco_table(path: Path, dataset: str, busco_ids: list[str]) -> None:
    lines = [
        "# BUSCO version is: 6.0.0",
        f"# The lineage dataset is: {dataset} (Creation date: test)",
        "# Busco id\tStatus\tSequence\tGene Start\tGene End",
    ]
    for index, busco_id in enumerate(busco_ids):
        start = 1000 + index * 1000
        lines.append(f"{busco_id}\tComplete\tchr1\t{start}\t{start + 100}")
    path.write_text("\n".join(lines) + "\n")


class PainterTests(unittest.TestCase):
    @staticmethod
    def sequence_report() -> list[dict]:
        return [
            {
                "role": "assembled-molecule",
                "assigned_molecule_location_type": "Chromosome",
                "chr_name": "1",
                "sequence_name": "SUPER_1",
                "genbank_accession": "CM1.1",
                "length": 1_000,
                "sort_order": 1,
            },
            {
                "role": "assembled-molecule",
                "assigned_molecule_location_type": "Chromosome",
                "chr_name": "2",
                "sequence_name": "SUPER_2",
                "genbank_accession": "CM2.1",
                "length": 2_000,
                "sort_order": 2,
            },
            {
                "role": "unlocalized-scaffold",
                "assigned_molecule_location_type": "Chromosome",
                "chr_name": "1",
                "sequence_name": "unloc_late",
                "genbank_accession": "JAAA2.1",
                "length": 100,
                "sort_order": 4,
            },
            {
                "role": "unlocalized-scaffold",
                "assigned_molecule_location_type": "Chromosome",
                "chr_name": "1",
                "sequence_name": "unloc_early",
                "genbank_accession": "JAAA1.1",
                "length": 200,
                "sort_order": 3,
            },
            {
                "role": "unplaced-scaffold",
                "chr_name": "Un",
                "sequence_name": "unplaced_1",
                "genbank_accession": "JAAAU.1",
                "length": 300,
                "sort_order": 5,
            },
        ]

    def test_coleoptera_auto_paint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            query = root / "full_table.tsv"
            write_busco_table(
                query,
                "coleoptera_odb12",
                ["100357at7041", "10053at7041", "not_in_reference"],
            )
            outputs = paint_buscos(
                query_table=query,
                prefix=root / "output",
                profile_name="auto",
                write_summary=True,
            )

            self.assertEqual(outputs.profile.id, "coleoptera")
            self.assertEqual(outputs.mapped_buscos, 2)
            text = outputs.all_locations.read_text()
            self.assertIn("\tC1\tassigned", text)
            self.assertIn("\tCX\tassigned", text)
            self.assertIn("\tNA\tunassigned", text)
            self.assertTrue(outputs.summary.is_file())

    def test_brachycera_selected_from_taxonomy_text(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            query = Path(tmp) / "full_table.tsv"
            write_busco_table(query, "diptera_odb12", ["100497at7147"])
            selected, _ = choose_profile(
                query,
                profile_name="auto",
                taxon_lineage="Eukaryota; Arthropoda; Diptera; Brachycera",
            )
            self.assertEqual(selected.id, "brachycera")

    def test_brachycera_reference_remains_an_explicit_choice(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            query = Path(tmp) / "full_table.tsv"
            write_busco_table(query, "diptera_odb12", ["100497at7147"])
            selected, _ = choose_profile(query, profile_name="brachycera")
            self.assertEqual(selected.id, "brachycera")

    def test_dataset_mismatch_is_an_error_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            query = Path(tmp) / "full_table.tsv"
            write_busco_table(query, "coleoptera_odb12", ["100357at7041"])
            with self.assertRaisesRegex(ValueError, "expects lepidoptera_odb10"):
                choose_profile(query, profile_name="merian")

    def test_missing_header_requires_explicit_profile(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            query = Path(tmp) / "full_table.tsv"
            query.write_text("100357at7041\tComplete\tchr1\t1\t100\n")
            with self.assertRaisesRegex(ValueError, "header is missing"):
                choose_profile(query, profile_name="auto")
            selected, dataset = choose_profile(query, profile_name="coleoptera")
            self.assertEqual(selected.id, "coleoptera")
            self.assertIsNone(dataset)

    def test_complete_and_every_duplicated_hit_are_retained(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            query = Path(tmp) / "full_table.tsv"
            query.write_text(
                "# The lineage dataset is: coleoptera_odb12\n"
                "# Busco id\tStatus\tSequence\tGene Start\tGene End\n"
                "100357at7041\tComplete\tSUPER_1\t100\t200\n"
                "10053at7041\tDuplicated\tunloc_early\t10\t20\n"
                "10053at7041\tDuplicated\tunloc_late\t30\t40\n"
                "100520at7041\tFragmented\tSUPER_1\t300\t400\n"
            )

            rows, _ = parse_busco_table(query)
            self.assertEqual(
                [row[0] for row in rows],
                [
                    "100357at7041",
                    "10053at7041",
                    "10053at7041",
                ],
            )

    def test_unlocalized_buscos_are_appended_and_unplaced_are_explicit(self) -> None:
        spans, placements, layout = build_sequence_layout(self.sequence_report())
        spans_by_chrom = {span.chrom: span for span in spans}
        self.assertEqual(spans_by_chrom["CM1.1"].localized_length_bp, 1_000)
        self.assertEqual(spans_by_chrom["CM1.1"].unlocalized_length_bp, 300)
        self.assertEqual(spans_by_chrom["CM1.1"].length_bp, 1_300)
        self.assertEqual(spans_by_chrom["CM2.1"].length_bp, 2_000)

        query_rows = [
            ("A", "SUPER_1", 100, 200),
            ("B", "unloc_early", 10, 20),
            ("B", "unloc_late", 30, 40),
            ("C", "unplaced_1", 50, 60),
        ]
        remapped, metadata, changed = place_query_buscos(query_rows, placements)

        self.assertEqual(changed, 4)
        self.assertEqual(remapped[0], ("A", "CM1.1", 100, 200))
        self.assertEqual(remapped[1], ("B", "CM1.1", 1_010, 1_020))
        self.assertEqual(remapped[2], ("B", "CM1.1", 1_230, 1_240))
        self.assertEqual(remapped[3], ("C", "JAAAU.1", 50, 60))
        self.assertEqual(metadata[1].position_status, "unlocalized-appended")
        self.assertEqual(metadata[3].position_status, "unplaced-excluded")
        self.assertEqual(
            [item.source_sequence for item in layout[:4]],
            ["CM1.1", "JAAA1.1", "JAAA2.1", "CM2.1"],
        )

        location_lines, mapped = build_location_rows(
            {"A": "C1", "B": "CX", "C": "C2"},
            remapped,
            placement_metadata=metadata,
        )
        self.assertEqual(mapped, 4)
        self.assertEqual(sum("\tCX\tassigned\t" in row for row in location_lines), 2)
        self.assertIn("\tunlocalized-scaffold\tunlocalized-appended", location_lines[2])
        self.assertIn("\tunplaced-scaffold\tunplaced-excluded", location_lines[4])

    def test_accession_outputs_exact_lengths_layout_and_empty_chromosome(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            query = root / "full_table.tsv"
            query.write_text(
                "# The lineage dataset is: coleoptera_odb12\n"
                "# Busco id\tStatus\tSequence\tGene Start\tGene End\n"
                "100357at7041\tComplete\tSUPER_1\t100\t200\n"
                "10053at7041\tDuplicated\tunloc_early\t10\t20\n"
                "100520at7041\tComplete\tunplaced_1\t50\t60\n"
            )
            with patch(
                "busco_alg_painter.painter.fetch_sequence_report",
                return_value=self.sequence_report(),
            ):
                outputs = paint_buscos(
                    query_table=query,
                    prefix=root / "output",
                    profile_name="coleoptera",
                    accession="GCA_TEST.1",
                    write_summary=True,
                )

            self.assertTrue(outputs.wrote_sequence_layout)
            lengths = outputs.chrom_lengths.read_text()
            self.assertIn(
                "CM1.1\t0.001300\t1300\t1000\t300\t2",
                lengths,
            )
            locations = outputs.all_locations.read_text()
            self.assertIn(
                "10053at7041\tCM1.1\t1015.0\tCX\tassigned\tJAAA1.1\t"
                "unlocalized-scaffold\tunlocalized-appended",
                locations,
            )
            self.assertIn(
                "NA\tCM2.1\tNA\tNA\tunassigned\tCM2.1\t" "assembled-molecule\tno-busco",
                locations,
            )
            self.assertIn("unplaced_1", outputs.sequence_layout.read_text())
            summary = outputs.summary.read_text()
            self.assertIn("CM1.1\t2\t1\t1", summary)
            self.assertIn("CM2.1\t0\t0\t0", summary)


if __name__ == "__main__":
    unittest.main()
