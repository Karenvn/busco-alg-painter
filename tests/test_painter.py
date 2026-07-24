from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from busco_alg_painter.painter import choose_profile, paint_buscos


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


if __name__ == "__main__":
    unittest.main()
