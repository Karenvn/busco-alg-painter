from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import pandas as pd

from busco_alg_painter.plotter import (
    calculate_alg_labels,
    canonicalize_alg_labels,
    normalize_location_columns,
    plot_locations,
)
from busco_alg_painter.profiles import load_profile


class PlotterTests(unittest.TestCase):
    def test_old_assigned_chr_column_is_supported(self) -> None:
        old = pd.DataFrame(
            {
                "buscoID": ["a", "b"],
                "query_chr": ["chr1", "chr1"],
                "position": [100, 200],
                "assigned_chr": ["m1", "mz"],
            }
        )
        normalized = normalize_location_columns(old)
        canonical = canonicalize_alg_labels(normalized, load_profile("merian"))
        self.assertEqual(canonical["assigned_alg"].tolist(), ["M1", "MZ"])

    def test_labels_follow_genomic_position(self) -> None:
        profile = load_profile("coleoptera")
        locations = pd.DataFrame(
            {
                "query_chr": ["chr1", "chr1"],
                "assigned_alg": ["CX", "C1"],
                "position": [100, 200],
            }
        )
        labels = calculate_alg_labels(locations, profile=profile, threshold=1, wrap=0)
        self.assertEqual(labels["chr1"], "CX; C1")

    def test_plot_smoke_with_auto_profile(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            location = root / "all_location.tsv"
            location.write_text(
                "buscoID\tquery_chr\tposition\tassigned_alg\tstatus\n"
                "a\tchr1\t100000\tC1\tassigned\n"
                "b\tchr1\t500000\tCX\tassigned\n"
                "c\tchr2\t250000\tC2\tassigned\n"
            )
            lengths = root / "assembly.fai"
            lengths.write_text("chr1\t1000000\nchr2\t800000\n")
            output = root / "figures" / "coleoptera"

            plot_locations(
                location_file=location,
                lengths_file=lengths,
                assembly_mode="draft",
                output_prefix=str(output),
                profile_name="auto",
                label_threshold=1,
            )
            self.assertTrue(output.with_suffix(".png").is_file())
            self.assertTrue(output.with_suffix(".svg").is_file())


if __name__ == "__main__":
    unittest.main()
