from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("MPLBACKEND", "Agg")

import pandas as pd

from busco_alg_painter.plotter import (
    calculate_alg_labels,
    calculate_windowed_alg_labels,
    canonicalize_alg_labels,
    load_data,
    load_lengths,
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
            self.assertIn("<!-- chr2 -->", output.with_suffix(".svg").read_text())

    def test_bar_alg_labels_can_be_hidden(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            location = root / "all_location.tsv"
            location.write_text(
                "buscoID\tquery_chr\tposition\tassigned_alg\tstatus\n"
                "a\tchr1\t100000\tC1\tassigned\n"
            )
            lengths = root / "assembly.fai"
            lengths.write_text("chr1\t1000000\n")
            output = root / "figures" / "unlabelled"

            with (
                patch(
                    "busco_alg_painter.plotter.calculate_alg_labels"
                ) as chromosome_labels,
                patch(
                    "busco_alg_painter.plotter.calculate_windowed_alg_labels"
                ) as window_labels,
            ):
                plot_locations(
                    location_file=location,
                    lengths_file=lengths,
                    assembly_mode="draft",
                    output_prefix=str(output),
                    profile_name="coleoptera",
                    show_bar_alg_labels=False,
                    label_window_mb=10,
                )

            chromosome_labels.assert_not_called()
            window_labels.assert_not_called()
            self.assertTrue(output.with_suffix(".png").is_file())
            self.assertTrue(output.with_suffix(".svg").is_file())

    def test_bar_alg_labels_default_to_merian_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            lengths = root / "assembly.fai"
            lengths.write_text("chr1\t1000000\n")

            cases = (
                ("M1", "merian", 1),
                ("db1a", "brachycera", 0),
                ("d1", "diptera", 0),
                ("C1", "coleoptera", 0),
            )
            for label, expected_profile, expected_calls in cases:
                with self.subTest(profile=expected_profile):
                    location = root / f"{expected_profile}.tsv"
                    location.write_text(
                        "buscoID\tquery_chr\tposition\tassigned_alg\tstatus\n"
                        f"a\tchr1\t100000\t{label}\tassigned\n"
                    )

                    with patch(
                        "busco_alg_painter.plotter.calculate_alg_labels",
                        wraps=calculate_alg_labels,
                    ) as chromosome_labels:
                        plot_locations(
                            location_file=location,
                            lengths_file=lengths,
                            assembly_mode="draft",
                            output_prefix=str(root / expected_profile),
                            profile_name="auto",
                            label_threshold=1,
                        )

                    self.assertEqual(chromosome_labels.call_count, expected_calls)

    def test_numeric_sequence_names_are_normalized_as_strings(self) -> None:
        locations = pd.DataFrame(
            {
                "buscoID": ["a", "b"],
                "query_chr": [1, 2],
                "position": [100, 200],
                "assigned_alg": ["C1", "C2"],
            }
        )

        normalized = normalize_location_columns(locations)

        self.assertEqual(normalized["query_chr"].tolist(), ["1", "2"])

    def test_final_lengths_prefer_exact_basepair_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            lengths_file = Path(tmp) / "chrom_lengths.tsv"
            lengths_file.write_text(
                "Chrom\tLength_Mb\tLength_bp\tLocalized_Length_bp\t"
                "Unlocalized_Length_bp\tUnlocalized_scaffold_count\n"
                "1\t1.235\t1234567\t1200000\t34567\t2\n"
            )

            lengths = load_lengths(lengths_file, assembly_mode="final")

            self.assertEqual(lengths.loc[0, "query_chr"], "1")
            self.assertEqual(lengths.loc[0, "length"], 1_234_567)
            self.assertEqual(lengths.loc[0, "localized_length"], 1_200_000)
            self.assertEqual(lengths.loc[0, "unlocalized_length"], 34_567)

    def test_lengths_are_required_unless_estimation_is_explicit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            location = Path(tmp) / "all_location.tsv"
            location.write_text(
                "buscoID\tquery_chr\tposition\tassigned_alg\tstatus\n"
                "a\tchr1\t100000\tC1\tassigned\n"
            )

            with self.assertRaisesRegex(ValueError, "lengths file is required"):
                load_data(location)
            _, estimated = load_data(location, allow_estimated_lengths=True)
            self.assertEqual(estimated.loc[0, "length"], 105_000)

    def test_missing_length_for_assigned_busco_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            location = root / "all_location.tsv"
            location.write_text(
                "buscoID\tquery_chr\tposition\tassigned_alg\tstatus\n"
                "a\tmissing_chr\t100000\tC1\tassigned\n"
            )
            lengths = root / "assembly.fai"
            lengths.write_text("chr1\t1000000\n")

            with self.assertRaisesRegex(ValueError, "no matching.*length"):
                plot_locations(
                    location_file=location,
                    lengths_file=lengths,
                    assembly_mode="draft",
                    output_prefix=str(root / "missing"),
                    profile_name="coleoptera",
                )

    def test_unplaced_busco_is_reported_but_not_treated_as_join_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            location = root / "all_location.tsv"
            location.write_text(
                "buscoID\tquery_chr\tposition\tassigned_alg\tstatus\t"
                "source_sequence\tsequence_role\tposition_status\n"
                "a\tchr1\t100000\tC1\tassigned\tchr1\t"
                "assembled-molecule\tchromosome\n"
                "b\tunplaced1\t100\tC2\tassigned\tunplaced1\t"
                "unplaced-scaffold\tunplaced-excluded\n"
            )
            lengths = root / "assembly.fai"
            lengths.write_text("chr1\t1000000\n")
            output = root / "unplaced"

            plot_locations(
                location_file=location,
                lengths_file=lengths,
                assembly_mode="draft",
                output_prefix=str(output),
                profile_name="coleoptera",
            )

            self.assertTrue(output.with_suffix(".svg").is_file())

    def test_unlocalized_tail_and_busco_render_in_explicit_appended_region(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            location = root / "all_location.tsv"
            location.write_text(
                "buscoID\tquery_chr\tposition\tassigned_alg\tstatus\t"
                "source_sequence\tsequence_role\tposition_status\n"
                "a\tCM1.1\t1100000\tC1\tassigned\tJAAA1.1\t"
                "unlocalized-scaffold\tunlocalized-appended\n"
                "NA\tCM2.1\tNA\tNA\tunassigned\tCM2.1\t"
                "assembled-molecule\tno-busco\n"
            )
            lengths = root / "chrom_lengths.tsv"
            lengths.write_text(
                "Chrom\tLength_Mb\tLength_bp\tLocalized_Length_bp\t"
                "Unlocalized_Length_bp\tUnlocalized_scaffold_count\n"
                "CM1.1\t1.300000\t1300000\t1000000\t300000\t1\n"
                "CM2.1\t0.800000\t800000\t800000\t0\t0\n"
            )
            output = root / "unlocalized"

            plot_locations(
                location_file=location,
                lengths_file=lengths,
                assembly_mode="final",
                output_prefix=str(output),
                profile_name="coleoptera",
                label_threshold=1,
            )

            svg = output.with_suffix(".svg").read_text()
            self.assertIn("<!-- CM2.1 -->", svg)
            self.assertIn("<!-- Unlocalized scaffolds (appended) -->", svg)

    def test_every_duplicate_row_gets_an_exact_marker_line(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            location = root / "all_location.tsv"
            location.write_text(
                "buscoID\tquery_chr\tposition\tassigned_alg\tstatus\n"
                "a\tchr1\t100000\tC1\tassigned\n"
                "b\tchr1\t100100\tC2\tassigned\n"
                "b\tchr1\t100200\tC2\tassigned\n"
            )
            lengths = root / "assembly.fai"
            lengths.write_text("chr1\t1000000\n")

            with patch(
                "matplotlib.axes._axes.Axes.vlines", autospec=True
            ) as marker_lines:
                plot_locations(
                    location_file=location,
                    lengths_file=lengths,
                    assembly_mode="draft",
                    output_prefix=str(root / "duplicates"),
                    profile_name="coleoptera",
                    label_threshold=1,
                )

            self.assertEqual(marker_lines.call_count, 3)

    def test_window_labels_include_unassigned_denominator_and_reject_ties(self) -> None:
        profile = load_profile("diptera")
        locations = pd.DataFrame(
            {
                "query_chr": ["chr1"] * 12,
                "assigned_alg": ["d1"] * 5 + ["NA"] * 7,
                "position": list(range(1, 13)),
            }
        )
        labels = calculate_windowed_alg_labels(
            locations,
            profile=profile,
            window_mb=1,
            min_buscos=5,
            min_fraction=0.5,
        )
        self.assertEqual(labels, {})

        tied = pd.DataFrame(
            {
                "query_chr": ["chr1"] * 10,
                "assigned_alg": ["d1"] * 5 + ["d2"] * 5,
                "position": list(range(1, 11)),
            }
        )
        labels = calculate_windowed_alg_labels(
            tied,
            profile=profile,
            window_mb=1,
            min_buscos=5,
            min_fraction=0.5,
        )
        self.assertEqual(labels, {})


if __name__ == "__main__":
    unittest.main()
