from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import pandas as pd

from busco_alg_painter.painter import build_ref_map, paint_buscos
from busco_alg_painter.plotter import (
    assigned_sex_chromosome_labels,
    load_lengths,
    plot_locations,
)
from busco_alg_painter.profiles import load_profile


class SexChromosomeTests(unittest.TestCase):
    def test_only_explicit_assigned_sex_chromosome_names_are_used(self) -> None:
        names = [
            "X",
            "Y",
            "x1",
            " Y2 ",
            "1",
            "chrX",
            "SUPER_Y",
            "CX",
            "X-arm",
            "",
            None,
        ]
        lengths = pd.DataFrame(
            {
                "query_chr": [f"CM{i}.1" for i in range(len(names))],
                "assigned_molecule": names,
            }
        )
        expected = {"CM0.1": "X", "CM1.1": "Y", "CM2.1": "X1", "CM3.1": "Y2"}
        for profile_name in ("diptera", "brachycera", "coleoptera"):
            with self.subTest(profile=profile_name):
                self.assertEqual(
                    assigned_sex_chromosome_labels(lengths, load_profile(profile_name)),
                    expected,
                )
        self.assertEqual(
            assigned_sex_chromosome_labels(lengths, load_profile("merian")), {}
        )

    def test_old_lengths_and_fai_do_not_infer_sex_from_sequence_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for filename, contents in (
                ("old.tsv", "Chrom\tLength_Mb\nX\t1\nY\t0.5\n"),
                ("assembly.fai", "X\t1000000\nY\t500000\n"),
            ):
                with self.subTest(filename=filename):
                    path = root / filename
                    path.write_text(contents)
                    lengths = load_lengths(path)
                    self.assertEqual(lengths["assigned_molecule"].tolist(), ["", ""])
                    self.assertEqual(
                        assigned_sex_chromosome_labels(
                            lengths, load_profile("coleoptera")
                        ),
                        {},
                    )

    @staticmethod
    def sequence_report() -> list[dict]:
        return [
            {
                "role": "assembled-molecule",
                "assigned_molecule_location_type": "Chromosome",
                "chr_name": name,
                "sequence_name": sequence,
                "genbank_accession": accession,
                "refseq_accession": refseq,
                "length": length,
            }
            for name, sequence, accession, refseq, length in (
                ("1", "Y", "CMA.1", "NCA.1", 2_000_000),
                ("X", "SUPER_X", "CMX.1", "NCX.1", 1_000_000),
                ("Y", "SUPER_Y", "CMY.1", "NCY.1", 500_000),
            )
        ] + [
            {
                "role": "unlocalized-scaffold",
                "assigned_molecule_location_type": "Chromosome",
                "chr_name": "X",
                "sequence_name": "unloc_X",
                "genbank_accession": "JAAX.1",
                "length": 10_000,
            },
            {
                "role": "unplaced-scaffold",
                "chr_name": "X1",
                "sequence_name": "unplaced",
                "genbank_accession": "JAAU.1",
                "length": 10_000,
            },
        ]

    def test_ncbi_metadata_survives_mapping_and_renders_beside_correct_bars(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for profile_name in ("diptera", "brachycera", "coleoptera"):
                with self.subTest(profile=profile_name):
                    profile = load_profile(profile_name)
                    ref_map = build_ref_map(profile.reference_table, profile)
                    # A dark ancestral element on an autosome must not cause X labelling.
                    busco_id = next(
                        key
                        for key, value in ref_map.items()
                        if value == profile.alg_order[-1]
                    )
                    query = root / f"{profile_name}_full_table.tsv"
                    query.write_text(
                        f"# The lineage dataset is: {profile.busco_dataset}\n"
                        + "\n".join(
                            f"{busco_id}\tDuplicated\t{sequence}\t100\t200"
                            for sequence in (
                                "CMX.1",
                                "NCX.1",
                                "SUPER_X",
                                "unloc_X",
                                "Y",
                            )
                        )
                        + "\n"
                    )
                    with patch(
                        "busco_alg_painter.painter.fetch_sequence_report",
                        return_value=self.sequence_report(),
                    ):
                        outputs = paint_buscos(
                            query,
                            root / profile_name,
                            profile_name=profile_name,
                            accession="GCA_TEST.1",
                        )
                    raw_lengths = pd.read_csv(outputs.chrom_lengths, sep="\t")
                    self.assertEqual(
                        raw_lengths.set_index("Chrom")["Assigned_Molecule"].to_dict(),
                        {"CMA.1": "1", "CMX.1": "X", "CMY.1": "Y"},
                    )
                    locations = pd.read_csv(
                        outputs.all_locations, sep="\t", keep_default_na=False
                    )
                    self.assertEqual((locations["query_chr"] == "CMX.1").sum(), 4)
                    self.assertEqual(
                        locations.loc[
                            locations["query_chr"] == "CMY.1", "position_status"
                        ].tolist(),
                        ["no-busco"],
                    )
                    prefix = root / f"{profile_name}_plot"
                    with patch("busco_alg_painter.plotter.plt.close"):
                        plot_locations(
                            outputs.all_locations,
                            str(prefix),
                            outputs.chrom_lengths,
                            profile=profile,
                        )
                    try:
                        axis = plt.gcf().axes[0]
                        self.assertEqual(
                            [
                                (text.get_text(), text.get_position())
                                for text in axis.texts
                            ],
                            [("X", (1_010_000 * 1.02, 1)), ("Y", (500_000 * 1.02, 0))],
                        )
                    finally:
                        plt.close()
                    self.assertTrue(prefix.with_suffix(".png").is_file())
                    self.assertIn("<!-- Y -->", prefix.with_suffix(".svg").read_text())

    def test_sex_labels_are_independent_of_alg_controls(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            lengths = root / "chrom_lengths.tsv"
            lengths.write_text(
                "Chrom\tLength_Mb\tAssigned_Molecule\nCMX.1\t1\tX\nCMY.1\t0.5\tY\n"
            )
            for profile_name, alg in (
                ("diptera", "d1"),
                ("brachycera", "db1a"),
                ("coleoptera", "CX"),
            ):
                location = root / f"{profile_name}.tsv"
                location.write_text(
                    "buscoID\tquery_chr\tposition\tassigned_alg\n"
                    f"a\tCMX.1\t100000\t{alg}\n"
                )
                for show_algs in (None, False, True):
                    with self.subTest(profile=profile_name, show_algs=show_algs):
                        with (
                            patch("busco_alg_painter.plotter.plt.savefig"),
                            patch("busco_alg_painter.plotter.plt.close"),
                        ):
                            plot_locations(
                                location,
                                str(root / "plot"),
                                lengths,
                                profile_name=profile_name,
                                show_bar_alg_labels=show_algs,
                                label_window_mb=10,
                                label_window_min_buscos=1,
                            )
                        try:
                            self.assertEqual(
                                [text.get_text() for text in plt.gcf().axes[0].texts],
                                [f"X ({alg})" if show_algs else "X", "Y"],
                            )
                        finally:
                            plt.close()

    def test_merian_render_is_unchanged_by_assigned_molecule_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            location = root / "locations.tsv"
            location.write_text(
                "buscoID\tquery_chr\tposition\tassigned_alg\n"
                "a\tCM1.1\t100000\tM1\nb\tCM2.1\t100000\tMZ\n"
            )
            for show_algs in (None, False, True):
                images = []
                for metadata in (False, True):
                    lengths = root / "lengths.tsv"
                    lengths.write_text(
                        "Chrom\tLength_Mb"
                        + ("\tAssigned_Molecule" if metadata else "")
                        + "\n"
                        + "CM1.1\t1"
                        + ("\tX" if metadata else "")
                        + "\n"
                        + "CM2.1\t0.5"
                        + ("\tY" if metadata else "")
                        + "\n"
                    )
                    output = root / "merian"
                    plot_locations(
                        location,
                        str(output),
                        lengths,
                        profile_name="merian",
                        show_bar_alg_labels=show_algs,
                        label_threshold=1,
                    )
                    images.append(output.with_suffix(".png").read_bytes())
                with self.subTest(show_algs=show_algs):
                    self.assertEqual(images[0], images[1])


if __name__ == "__main__":
    unittest.main()
