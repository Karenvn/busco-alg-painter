from __future__ import annotations

import hashlib
import unittest

from busco_alg_painter.painter import build_ref_map
from busco_alg_painter.profiles import (
    infer_profile_from_labels,
    load_bundled_profiles,
    load_profile,
    profile_for_dataset,
    profile_from_taxonomy_ids,
)


class ProfileTests(unittest.TestCase):
    def test_all_bundled_profiles_are_valid(self) -> None:
        profiles = load_bundled_profiles()
        self.assertEqual(
            [profile.id for profile in profiles],
            ["merian", "diptera", "brachycera", "coleoptera"],
        )
        for profile in profiles:
            self.assertTrue(profile.reference_table.is_file())
            self.assertEqual(
                set(profile.alg_order),
                set(profile.palette()),
                f"incomplete default palette for {profile.id}",
            )

    def test_reference_assignment_counts(self) -> None:
        expected = {
            "merian": 4112,
            "diptera": 3804,
            "brachycera": 4785,
            "coleoptera": 3603,
        }
        for name, count in expected.items():
            profile = load_profile(name)
            assignments = build_ref_map(profile.reference_table, profile)
            self.assertEqual(len(assignments), count)

    def test_reference_checksums(self) -> None:
        reference_dir = load_profile("coleoptera").reference_table.parent
        checksum_file = reference_dir / "REFERENCE_CHECKSUMS.sha256"
        for line in checksum_file.read_text().splitlines():
            expected, filename = line.split(maxsplit=1)
            digest = hashlib.sha256((reference_dir / filename).read_bytes()).hexdigest()
            self.assertEqual(digest, expected, filename)

    def test_coleoptera_csv_bom_case_and_unassigned_rows(self) -> None:
        profile = load_profile("coleoptera")
        self.assertTrue(
            profile.reference_table.read_bytes().startswith(b"\xef\xbb\xbf")
        )
        assignments = build_ref_map(profile.reference_table, profile)
        self.assertEqual(assignments["100357at7041"], "C1")
        self.assertEqual(assignments["10053at7041"], "CX")
        self.assertNotIn("101435at7041", assignments)

    def test_auto_dataset_and_taxonomy_resolution(self) -> None:
        default = profile_for_dataset("diptera_odb12")
        self.assertIsNotNone(default)
        self.assertEqual(default.id, "diptera")

        specific = profile_from_taxonomy_ids("diptera_odb12", {7147, 7203, 999999})
        self.assertIsNotNone(specific)
        self.assertEqual(specific.id, "brachycera")

    def test_infer_profile_from_existing_labels(self) -> None:
        self.assertEqual(infer_profile_from_labels(["M1", "MZ"]).id, "merian")
        self.assertEqual(infer_profile_from_labels(["db1a", "db6"]).id, "brachycera")
        self.assertEqual(infer_profile_from_labels(["C1", "CX"]).id, "coleoptera")


if __name__ == "__main__":
    unittest.main()
