from __future__ import annotations

import unittest

from busco_alg_painter.cli import build_parser


class CliTests(unittest.TestCase):
    def test_bar_alg_labels_boolean_option(self) -> None:
        parser = build_parser()

        default_args = parser.parse_args(["plot", "--file", "locations.tsv"])
        shown_args = parser.parse_args(
            ["plot", "--file", "locations.tsv", "--bar-alg-labels"]
        )
        hidden_args = parser.parse_args(
            ["plot", "--file", "locations.tsv", "--no-bar-alg-labels"]
        )

        self.assertIsNone(default_args.bar_alg_labels)
        self.assertTrue(shown_args.bar_alg_labels)
        self.assertFalse(hidden_args.bar_alg_labels)

    def test_estimated_lengths_require_explicit_opt_in(self) -> None:
        parser = build_parser()

        default_args = parser.parse_args(["plot", "--file", "locations.tsv"])
        estimated_args = parser.parse_args(
            [
                "plot",
                "--file",
                "locations.tsv",
                "--allow-estimated-lengths",
            ]
        )

        self.assertFalse(default_args.allow_estimated_lengths)
        self.assertTrue(estimated_args.allow_estimated_lengths)


if __name__ == "__main__":
    unittest.main()
