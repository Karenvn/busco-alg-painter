"""Command-line interface for busco-alg-painter."""

from __future__ import annotations

import argparse
from pathlib import Path

from busco_alg_painter import __version__
from busco_alg_painter.painter import PROFILE_CHOICES, paint_buscos
from busco_alg_painter.plotter import (
    ASSEMBLY_MODES,
    DEFAULT_LABEL_WRAP,
    DEFAULT_LABEL_WINDOW_MIN_FRACTION,
    DEFAULT_PANEL_SIZE,
    MAX_PANEL_COLUMNS,
    plot_locations,
)


def add_profile_arguments(
    parser: argparse.ArgumentParser, *, default_profile: str
) -> None:
    parser.add_argument(
        "--profile",
        "--lineage",
        dest="profile",
        choices=PROFILE_CHOICES,
        default=default_profile,
        help=(
            "Bundled ALG profile. 'auto' uses the BUSCO header and taxonomy "
            f"(default: {default_profile})"
        ),
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Custom profile TOML file; overrides --profile",
    )


def configure_paint_parser(
    parser: argparse.ArgumentParser,
    *,
    set_func: bool = True,
    prefix_help: str = "Output directory-like prefix or filename stem",
    default_profile: str = "auto",
) -> None:
    add_profile_arguments(parser, default_profile=default_profile)
    parser.add_argument(
        "--reference-table",
        "--reference_table",
        "-r",
        dest="reference_table",
        type=Path,
        default=None,
        help="Custom BUSCO-to-ALG table; uses the selected profile's schema",
    )
    parser.add_argument(
        "--allow-lineage-mismatch",
        action="store_true",
        help="Allow a BUSCO dataset that differs from the selected profile",
    )
    parser.add_argument(
        "--taxid",
        type=int,
        default=None,
        help="NCBI taxid used by --profile auto",
    )
    parser.add_argument(
        "--taxon-lineage",
        "--taxon_lineage",
        dest="taxon_lineage",
        default=None,
        help="Taxonomy text used by --profile auto",
    )
    parser.add_argument(
        "--query-table",
        "--query_table",
        "-q",
        dest="query_table",
        type=Path,
        required=True,
        help="BUSCO full_table.tsv for the assembly to paint",
    )
    parser.add_argument("--prefix", "-p", default="buscopainter", help=prefix_help)
    parser.add_argument(
        "--accession",
        "-a",
        help="Assembly accession for NCBI chromosome lengths and auto profile",
    )
    parser.add_argument(
        "--write-summary",
        "--write_summary",
        dest="write_summary",
        action="store_true",
        help="Write per-chromosome BUSCO count summary",
    )
    if set_func:
        parser.set_defaults(func=run_paint)


def add_plot_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--palette",
        default=None,
        help="Named palette from the selected profile (default: profile setting)",
    )
    parser.add_argument(
        "--label-threshold",
        type=int,
        default=5,
        help="Minimum BUSCOs for an ALG label to appear (default: 5)",
    )
    parser.add_argument(
        "--label-wrap",
        type=int,
        default=DEFAULT_LABEL_WRAP,
        help=(
            "Wrap ALG labels after this many labels per line "
            f"(default: {DEFAULT_LABEL_WRAP}; use 0 to disable)"
        ),
    )
    parser.add_argument(
        "--label-window-mb",
        type=float,
        default=0,
        help=(
            "Use dominant ALGs in non-overlapping windows of this size in Mb; "
            "0 keeps chromosome-wide labels"
        ),
    )
    parser.add_argument(
        "--label-window-min-buscos",
        type=int,
        default=5,
        help="Minimum assigned BUSCOs required to label a window (default: 5)",
    )
    parser.add_argument(
        "--label-window-min-fraction",
        type=float,
        default=DEFAULT_LABEL_WINDOW_MIN_FRACTION,
        help=(
            "Minimum fraction of assigned BUSCOs in a window that must belong "
            f"to the dominant ALG (default: {DEFAULT_LABEL_WINDOW_MIN_FRACTION})"
        ),
    )
    parser.add_argument(
        "--panel-size",
        type=int,
        default=DEFAULT_PANEL_SIZE,
        help=(
            "Target chromosomes/scaffolds per panel before splitting into "
            f"columns (default: {DEFAULT_PANEL_SIZE})"
        ),
    )
    parser.add_argument(
        "--max-columns",
        type=int,
        default=MAX_PANEL_COLUMNS,
        help=f"Maximum number of plot columns (default: {MAX_PANEL_COLUMNS})",
    )


def add_lengths_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-l",
        "--lengths",
        type=Path,
        default=None,
        help="Final chrom_lengths.tsv or draft .fai lengths file",
    )
    parser.add_argument(
        "--assembly-mode",
        choices=ASSEMBLY_MODES,
        default="auto",
        help=(
            "auto detects the lengths format; final expects Chrom/Length_Mb TSV; "
            "draft expects .fai"
        ),
    )


def add_paint_parser(
    subparsers: argparse._SubParsersAction, *, default_profile: str
) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "paint", help="Map BUSCO full_table rows to ancestral linkage groups"
    )
    configure_paint_parser(parser, default_profile=default_profile)
    return parser


def configure_plot_parser(
    parser: argparse.ArgumentParser, *, default_profile: str = "auto"
) -> None:
    add_profile_arguments(parser, default_profile=default_profile)
    parser.add_argument(
        "-f",
        "--file",
        type=Path,
        required=True,
        help="Painter all_location.tsv output",
    )
    add_lengths_options(parser)
    parser.add_argument(
        "-p",
        "--prefix",
        default="buscopainter",
        help="Output prefix for PNG/SVG plot files",
    )
    add_plot_options(parser)
    parser.set_defaults(func=run_plot)


def add_plot_parser(
    subparsers: argparse._SubParsersAction, *, default_profile: str
) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "plot", help="Plot ALG assignments for final or draft assemblies"
    )
    configure_plot_parser(parser, default_profile=default_profile)
    return parser


def configure_run_parser(
    parser: argparse.ArgumentParser, *, default_profile: str
) -> None:
    configure_paint_parser(
        parser,
        set_func=False,
        prefix_help=(
            "Output prefix for final plots and intermediate TSV files. "
            "If this is a directory, the plot stem is the directory name."
        ),
        default_profile=default_profile,
    )
    add_lengths_options(parser)
    add_plot_options(parser)
    parser.set_defaults(func=run_all)


def add_run_parser(
    subparsers: argparse._SubParsersAction, *, default_profile: str
) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "run", help="Map BUSCOs to ALGs and generate chromosome plots"
    )
    configure_run_parser(parser, default_profile=default_profile)
    return parser


def build_parser(
    *, prog: str = "busco-alg-painter", default_profile: str = "auto"
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Paint BUSCO full_table.tsv files with taxon-specific ancestral "
            "linkage-group assignments."
        ),
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    add_run_parser(subparsers, default_profile=default_profile)
    add_paint_parser(subparsers, default_profile=default_profile)
    add_plot_parser(subparsers, default_profile=default_profile)
    return parser


def run_paint(args: argparse.Namespace) -> None:
    outputs = paint_buscos(
        query_table=args.query_table,
        prefix=args.prefix,
        reference_table=args.reference_table,
        profile_name=args.profile,
        config_path=args.config,
        taxid=args.taxid,
        taxon_lineage=args.taxon_lineage,
        accession=args.accession,
        write_summary=args.write_summary,
        allow_lineage_mismatch=args.allow_lineage_mismatch,
    )
    print(f"[INFO] Profile: {outputs.profile.id}")
    print("[INFO] Outputs written:")
    print(f"[INFO]   {outputs.all_locations}")
    if outputs.wrote_lengths:
        print(f"[INFO]   {outputs.chrom_lengths}")
    if outputs.wrote_summary:
        print(f"[INFO]   {outputs.summary}")


def run_plot(args: argparse.Namespace) -> None:
    plot_locations(
        location_file=args.file,
        lengths_file=args.lengths,
        assembly_mode=args.assembly_mode,
        output_prefix=args.prefix,
        profile_name=args.profile,
        config_path=args.config,
        palette=args.palette,
        label_threshold=args.label_threshold,
        label_wrap=args.label_wrap,
        label_window_mb=args.label_window_mb,
        label_window_min_buscos=args.label_window_min_buscos,
        label_window_min_fraction=args.label_window_min_fraction,
        panel_size=args.panel_size,
        max_columns=args.max_columns,
    )


def plot_prefix_from_run_prefix(prefix: str | Path) -> str:
    prefix_text = str(prefix)
    prefix_path = Path(prefix)
    if prefix_text.endswith(("/", "\\")) or prefix_path.is_dir():
        plot_name = prefix_path.name or "buscopainter"
        return str(prefix_path / plot_name)
    return prefix_text


def run_all(args: argparse.Namespace) -> None:
    outputs = paint_buscos(
        query_table=args.query_table,
        prefix=args.prefix,
        reference_table=args.reference_table,
        profile_name=args.profile,
        config_path=args.config,
        taxid=args.taxid,
        taxon_lineage=args.taxon_lineage,
        accession=args.accession,
        write_summary=args.write_summary,
        allow_lineage_mismatch=args.allow_lineage_mismatch,
    )

    lengths_file = args.lengths
    assembly_mode = args.assembly_mode
    if lengths_file is None and outputs.wrote_lengths:
        lengths_file = outputs.chrom_lengths
        if assembly_mode == "auto":
            assembly_mode = "final"

    plot_prefix = plot_prefix_from_run_prefix(args.prefix)
    plot_locations(
        location_file=outputs.all_locations,
        lengths_file=lengths_file,
        assembly_mode=assembly_mode,
        output_prefix=plot_prefix,
        profile=outputs.profile,
        palette=args.palette,
        label_threshold=args.label_threshold,
        label_wrap=args.label_wrap,
        label_window_mb=args.label_window_mb,
        label_window_min_buscos=args.label_window_min_buscos,
        label_window_min_fraction=args.label_window_min_fraction,
        panel_size=args.panel_size,
        max_columns=args.max_columns,
    )
    print(f"[INFO] Profile: {outputs.profile.id}")
    print("[INFO] Run outputs written:")
    print(f"[INFO]   {outputs.all_locations}")
    if outputs.wrote_lengths:
        print(f"[INFO]   {outputs.chrom_lengths}")
    if outputs.wrote_summary:
        print(f"[INFO]   {outputs.summary}")
    print(f"[INFO]   {plot_prefix}.png")
    print(f"[INFO]   {plot_prefix}.svg")


def _main(
    argv: list[str] | None,
    *,
    prog: str,
    default_profile: str,
) -> None:
    parser = build_parser(prog=prog, default_profile=default_profile)
    args = parser.parse_args(argv)
    args.func(args)


def main(argv: list[str] | None = None) -> None:
    _main(argv, prog="busco-alg-painter", default_profile="auto")


def diptera_main(argv: list[str] | None = None) -> None:
    """Compatibility entry point for diptera-busco-painter and dbp."""
    _main(argv, prog="diptera-busco-painter", default_profile="auto")


def merian_main(argv: list[str] | None = None) -> None:
    """Compatibility entry point for merian-busco-painter and mbp."""
    _main(argv, prog="merian-busco-painter", default_profile="merian")


def paint_main(argv: list[str] | None = None) -> None:
    """Compatibility entry point for the former standalone mapper."""
    parser = argparse.ArgumentParser(description="BUSCO to ALG mapper")
    configure_paint_parser(parser, default_profile="auto")
    args = parser.parse_args(argv)
    run_paint(args)


def plot_main(argv: list[str] | None = None) -> None:
    """Compatibility entry point for the former standalone plotter."""
    parser = argparse.ArgumentParser(description="Plot BUSCO ALG assignments")
    configure_plot_parser(parser, default_profile="auto")
    args = parser.parse_args(argv)
    run_plot(args)
