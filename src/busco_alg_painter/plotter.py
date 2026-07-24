"""Plot BUSCO locations colored by configured ALG assignments."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager

from busco_alg_painter.profiles import (
    Profile,
    infer_profile_from_labels,
    load_profile,
)

LABEL_OFFSET_FACTOR = 0.02
LABEL_PADDING_FACTOR = 0.06
DEFAULT_PANEL_SIZE = 20
PLOT_BODY_WIDTH_CM = 22
MAX_PANEL_COLUMNS = 3
COMPACT_LABEL_THRESHOLD = 80
ASSEMBLY_MODES = ("auto", "final", "draft")
DEFAULT_LABEL_WRAP = 4
DEFAULT_LABEL_WINDOW_MIN_FRACTION = 0.5


def resolve_open_sans_font(env_var: str = "GENOMENOTES_FONT") -> str | None:
    """Locate a regular Open Sans font file."""
    import os

    explicit = os.environ.get(env_var)
    if explicit and Path(explicit).is_file():
        return explicit

    def pick_upright(paths: list[Path]) -> str | None:
        regular = [path for path in paths if "Regular" in str(path)]
        if regular:
            return str(sorted(regular)[0])
        upright = [path for path in paths if "italic" not in path.name.lower()]
        if upright:
            return str(sorted(upright)[0])
        return str(sorted(paths)[0]) if paths else None

    user_fonts = Path.home() / "Library" / "Fonts"
    chosen = pick_upright(list(user_fonts.glob("OpenSans*.ttf")))
    if chosen:
        return chosen

    package_root = Path(__file__).resolve().parent.parent
    font_dir = package_root / "assets" / "fonts"
    if font_dir.is_dir():
        hits: list[Path] = []
        for pattern in ("OpenSans-Regular.ttf", "OpenSans*.ttf", "open-sans*.ttf"):
            hits.extend(font_dir.glob(pattern))
        chosen = pick_upright(hits)
        if chosen:
            return chosen

    return None


def setup_font() -> None:
    """Configure Open Sans when available, otherwise use matplotlib defaults."""
    try:
        font_path = resolve_open_sans_font()
        if font_path:
            font_manager.fontManager.addfont(font_path)
            plt.rcParams["font.family"] = "Open Sans"
            plt.rcParams["font.style"] = "normal"
            plt.rcParams["font.weight"] = "normal"
            print(f"[INFO] Using Open Sans font: {font_path}")
        else:
            print("[WARN] Open Sans not found, using default font")
    except Exception as exc:
        print(f"[WARN] Could not load Open Sans: {exc}")


def detect_lengths_format(lengths_file: Path) -> str:
    """Return 'final' for Chrom/Length_Mb tables, or 'draft' for .fai files."""
    with open(lengths_file) as fh:
        first_line = fh.readline().strip().split("\t")

    if first_line[:2] == ["Chrom", "Length_Mb"]:
        return "final"
    if len(first_line) >= 2:
        return "draft"
    raise ValueError(f"Could not detect lengths file format: {lengths_file}")


def load_lengths(lengths_file: Path, assembly_mode: str = "auto") -> pd.DataFrame:
    """Load either final chrom_lengths.tsv or draft .fai lengths into bp."""
    if assembly_mode not in ASSEMBLY_MODES:
        raise ValueError(f"Unknown assembly mode: {assembly_mode}")

    detected = detect_lengths_format(lengths_file)
    if assembly_mode != "auto" and assembly_mode != detected:
        expected = "Chrom/Length_Mb TSV" if assembly_mode == "final" else ".fai"
        raise ValueError(
            f"{lengths_file} looks like {detected!r} lengths, but "
            f"--assembly-mode {assembly_mode} expects {expected}"
        )

    if detected == "final":
        chrom_lengths = pd.read_csv(lengths_file, sep="\t")
        chrom_lengths["length"] = chrom_lengths["Length_Mb"] * 1e6
        return chrom_lengths.rename(columns={"Chrom": "query_chr"})[
            ["query_chr", "length"]
        ]

    return pd.read_csv(
        lengths_file,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["query_chr", "length"],
    )


def normalize_location_columns(locations: pd.DataFrame) -> pd.DataFrame:
    """Support new assigned_alg output and older assigned_chr files."""
    if "assigned_alg" not in locations.columns and "assigned_chr" in locations.columns:
        locations = locations.rename(columns={"assigned_chr": "assigned_alg"})
    if "assigned_alg" not in locations.columns:
        raise ValueError("Location table must contain assigned_alg or assigned_chr")

    locations["query_chr"] = locations["query_chr"].str.replace(":.*", "", regex=True)
    locations["position"] = pd.to_numeric(locations["position"], errors="coerce")
    locations["assigned_alg"] = locations["assigned_alg"].astype(str)
    return locations


def canonicalize_alg_labels(locations: pd.DataFrame, profile: Profile) -> pd.DataFrame:
    """Use the exact label spelling and case declared by a profile."""
    locations = locations.copy()
    locations["assigned_alg"] = locations["assigned_alg"].map(
        lambda value: profile.normalize_label(value) or value
    )
    return locations


def load_data(
    location_file: Path, lengths_file: Path | None = None, assembly_mode: str = "auto"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load BUSCO locations and chromosome/scaffold lengths."""
    locations = pd.read_csv(location_file, sep="\t", keep_default_na=False)
    locations = normalize_location_columns(locations)

    if lengths_file:
        chrom_lengths = load_lengths(lengths_file, assembly_mode=assembly_mode)
    else:
        if assembly_mode in {"final", "draft"}:
            print(
                "[WARN] No lengths file supplied; estimating lengths from BUSCO "
                f"positions despite --assembly-mode {assembly_mode}"
            )
        chrom_lengths = locations.groupby("query_chr")["position"].max().reset_index()
        chrom_lengths["length"] = chrom_lengths["position"] * 1.05

    return locations, chrom_lengths


def format_alg_label(algs: list[str], wrap: int = DEFAULT_LABEL_WRAP) -> str:
    """Return an ALG label, optionally wrapped by label count."""
    if wrap <= 0:
        return "; ".join(algs)

    wrapped_lines = [
        "; ".join(algs[index : index + wrap]) for index in range(0, len(algs), wrap)
    ]
    return "\n".join(wrapped_lines)


def split_balanced(values: list[str], n_groups: int) -> list[list[str]]:
    """Split values into ordered groups whose sizes differ by at most one."""
    n_groups = max(1, n_groups)
    base_size, extra = divmod(len(values), n_groups)
    groups: list[list[str]] = []
    start = 0
    for group_index in range(n_groups):
        group_size = base_size + (1 if group_index < extra else 0)
        end = start + group_size
        groups.append(values[start:end])
        start = end
    return groups


def calculate_alg_labels(
    locations: pd.DataFrame,
    profile: Profile,
    threshold: int = 5,
    wrap: int = DEFAULT_LABEL_WRAP,
) -> dict[str, str]:
    """Return ALG labels ordered by median BUSCO position per chromosome."""
    counts = (
        locations.groupby(["query_chr", "assigned_alg"])
        .agg(n=("assigned_alg", "size"), median_position=("position", "median"))
        .reset_index()
    )
    counts = counts[counts["n"] >= threshold]

    order_index = {label: index for index, label in enumerate(profile.alg_order)}
    labels: dict[str, str] = {}
    for chrom in counts["query_chr"].unique():
        chrom_counts = counts[counts["query_chr"] == chrom].copy()
        chrom_counts["alg_sort"] = chrom_counts["assigned_alg"].map(
            lambda alg: order_index.get(alg, len(order_index))
        )
        chrom_counts = chrom_counts.sort_values(
            ["median_position", "alg_sort"], kind="stable"
        )
        labels[chrom] = format_alg_label(
            chrom_counts["assigned_alg"].tolist(), wrap=wrap
        )

    return labels


def collapse_consecutive(values: list[str]) -> list[str]:
    collapsed: list[str] = []
    for value in values:
        if not collapsed or collapsed[-1] != value:
            collapsed.append(value)
    return collapsed


def calculate_windowed_alg_labels(
    locations: pd.DataFrame,
    profile: Profile,
    window_mb: float,
    min_buscos: int = 5,
    min_fraction: float = DEFAULT_LABEL_WINDOW_MIN_FRACTION,
    wrap: int = DEFAULT_LABEL_WRAP,
) -> dict[str, str]:
    """Return labels from dominant ALGs in non-overlapping genomic windows."""
    if window_mb <= 0:
        raise ValueError("--label-window-mb must be greater than 0")
    if not 0 < min_fraction <= 1:
        raise ValueError("--label-window-min-fraction must be > 0 and <= 1")

    window_bp = window_mb * 1_000_000
    windowed = locations.dropna(subset=["position"]).copy()
    windowed["window"] = (windowed["position"] // window_bp).astype(int)

    counts = (
        windowed.groupby(["query_chr", "window", "assigned_alg"])
        .size()
        .reset_index(name="n")
    )
    if counts.empty:
        return {}

    order_index = {label: index for index, label in enumerate(profile.alg_order)}
    counts["alg_sort"] = counts["assigned_alg"].map(
        lambda alg: order_index.get(alg, len(order_index))
    )
    counts = counts.sort_values(
        ["query_chr", "window", "n", "alg_sort"],
        ascending=[True, True, False, True],
        kind="stable",
    )
    dominant = counts.drop_duplicates(["query_chr", "window"], keep="first").copy()
    totals = (
        windowed.groupby(["query_chr", "window"]).size().reset_index(name="window_n")
    )
    dominant = dominant.merge(totals, on=["query_chr", "window"])
    dominant["fraction"] = dominant["n"] / dominant["window_n"]
    dominant = dominant[
        (dominant["n"] >= min_buscos) & (dominant["fraction"] >= min_fraction)
    ]

    labels: dict[str, str] = {}
    for chrom in dominant["query_chr"].unique():
        chrom_windows = dominant[dominant["query_chr"] == chrom].sort_values("window")
        algs = collapse_consecutive(chrom_windows["assigned_alg"].tolist())
        if algs:
            labels[chrom] = format_alg_label(algs, wrap=wrap)
    return labels


def plot_alg_chromosomes(
    locations: pd.DataFrame,
    chrom_lengths: pd.DataFrame,
    output_prefix: str,
    profile: Profile,
    palette: str | None = None,
    label_threshold: int = 5,
    panel_size: int = DEFAULT_PANEL_SIZE,
    max_columns: int = MAX_PANEL_COLUMNS,
    label_wrap: int = DEFAULT_LABEL_WRAP,
    label_window_mb: float = 0,
    label_window_min_buscos: int = 5,
    label_window_min_fraction: float = DEFAULT_LABEL_WINDOW_MIN_FRACTION,
) -> None:
    """Create the main ALG plot with chromosome labels."""
    setup_font()
    Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)

    valid_labels = profile.valid_labels
    colors = profile.palette(palette)
    is_valid_alg = locations["assigned_alg"].isin(valid_labels)

    plotted_chroms = set(chrom_lengths["query_chr"].dropna().unique())
    chrom_lengths = chrom_lengths[
        chrom_lengths["query_chr"].isin(plotted_chroms)
    ].copy()
    if chrom_lengths.empty:
        raise ValueError("No plotted chromosomes/scaffolds have matching lengths")

    label_locations = locations[(locations["buscoID"] != "NA") & is_valid_alg].copy()
    if label_window_mb > 0:
        print(
            "[INFO] Labelling dominant ALGs in "
            f"{label_window_mb:g} Mb windows "
            f"(min BUSCOs: {label_window_min_buscos}, "
            f"min fraction: {label_window_min_fraction:g})"
        )
        alg_labels = calculate_windowed_alg_labels(
            label_locations,
            profile=profile,
            window_mb=label_window_mb,
            min_buscos=label_window_min_buscos,
            min_fraction=label_window_min_fraction,
            wrap=label_wrap,
        )
    else:
        alg_labels = calculate_alg_labels(
            label_locations,
            profile=profile,
            threshold=label_threshold,
            wrap=label_wrap,
        )
    chrom_order = chrom_lengths.sort_values("length", ascending=False)[
        "query_chr"
    ].tolist()

    n_chroms = len(chrom_order)
    print(
        f"[INFO] Plotting {len(locations)} BUSCOs across "
        f"{n_chroms} chromosomes/scaffolds..."
    )

    panel_size = max(1, int(panel_size))
    max_columns = max(1, int(max_columns))
    ncols = min(max_columns, max(1, math.ceil(n_chroms / panel_size)))
    panel_chroms_list = split_balanced(chrom_order, ncols)
    chroms_per_panel = max(
        (len(panel_chroms) for panel_chroms in panel_chroms_list), default=0
    )
    print(
        f"[INFO] Layout: {ncols} column(s), up to "
        f"{chroms_per_panel} chromosomes/scaffolds per column"
    )

    compact_layout = n_chroms > COMPACT_LABEL_THRESHOLD
    label_fontsize = 8 if compact_layout else 10
    panel_height = max(
        profile.plot.min_plot_height_cm,
        profile.plot.row_height_cm * min(chroms_per_panel, max(1, n_chroms)),
    )

    panel_limits: list[float] = []
    for panel_chroms in panel_chroms_list:
        if panel_chroms:
            panel_max_length = chrom_lengths[
                chrom_lengths["query_chr"].isin(panel_chroms)
            ]["length"].max()
            panel_limits.append(panel_max_length * (1 + LABEL_PADDING_FACTOR))
        else:
            panel_limits.append(1.0)

    fig_width = PLOT_BODY_WIDTH_CM / 2.54
    panel_widths = [1] * ncols if compact_layout else panel_limits
    fig, axes = plt.subplots(
        nrows=1,
        ncols=ncols,
        figsize=(fig_width, panel_height / 2.54),
        squeeze=False,
        gridspec_kw={"width_ratios": panel_widths},
    )
    axes = axes[0]

    for col_idx in range(ncols):
        ax = axes[col_idx]
        panel_chroms = panel_chroms_list[col_idx]
        if not panel_chroms:
            ax.axis("off")
            continue

        y_positions = {chrom: i for i, chrom in enumerate(reversed(panel_chroms))}
        panel_limit = panel_limits[col_idx]

        for chrom in panel_chroms:
            y = y_positions[chrom]
            length = chrom_lengths[chrom_lengths["query_chr"] == chrom][
                "length"
            ].values[0]
            bar_bottom = y - profile.plot.bar_height / 2

            rect = patches.Rectangle(
                (0, bar_bottom),
                length,
                profile.plot.bar_height,
                facecolor="white",
                edgecolor="black",
                linewidth=0.5,
            )
            ax.add_patch(rect)

            chrom_buscos = locations[
                (locations["query_chr"] == chrom)
                & locations["assigned_alg"].isin(valid_labels)
            ]
            for _, busco in chrom_buscos.iterrows():
                if pd.notna(busco["position"]):
                    color = colors.get(busco["assigned_alg"], (0.85, 0.85, 0.85))
                    tile = patches.Rectangle(
                        (
                            busco["position"] - profile.plot.tile_width_bp / 2,
                            bar_bottom,
                        ),
                        profile.plot.tile_width_bp,
                        profile.plot.bar_height,
                        facecolor=color,
                        edgecolor="none",
                    )
                    ax.add_patch(tile)

            if chrom in alg_labels:
                ax.text(
                    length * (1 + LABEL_OFFSET_FACTOR),
                    y,
                    alg_labels[chrom],
                    va="center",
                    ha="left",
                    fontsize=label_fontsize,
                    color="#333333",
                    linespacing=0.85,
                )

        ax.set_xlim(0, panel_limit)
        ax.set_ylim(-0.6, len(panel_chroms) - 0.4)
        ax.set_xlabel("")
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x / 1e6:.0f}"))
        ax.set_yticks([y_positions[chrom] for chrom in panel_chroms])
        ax.set_yticklabels(panel_chroms, fontsize=label_fontsize)
        ax.set_ylabel("")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    fig.supxlabel("Position (Mb)", fontsize=11)
    plt.tight_layout()

    labels = list(profile.alg_order)
    if is_valid_alg.any():
        legend_elements = [
            patches.Patch(facecolor=colors[label], label=label) for label in labels
        ]
        fig.legend(
            handles=legend_elements,
            title=profile.legend_title,
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            fontsize=10,
            title_fontsize=11,
            ncol=profile.legend_columns,
            columnspacing=1.2,
        )
    else:
        print(f"[WARN] No BUSCOs were assigned to known {profile.id} labels")

    for ext in ("png", "svg"):
        output_file = f"{output_prefix}.{ext}"
        dpi = 300 if ext == "png" else None
        plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
        print(f"[INFO] Saved: {output_file}")

    plt.close()


def plot_locations(
    location_file: Path,
    output_prefix: str,
    lengths_file: Path | None = None,
    assembly_mode: str = "auto",
    profile_name: str = "auto",
    config_path: Path | None = None,
    profile: Profile | None = None,
    palette: str | None = None,
    label_threshold: int = 5,
    panel_size: int = DEFAULT_PANEL_SIZE,
    max_columns: int = MAX_PANEL_COLUMNS,
    label_wrap: int = DEFAULT_LABEL_WRAP,
    label_window_mb: float = 0,
    label_window_min_buscos: int = 5,
    label_window_min_fraction: float = DEFAULT_LABEL_WINDOW_MIN_FRACTION,
) -> None:
    print("[INFO] Loading data...")
    locations, chrom_lengths = load_data(
        location_file, lengths_file, assembly_mode=assembly_mode
    )
    if profile is None:
        if config_path is not None:
            profile = load_profile(config_path=config_path)
        elif profile_name == "auto":
            profile = infer_profile_from_labels(locations["assigned_alg"])
            print(f"[INFO] Inferred plot profile: {profile.id}")
        else:
            profile = load_profile(profile_name)
    locations = canonicalize_alg_labels(locations, profile)
    print(
        f"[INFO] Loaded {len(locations)} BUSCOs and "
        f"{len(chrom_lengths)} chromosome/scaffold lengths..."
    )
    plot_alg_chromosomes(
        locations,
        chrom_lengths,
        output_prefix,
        profile=profile,
        palette=palette,
        label_threshold=label_threshold,
        panel_size=panel_size,
        max_columns=max_columns,
        label_wrap=label_wrap,
        label_window_mb=label_window_mb,
        label_window_min_buscos=label_window_min_buscos,
        label_window_min_fraction=label_window_min_fraction,
    )
    print("[INFO] Done.")
