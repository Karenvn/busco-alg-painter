"""Load and validate bundled or user-supplied ALG profiles."""

from __future__ import annotations

from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Any, Iterable

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 only
    import tomli as tomllib


BUNDLED_PROFILE_NAMES = ("merian", "diptera", "brachycera", "coleoptera")


@dataclass(frozen=True)
class PlotDefaults:
    bar_height: float = 0.45
    row_height_cm: float = 0.85
    min_plot_height_cm: float = 12.0
    tile_width_bp: int = 50_000


@dataclass(frozen=True)
class Profile:
    id: str
    title: str
    description: str
    paper_doi: str
    busco_dataset: str
    taxid: int
    taxonomy_names: tuple[str, ...]
    auto_priority: int
    auto_default: bool
    reference_table: Path
    reference_source: str
    delimiter: str
    busco_column: int
    alg_column: int
    alg_order: tuple[str, ...]
    legend_title: str
    legend_columns: int
    default_palette: str
    palettes: dict[str, dict[str, str]]
    plot: PlotDefaults
    config_path: Path

    @property
    def valid_labels(self) -> set[str]:
        return set(self.alg_order)

    @property
    def canonical_labels(self) -> dict[str, str]:
        return {label.casefold(): label for label in self.alg_order}

    def normalize_label(self, value: str) -> str | None:
        return self.canonical_labels.get(value.strip().casefold())

    def palette(self, name: str | None = None) -> dict[str, str]:
        palette_name = name or self.default_palette
        try:
            return self.palettes[palette_name]
        except KeyError as exc:
            choices = ", ".join(sorted(self.palettes))
            raise ValueError(
                f"Unknown palette {palette_name!r} for profile {self.id!r}; "
                f"choose from: {choices}"
            ) from exc


def _read_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as fh:
        return tomllib.load(fh)


def bundled_profile_path(name: str) -> Path:
    if name not in BUNDLED_PROFILE_NAMES:
        choices = ", ".join(BUNDLED_PROFILE_NAMES)
        raise ValueError(f"Unknown profile {name!r}; choose from: {choices}")
    resource = files("busco_alg_painter").joinpath(f"data/profiles/{name}.toml")
    return Path(str(resource))


def load_profile(
    name: str | None = None, *, config_path: str | Path | None = None
) -> Profile:
    """Load one bundled profile or a custom TOML profile."""
    if config_path is not None:
        path = Path(config_path).expanduser().resolve()
    else:
        if not name or name == "auto":
            raise ValueError("A concrete profile name or --config is required")
        path = bundled_profile_path(name)

    data = _read_toml(path)
    reference = Path(data["reference"]["table"])
    if not reference.is_absolute():
        reference = (path.parent / reference).resolve()

    alg_order = tuple(str(label) for label in data["labels"]["order"])
    if not alg_order or len(set(alg_order)) != len(alg_order):
        raise ValueError(f"{path}: labels.order must contain unique labels")

    palettes = {
        palette_name: {str(label): str(color) for label, color in colors.items()}
        for palette_name, colors in data["palettes"].items()
    }
    default_palette = str(data["profile"]["default_palette"])
    if default_palette not in palettes:
        raise ValueError(f"{path}: default palette {default_palette!r} is not defined")
    for palette_name, colors in palettes.items():
        missing = set(alg_order) - set(colors)
        if missing:
            raise ValueError(
                f"{path}: palette {palette_name!r} is missing: "
                + ", ".join(sorted(missing))
            )

    plot_data = data.get("plot", {})
    profile = Profile(
        id=str(data["profile"]["id"]),
        title=str(data["profile"]["title"]),
        description=str(data["profile"].get("description", "")),
        paper_doi=str(data["profile"].get("paper_doi", "")),
        busco_dataset=str(data["profile"]["busco_dataset"]),
        taxid=int(data["profile"]["taxid"]),
        taxonomy_names=tuple(
            str(value).casefold() for value in data["profile"].get("taxonomy_names", [])
        ),
        auto_priority=int(data["profile"].get("auto_priority", 0)),
        auto_default=bool(data["profile"].get("auto_default", False)),
        reference_table=reference,
        reference_source=str(data["reference"].get("source", "")),
        delimiter=str(data["reference"]["delimiter"]),
        busco_column=int(data["reference"]["busco_column"]),
        alg_column=int(data["reference"]["alg_column"]),
        alg_order=alg_order,
        legend_title=str(data["labels"]["legend_title"]),
        legend_columns=int(data["labels"].get("legend_columns", 1)),
        default_palette=default_palette,
        palettes=palettes,
        plot=PlotDefaults(
            bar_height=float(plot_data.get("bar_height", 0.45)),
            row_height_cm=float(plot_data.get("row_height_cm", 0.85)),
            min_plot_height_cm=float(plot_data.get("min_plot_height_cm", 12)),
            tile_width_bp=int(plot_data.get("tile_width_bp", 50_000)),
        ),
        config_path=path,
    )
    if name and name != "auto" and config_path is None and profile.id != name:
        raise ValueError(f"{path}: profile id must be {name!r}")
    return profile


def load_bundled_profiles() -> tuple[Profile, ...]:
    return tuple(load_profile(name) for name in BUNDLED_PROFILE_NAMES)


def profile_for_dataset(
    dataset: str, profiles: Iterable[Profile] | None = None
) -> Profile | None:
    candidates = [
        profile
        for profile in (profiles or load_bundled_profiles())
        if profile.busco_dataset == dataset and profile.auto_default
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda profile: profile.auto_priority)


def profile_from_taxonomy_ids(
    dataset: str,
    lineage_ids: set[int],
    profiles: Iterable[Profile] | None = None,
) -> Profile | None:
    candidates = [
        profile
        for profile in (profiles or load_bundled_profiles())
        if profile.busco_dataset == dataset and profile.taxid in lineage_ids
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda profile: profile.auto_priority)


def profile_from_taxonomy_text(
    dataset: str,
    taxonomy_text: str,
    profiles: Iterable[Profile] | None = None,
) -> Profile | None:
    text = taxonomy_text.casefold()
    candidates = [
        profile
        for profile in (profiles or load_bundled_profiles())
        if profile.busco_dataset == dataset
        and any(name in text for name in profile.taxonomy_names)
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda profile: profile.auto_priority)


def infer_profile_from_labels(labels: Iterable[str]) -> Profile:
    """Infer a bundled profile from ALG labels in an existing location table."""
    observed = {str(label).strip().casefold() for label in labels}
    scored = []
    for profile in load_bundled_profiles():
        valid = {label.casefold() for label in profile.alg_order}
        scored.append(
            (len(observed.intersection(valid)), profile.auto_priority, profile)
        )
    best_score, _, best_profile = max(scored, key=lambda item: (item[0], item[1]))
    if best_score == 0:
        raise ValueError(
            "Could not infer an ALG profile from assigned labels; pass --profile"
        )
    return best_profile
