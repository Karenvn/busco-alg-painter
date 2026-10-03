"""Map BUSCO full_table rows to configured ancestral linkage groups."""

from __future__ import annotations

import csv
import os
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import requests

from busco_alg_painter.profiles import (
    Profile,
    load_bundled_profiles,
    load_profile,
    profile_for_dataset,
    profile_from_taxonomy_ids,
    profile_from_taxonomy_text,
)

API_KEY = os.getenv("NCBI_API_KEY")
PROFILE_CHOICES = ("auto", "merian", "diptera", "brachycera", "coleoptera")
BUSCO_SEQUENCE_COORDINATES_RE = re.compile(r"^(?P<sequence>.+):\d+-\d+$")


@dataclass(frozen=True)
class PaintOutputs:
    all_locations: Path
    chrom_lengths: Path
    sequence_layout: Path
    summary: Path
    wrote_lengths: bool
    wrote_sequence_layout: bool
    wrote_summary: bool
    profile: Profile
    reference_table: Path
    busco_dataset: str | None
    mapped_buscos: int
    total_buscos: int


@dataclass(frozen=True)
class ChromosomeSpan:
    """A chromosome sequence plus scaffolds assigned but not localized to it."""

    chrom: str
    localized_length_bp: int
    unlocalized_length_bp: int
    unlocalized_scaffold_count: int
    assigned_molecule: str = ""

    @property
    def length_bp(self) -> int:
        return self.localized_length_bp + self.unlocalized_length_bp


@dataclass(frozen=True)
class SequencePlacement:
    """How one assembly sequence is represented in the chromosome plot."""

    source_sequence: str
    sequence_name: str
    role: str
    assigned_chromosome: str
    target_chrom: str | None
    length_bp: int
    offset_bp: int
    position_status: str


def parse_busco_dataset(path: Path) -> str | None:
    """Return the BUSCO lineage dataset named in a full_table.tsv header."""
    marker = "# The lineage dataset is:"
    with path.open(encoding="utf-8-sig") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith(marker):
                return line.removeprefix(marker).strip().split()[0]
    return None


def parse_busco_table(path: Path) -> tuple[list[tuple[str, str, int, int]], list[str]]:
    """Return BUSCO ID, chromosome, start and stop for Complete/Duplicated rows."""
    table: list[tuple[str, str, int, int]] = []
    chromosomes: set[str] = set()
    keep_status = {"Complete", "Duplicated"}

    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#") or len(row) < 5:
                continue
            busco_id, status, chrom, start, stop = row[:5]
            if status not in keep_status:
                continue
            coordinate_match = BUSCO_SEQUENCE_COORDINATES_RE.fullmatch(chrom)
            if coordinate_match is not None:
                chrom = coordinate_match.group("sequence")
            try:
                start_coord, end_coord = int(start), int(stop)
            except ValueError:
                continue
            table.append((busco_id, chrom, start_coord, end_coord))
            chromosomes.add(chrom)
    return table, sorted(chromosomes)


def _reference_rows(path: Path, delimiter: str):
    with path.open(newline="", encoding="utf-8-sig") as fh:
        if delimiter == "whitespace":
            for line in fh:
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    yield stripped.split()
            return

        reader = csv.reader(fh, delimiter=delimiter)
        for row in reader:
            if row and not row[0].startswith("#"):
                yield row


def build_ref_map(ref_path: Path, profile: Profile) -> dict[str, str]:
    """Return a BUSCO-ID to canonical ALG-label mapping."""
    if not ref_path.is_file():
        raise FileNotFoundError(f"Reference table not found: {ref_path}")

    ref_map: dict[str, str] = {}
    largest_column = max(profile.busco_column, profile.alg_column)
    for row in _reference_rows(ref_path, profile.delimiter):
        if len(row) <= largest_column:
            continue
        busco_id = row[profile.busco_column].strip()
        label = profile.normalize_label(row[profile.alg_column])
        if not busco_id or label is None:
            continue
        previous = ref_map.get(busco_id)
        if previous is not None and previous != label:
            raise ValueError(
                f"{ref_path}: BUSCO {busco_id!r} maps to both "
                f"{previous!r} and {label!r}"
            )
        ref_map[busco_id] = label

    if not ref_map:
        raise ValueError(
            f"No {profile.id} ALG assignments were read from {ref_path}; "
            "check the profile column and delimiter settings"
        )
    return ref_map


def build_location_rows(
    ref_map: dict[str, str],
    query_table: list[tuple[str, str, int, int]],
    placement_metadata: list[SequencePlacement] | None = None,
) -> tuple[list[str], int]:
    if placement_metadata is not None and len(placement_metadata) != len(query_table):
        raise ValueError("BUSCO rows and sequence-placement metadata are misaligned")

    rows = [
        "buscoID\tquery_chr\tposition\tassigned_alg\tstatus\t"
        "source_sequence\tsequence_role\tposition_status"
    ]
    mapped = 0
    for index, (busco_id, query_chr, start, end) in enumerate(query_table):
        position = (start + end) / 2
        assigned = ref_map.get(busco_id, "NA")
        status = "assigned" if assigned != "NA" else "unassigned"
        if assigned != "NA":
            mapped += 1
        if placement_metadata is None:
            source_sequence = query_chr
            sequence_role = "input-sequence"
            position_status = "native"
        else:
            placement = placement_metadata[index]
            source_sequence = placement.source_sequence
            sequence_role = placement.role
            position_status = placement.position_status
        rows.append(
            f"{busco_id}\t{query_chr}\t{position}\t{assigned}\t{status}\t"
            f"{source_sequence}\t{sequence_role}\t{position_status}"
        )
    return rows, mapped


def ncbi_json(url: str) -> dict:
    headers = {"accept": "application/json", "User-Agent": "busco-alg-painter"}
    params = {}
    if API_KEY:
        params["api_key"] = API_KEY
    response = requests.get(url, headers=headers, params=params, timeout=60)
    response.raise_for_status()
    return response.json()


def fetch_assembly_taxid(accession: str) -> int | None:
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{accession}/dataset_report"
    )
    payload = ncbi_json(url)
    reports = payload.get("reports", [])
    if not reports:
        return None
    tax_id = reports[0].get("organism", {}).get("tax_id")
    return int(tax_id) if tax_id else None


def fetch_taxonomy_lineage_ids(taxid: int) -> set[int]:
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{taxid}"
    payload = ncbi_json(url)
    nodes = payload.get("taxonomy_nodes", [])
    if not nodes:
        return {taxid}
    taxonomy = nodes[0].get("taxonomy", {})
    lineage = {int(item) for item in taxonomy.get("lineage", [])}
    lineage.add(int(taxonomy.get("tax_id", taxid)))
    return lineage


def choose_profile(
    query_table: Path,
    profile_name: str = "auto",
    config_path: Path | None = None,
    taxid: int | None = None,
    taxon_lineage: str | None = None,
    accession: str | None = None,
    allow_lineage_mismatch: bool = False,
) -> tuple[Profile, str | None]:
    """Resolve a profile and validate it against the BUSCO input header."""
    busco_dataset = parse_busco_dataset(query_table)

    if config_path is not None:
        selected = load_profile(config_path=config_path)
        print(f"[INFO] Using custom profile: {selected.config_path}")
    elif profile_name != "auto":
        selected = load_profile(profile_name)
        print(f"[INFO] Using requested profile: {selected.id}")
    else:
        if busco_dataset is None:
            raise ValueError(
                "Cannot select --profile auto because the BUSCO dataset header "
                "is missing; pass --profile explicitly"
            )

        profiles = load_bundled_profiles()
        selected = None
        if taxon_lineage:
            selected = profile_from_taxonomy_text(
                busco_dataset, taxon_lineage, profiles
            )
            if selected:
                print(
                    f"[INFO] Auto profile selected {selected.id} " "from taxonomy text"
                )

        resolved_taxid = taxid
        if selected is None and resolved_taxid is None and accession:
            try:
                resolved_taxid = fetch_assembly_taxid(accession)
            except requests.RequestException as exc:
                print(f"[WARN] Could not fetch NCBI taxonomy for {accession}: {exc}")

        if selected is None and resolved_taxid is not None:
            try:
                selected = profile_from_taxonomy_ids(
                    busco_dataset,
                    fetch_taxonomy_lineage_ids(resolved_taxid),
                    profiles,
                )
                if selected:
                    print(
                        f"[INFO] Auto profile selected {selected.id} "
                        f"from taxid {resolved_taxid}"
                    )
            except requests.RequestException as exc:
                print(
                    f"[WARN] Could not fetch NCBI taxonomy for "
                    f"taxid {resolved_taxid}: {exc}"
                )

        if selected is None:
            selected = profile_for_dataset(busco_dataset, profiles)
            if selected:
                print(
                    f"[INFO] Auto profile selected {selected.id} "
                    f"from BUSCO dataset {busco_dataset}"
                )

        if selected is None:
            raise ValueError(
                f"No bundled ALG profile supports BUSCO dataset {busco_dataset!r}"
            )

    if (
        busco_dataset
        and busco_dataset != selected.busco_dataset
        and not allow_lineage_mismatch
    ):
        raise ValueError(
            f"BUSCO table reports {busco_dataset}, but profile {selected.id!r} "
            f"expects {selected.busco_dataset}. Use --allow-lineage-mismatch "
            "only when the reference table is known to be compatible."
        )
    if busco_dataset is None:
        print("[WARN] BUSCO dataset header not found; compatibility was not checked")
    elif busco_dataset != selected.busco_dataset:
        print(
            f"[WARN] Allowing BUSCO dataset mismatch: input is {busco_dataset}, "
            f"profile expects {selected.busco_dataset}"
        )
    return selected, busco_dataset


def fetch_sequence_report(accession: str) -> list[dict]:
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{accession}/sequence_reports"
    )
    print(f"[INFO] Fetching chromosome info from NCBI for {accession}...")
    payload = ncbi_json(url)
    return payload.get("sequence_report", {}).get("records") or payload.get(
        "reports", []
    )


def sequence_name_to_genbank_map(records: list[dict]) -> dict[str, str]:
    """Map NCBI sequence names and accessions to GenBank accessions."""
    chrom_map: dict[str, str] = {}
    for rec in records:
        genbank = rec.get("genbank_accession")
        if not genbank:
            continue
        for field in ("genbank_accession", "refseq_accession", "sequence_name"):
            value = rec.get(field)
            if value:
                chrom_map[value] = genbank
    return chrom_map


def _sequence_id(record: dict) -> str | None:
    for field in ("genbank_accession", "refseq_accession", "sequence_name"):
        value = record.get(field)
        if value:
            return str(value)
    return None


def _sequence_length(record: dict, *, required: bool) -> int:
    try:
        length = int(record.get("length", 0))
    except (TypeError, ValueError) as exc:
        if required:
            raise ValueError(
                f"Invalid NCBI sequence length in record: {record}"
            ) from exc
        return 0
    if required and length <= 0:
        raise ValueError(f"Missing NCBI sequence length in record: {record}")
    return length


def _record_sort_key(record: dict) -> tuple[int, str]:
    try:
        sort_order = int(record.get("sort_order"))
    except (TypeError, ValueError):
        sort_order = 2**63 - 1
    return sort_order, _sequence_id(record) or ""


def _add_placement_aliases(
    aliases: dict[str, SequencePlacement],
    record: dict,
    placement: SequencePlacement,
) -> None:
    for field in ("genbank_accession", "refseq_accession", "sequence_name"):
        value = record.get(field)
        if not value:
            continue
        key = str(value)
        previous = aliases.get(key)
        if previous is not None and previous != placement:
            raise ValueError(f"NCBI sequence alias {key!r} refers to multiple records")
        aliases[key] = placement


def build_sequence_layout(
    records: list[dict],
) -> tuple[list[ChromosomeSpan], dict[str, SequencePlacement], list[SequencePlacement]]:
    """Build an explicit plot layout from an NCBI sequence report.

    The assembled chromosome sequence retains its native coordinates.
    Unlocalized scaffolds are appended after it in deterministic report/accession
    order. This keeps them in the EBP-style assigned chromosome length while
    making their lack of a known chromosomal coordinate explicit and auditable.
    Unplaced and other non-chromosome sequences are recorded but not assigned a
    chromosome bar.
    """
    main_by_name: dict[str, dict] = {}
    for record in records:
        if (
            record.get("role") == "assembled-molecule"
            and record.get("assigned_molecule_location_type") == "Chromosome"
        ):
            chrom_name = str(record.get("chr_name", ""))
            accession = _sequence_id(record)
            if not chrom_name or not accession:
                raise ValueError(
                    "NCBI assembled chromosome record lacks chr_name or accession"
                )
            if chrom_name in main_by_name:
                raise ValueError(
                    f"NCBI sequence report contains multiple records for "
                    f"chromosome {chrom_name!r}"
                )
            _sequence_length(record, required=True)
            main_by_name[chrom_name] = record

    if not main_by_name:
        raise ValueError("NCBI sequence report contains no assembled chromosomes")

    unlocalized_by_name: dict[str, list[dict]] = {
        chrom_name: [] for chrom_name in main_by_name
    }
    for record in records:
        if record.get("role") != "unlocalized-scaffold":
            continue
        chrom_name = str(record.get("chr_name", ""))
        if chrom_name in unlocalized_by_name:
            _sequence_length(record, required=True)
            unlocalized_by_name[chrom_name].append(record)

    aliases: dict[str, SequencePlacement] = {}
    layout: list[SequencePlacement] = []
    spans: list[ChromosomeSpan] = []
    represented_ids: set[str] = set()

    for chrom_name, main_record in main_by_name.items():
        target_chrom = _sequence_id(main_record)
        assert target_chrom is not None
        localized_length = _sequence_length(main_record, required=True)
        main_placement = SequencePlacement(
            source_sequence=target_chrom,
            sequence_name=str(main_record.get("sequence_name", target_chrom)),
            role="assembled-molecule",
            assigned_chromosome=chrom_name,
            target_chrom=target_chrom,
            length_bp=localized_length,
            offset_bp=0,
            position_status="chromosome",
        )
        layout.append(main_placement)
        represented_ids.add(target_chrom)
        _add_placement_aliases(aliases, main_record, main_placement)

        offset = localized_length
        unlocalized_records = sorted(
            unlocalized_by_name[chrom_name], key=_record_sort_key
        )
        for record in unlocalized_records:
            source_sequence = _sequence_id(record)
            if source_sequence is None:
                raise ValueError(
                    "NCBI unlocalized scaffold record lacks a sequence accession/name"
                )
            length = _sequence_length(record, required=True)
            placement = SequencePlacement(
                source_sequence=source_sequence,
                sequence_name=str(record.get("sequence_name", source_sequence)),
                role="unlocalized-scaffold",
                assigned_chromosome=chrom_name,
                target_chrom=target_chrom,
                length_bp=length,
                offset_bp=offset,
                position_status="unlocalized-appended",
            )
            layout.append(placement)
            represented_ids.add(source_sequence)
            _add_placement_aliases(aliases, record, placement)
            offset += length

        spans.append(
            ChromosomeSpan(
                chrom=target_chrom,
                localized_length_bp=localized_length,
                unlocalized_length_bp=offset - localized_length,
                unlocalized_scaffold_count=len(unlocalized_records),
                assigned_molecule=chrom_name,
            )
        )

    for record in records:
        source_sequence = _sequence_id(record)
        if not source_sequence or source_sequence in represented_ids:
            continue
        role = str(record.get("role", "unknown"))
        chrom_name = str(record.get("chr_name", ""))
        if role == "unplaced-scaffold":
            position_status = "unplaced-excluded"
        elif role == "unlocalized-scaffold":
            position_status = "unlocalized-without-parent"
        else:
            position_status = "non-chromosome-excluded"
        placement = SequencePlacement(
            source_sequence=source_sequence,
            sequence_name=str(record.get("sequence_name", source_sequence)),
            role=role,
            assigned_chromosome=chrom_name,
            target_chrom=None,
            length_bp=_sequence_length(record, required=False),
            offset_bp=0,
            position_status=position_status,
        )
        layout.append(placement)
        represented_ids.add(source_sequence)
        _add_placement_aliases(aliases, record, placement)

    spans.sort(key=lambda span: -span.length_bp)
    return spans, aliases, layout


def remap_query_chromosomes(
    query_rows: list[tuple[str, str, int, int]], chrom_map: dict[str, str]
) -> tuple[list[tuple[str, str, int, int]], int]:
    """Return query rows with sequence IDs converted to GenBank accessions."""
    remapped: list[tuple[str, str, int, int]] = []
    changed = 0
    for busco_id, query_chr, start, end in query_rows:
        mapped_chr = chrom_map.get(query_chr, query_chr)
        if mapped_chr != query_chr:
            changed += 1
        remapped.append((busco_id, mapped_chr, start, end))
    return remapped, changed


def place_query_buscos(
    query_rows: list[tuple[str, str, int, int]],
    placements: dict[str, SequencePlacement],
) -> tuple[list[tuple[str, str, int, int]], list[SequencePlacement], int]:
    """Map BUSCO coordinates into the explicit chromosome plot layout."""
    remapped: list[tuple[str, str, int, int]] = []
    metadata: list[SequencePlacement] = []
    changed = 0
    for busco_id, query_chr, start, end in query_rows:
        placement = placements.get(query_chr)
        if placement is None:
            placement = SequencePlacement(
                source_sequence=query_chr,
                sequence_name=query_chr,
                role="unknown",
                assigned_chromosome="",
                target_chrom=None,
                length_bp=0,
                offset_bp=0,
                position_status="unmatched",
            )
        mapped_chr = placement.target_chrom or placement.source_sequence
        mapped_start = start + placement.offset_bp
        mapped_end = end + placement.offset_bp
        if mapped_chr != query_chr or placement.offset_bp:
            changed += 1
        remapped.append((busco_id, mapped_chr, mapped_start, mapped_end))
        metadata.append(placement)
    return remapped, metadata, changed


def chrom_lengths_with_unloc(records: list[dict]) -> list[tuple[str, int]]:
    """Return main GenBank chromosome lengths including unlocalized scaffolds."""
    print("[INFO] Using NCBI GenBank accessions for chromosome labels")
    spans, _, _ = build_sequence_layout(records)
    return [(span.chrom, span.length_bp) for span in spans]


def sequence_layout_lines(layout: list[SequencePlacement]) -> list[str]:
    lines = [
        "source_sequence\tsequence_name\tsequence_role\tassigned_chromosome\t"
        "plot_chrom\tlength_bp\toffset_bp\tplotted_start_bp\tplotted_end_bp\t"
        "position_status"
    ]
    for placement in layout:
        plotted_start = (
            str(placement.offset_bp + 1) if placement.target_chrom is not None else "NA"
        )
        plotted_end = (
            str(placement.offset_bp + placement.length_bp)
            if placement.target_chrom is not None
            else "NA"
        )
        lines.append(
            f"{placement.source_sequence}\t{placement.sequence_name}\t"
            f"{placement.role}\t{placement.assigned_chromosome or 'NA'}\t"
            f"{placement.target_chrom or 'NA'}\t{placement.length_bp}\t"
            f"{placement.offset_bp}\t{plotted_start}\t{plotted_end}\t"
            f"{placement.position_status}"
        )
    return lines


def write_tsv(lines: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def resolve_output_paths(prefix: str | Path) -> tuple[Path, Path, Path, Path]:
    """Resolve output paths from either a directory-like prefix or file stem."""
    prefix_text = str(prefix)
    prefix_path = Path(prefix)
    if prefix_text.endswith(("/", "\\")) or prefix_path.is_dir():
        out_dir = prefix_path
        return (
            out_dir / "all_location.tsv",
            out_dir / "chrom_lengths.tsv",
            out_dir / "sequence_layout.tsv",
            out_dir / "summary.tsv",
        )

    out_dir = prefix_path.parent
    stem = prefix_path.name
    return (
        out_dir / f"{stem}_all_location.tsv",
        out_dir / f"{stem}_chrom_lengths.tsv",
        out_dir / f"{stem}_sequence_layout.tsv",
        out_dir / f"{stem}_summary.tsv",
    )


def paint_buscos(
    query_table: Path,
    prefix: str | Path,
    reference_table: Path | None = None,
    profile_name: str = "auto",
    config_path: Path | None = None,
    taxid: int | None = None,
    taxon_lineage: str | None = None,
    accession: str | None = None,
    write_summary: bool = False,
    allow_lineage_mismatch: bool = False,
) -> PaintOutputs:
    """Run the mapper workflow and return generated output paths."""
    query_table = Path(query_table)
    out_all, out_len, out_layout, out_sum = resolve_output_paths(prefix)
    out_all.parent.mkdir(parents=True, exist_ok=True)

    profile, busco_dataset = choose_profile(
        query_table=query_table,
        profile_name=profile_name,
        config_path=config_path,
        taxid=taxid,
        taxon_lineage=taxon_lineage,
        accession=accession,
        allow_lineage_mismatch=allow_lineage_mismatch,
    )
    ref_path = Path(reference_table) if reference_table else profile.reference_table
    print(f"[INFO] Using {profile.id} ALG table: {ref_path}")

    ref_map = build_ref_map(ref_path, profile)
    query_rows, query_chromosomes = parse_busco_table(query_table)
    chrom_order = query_chromosomes.copy()
    placement_metadata: list[SequencePlacement] | None = None

    wrote_len = False
    wrote_layout = False
    if accession:
        sequence_report = fetch_sequence_report(accession)
        spans, placements, layout = build_sequence_layout(sequence_report)
        query_rows, placement_metadata, remapped_chroms = place_query_buscos(
            query_rows, placements
        )
        if remapped_chroms:
            print(
                f"[INFO] Remapped {remapped_chroms} BUSCO rows "
                "into the NCBI chromosome layout"
            )
        print("[INFO] Using NCBI GenBank accessions for chromosome labels")
        length_lines = [
            "Chrom\tLength_Mb\tLength_bp\tLocalized_Length_bp\t"
            "Unlocalized_Length_bp\tUnlocalized_scaffold_count\tAssigned_Molecule"
        ] + [
            f"{span.chrom}\t{span.length_bp / 1e6:.6f}\t{span.length_bp}\t"
            f"{span.localized_length_bp}\t{span.unlocalized_length_bp}\t"
            f"{span.unlocalized_scaffold_count}\t{span.assigned_molecule}"
            for span in spans
        ]
        write_tsv(length_lines, out_len)
        write_tsv(sequence_layout_lines(layout), out_layout)
        wrote_len = True
        wrote_layout = True
        chrom_order = [span.chrom for span in spans]

    all_rows, mapped_buscos = build_location_rows(
        ref_map, query_rows, placement_metadata=placement_metadata
    )
    if query_rows:
        pct_mapped = mapped_buscos / len(query_rows) * 100
        print(
            f"[INFO] Mapped {mapped_buscos}/{len(query_rows)} BUSCO rows "
            f"to {profile.legend_title} ({pct_mapped:.1f}%)"
        )
    else:
        print("[WARN] No Complete or Duplicated BUSCO rows found")

    query_chroms = {chrom for _, chrom, _, _ in query_rows}
    missing = [chrom for chrom in chrom_order if chrom not in query_chroms]
    for chrom in missing:
        all_rows.append(
            f"NA\t{chrom}\tNA\tNA\tunassigned\t{chrom}\t" "assembled-molecule\tno-busco"
        )
    write_tsv(all_rows, out_all)

    if placement_metadata is not None:
        unlocalized_count = sum(
            placement.position_status == "unlocalized-appended"
            for placement in placement_metadata
        )
        excluded_count = sum(
            placement.target_chrom is None for placement in placement_metadata
        )
        if unlocalized_count:
            print(
                f"[INFO] Appended and retained {unlocalized_count} BUSCO rows "
                "from unlocalized scaffolds"
            )
        if excluded_count:
            print(
                f"[WARN] {excluded_count} BUSCO rows are on unplaced or other "
                "non-chromosome sequences and are recorded but excluded from "
                "the chromosome plot"
            )

    wrote_sum = False
    if write_summary:
        counts = Counter(chrom for _, chrom, _, _ in query_rows)
        counts.update({chrom: 0 for chrom in missing})
        localized_counts: Counter[str] = Counter()
        unlocalized_counts: Counter[str] = Counter()
        if placement_metadata is None:
            localized_counts.update(chrom for _, chrom, _, _ in query_rows)
        else:
            for row, placement in zip(query_rows, placement_metadata):
                chrom = row[1]
                if placement.position_status == "unlocalized-appended":
                    unlocalized_counts[chrom] += 1
                elif placement.target_chrom is not None:
                    localized_counts[chrom] += 1
        summary_lines = [
            "query_chr\tbusco_hits\tlocalized_busco_hits\t" "unlocalized_busco_hits"
        ] + [
            f"{chrom}\t{counts[chrom]}\t{localized_counts[chrom]}\t"
            f"{unlocalized_counts[chrom]}"
            for chrom in chrom_order
        ]
        write_tsv(summary_lines, out_sum)
        wrote_sum = True

    return PaintOutputs(
        all_locations=out_all,
        chrom_lengths=out_len,
        sequence_layout=out_layout,
        summary=out_sum,
        wrote_lengths=wrote_len,
        wrote_sequence_layout=wrote_layout,
        wrote_summary=wrote_sum,
        profile=profile,
        reference_table=ref_path,
        busco_dataset=busco_dataset,
        mapped_buscos=mapped_buscos,
        total_buscos=len(query_rows),
    )
