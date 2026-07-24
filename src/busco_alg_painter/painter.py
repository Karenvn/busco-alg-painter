"""Map BUSCO full_table rows to configured ancestral linkage groups."""

from __future__ import annotations

import csv
import os
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


@dataclass(frozen=True)
class PaintOutputs:
    all_locations: Path
    chrom_lengths: Path
    summary: Path
    wrote_lengths: bool
    wrote_summary: bool
    profile: Profile
    reference_table: Path
    busco_dataset: str | None
    mapped_buscos: int
    total_buscos: int


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
    ref_map: dict[str, str], query_table: list[tuple[str, str, int, int]]
) -> tuple[list[str], int]:
    rows = ["buscoID\tquery_chr\tposition\tassigned_alg\tstatus"]
    mapped = 0
    for busco_id, query_chr, start, end in query_table:
        position = (start + end) / 2
        assigned = ref_map.get(busco_id, "NA")
        status = "assigned" if assigned != "NA" else "unassigned"
        if assigned != "NA":
            mapped += 1
        rows.append(f"{busco_id}\t{query_chr}\t{position}\t{assigned}\t{status}")
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


def chrom_lengths_with_unloc(records: list[dict]) -> list[tuple[str, int]]:
    """Return main GenBank chromosome lengths including unlocalized scaffolds."""
    print("[INFO] Using NCBI GenBank accessions for chromosome labels")
    main_acc: dict[str, str] = {}
    for rec in records:
        if (
            rec.get("role") == "assembled-molecule"
            and rec.get("assigned_molecule_location_type") == "Chromosome"
        ):
            genbank = rec.get("genbank_accession")
            if genbank:
                main_acc[rec["chr_name"]] = genbank

    bp_tot: dict[str, int] = {acc: 0 for acc in main_acc.values()}
    for rec in records:
        role = rec.get("role")
        loc = rec.get("assigned_molecule_location_type", "")
        if role == "assembled-molecule" and loc == "Chromosome":
            acc = rec.get("genbank_accession")
            if acc in bp_tot:
                bp_tot[acc] += int(rec.get("length", 0))
        elif role == "unlocalized-scaffold":
            acc = main_acc.get(rec.get("chr_name"))
            if acc:
                bp_tot[acc] += int(rec.get("length", 0))

    return sorted(bp_tot.items(), key=lambda item: -item[1])


def write_tsv(lines: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def resolve_output_paths(prefix: str | Path) -> tuple[Path, Path, Path]:
    """Resolve output paths from either a directory-like prefix or file stem."""
    prefix_text = str(prefix)
    prefix_path = Path(prefix)
    if prefix_text.endswith(("/", "\\")) or prefix_path.is_dir():
        out_dir = prefix_path
        return (
            out_dir / "all_location.tsv",
            out_dir / "chrom_lengths.tsv",
            out_dir / "summary.tsv",
        )

    out_dir = prefix_path.parent
    stem = prefix_path.name
    return (
        out_dir / f"{stem}_all_location.tsv",
        out_dir / f"{stem}_chrom_lengths.tsv",
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
    out_all, out_len, out_sum = resolve_output_paths(prefix)
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

    wrote_len = False
    if accession:
        sequence_report = fetch_sequence_report(accession)
        query_rows, remapped_chroms = remap_query_chromosomes(
            query_rows, sequence_name_to_genbank_map(sequence_report)
        )
        if remapped_chroms:
            print(
                f"[INFO] Remapped {remapped_chroms} BUSCO rows "
                "to NCBI GenBank accessions"
            )
        pairs = chrom_lengths_with_unloc(sequence_report)
        length_lines = ["Chrom\tLength_Mb"] + [
            f"{chrom}\t{basepairs / 1e6:.3f}" for chrom, basepairs in pairs
        ]
        write_tsv(length_lines, out_len)
        wrote_len = True
        chrom_order = [chrom for chrom, _ in pairs]

    all_rows, mapped_buscos = build_location_rows(ref_map, query_rows)
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
        all_rows.append(f"NA\t{chrom}\tNA\tNA\tunassigned")
    write_tsv(all_rows, out_all)

    wrote_sum = False
    if write_summary:
        counts = Counter(chrom for _, chrom, _, _ in query_rows)
        counts.update({chrom: 0 for chrom in missing})
        summary_lines = ["query_chr\tbusco_hits"] + [
            f"{chrom}\t{counts[chrom]}" for chrom in chrom_order
        ]
        write_tsv(summary_lines, out_sum)
        wrote_sum = True

    return PaintOutputs(
        all_locations=out_all,
        chrom_lengths=out_len,
        summary=out_sum,
        wrote_lengths=wrote_len,
        wrote_summary=wrote_sum,
        profile=profile,
        reference_table=ref_path,
        busco_dataset=busco_dataset,
        mapped_buscos=mapped_buscos,
        total_buscos=len(query_rows),
    )
