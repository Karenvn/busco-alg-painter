# busco-alg-painter

`busco-alg-painter` maps BUSCO genes to taxon-specific ancestral linkage
groups (ALGs) and plots those assignments along chromosomes or scaffolds. One
shared mapper and plotter supports Lepidoptera, Diptera and Coleoptera through
versioned profile files.

## Included profiles

| Profile | Groups | Required BUSCO dataset | Reference assignments |
| --- | ---: | --- | ---: |
| `merian` | MZ, M1–M31 | `lepidoptera_odb10` | 4,112 |
| `diptera` | d1–d6 | `diptera_odb12` | 3,804 |
| `brachycera` | db1a–db6 | `diptera_odb12` | 4,785 |
| `coleoptera` | C1–C7, CX | `coleoptera_odb12` | 3,603 |

Profiles define the expected BUSCO dataset, NCBI taxon, reference-table
format, ordered ALG labels, palette, legend and taxon-specific plot defaults.
They are stored under
`src/busco_alg_painter/data/profiles/`.

## Scientific references

Please cite the paper associated with the profile used:

- **Lepidoptera / Merian elements:** Wright, C. J., Stevens, L.,
  Mackintosh, A., Lawniczak, M. & Blaxter, M. (2024).
  [Comparative genomics reveals the dynamics of chromosome evolution in
  Lepidoptera](https://doi.org/10.1038/s41559-024-02329-4).
  *Nature Ecology & Evolution*.
- **Diptera and Brachycera ALGs:** Gries, J. *et al.* (2026).
  [340 dipteran genomes reveal the origin of Muller elements and sex
  chromosomes in Diptera](https://doi.org/10.64898/2026.06.01.729285).
  bioRxiv preprint.
- **Coleoptera ALGs:** Maulana, A. *et al.* (2026).
  [338 coleopteran genomes reveal exceptional rearrangement variation compared
  to other insect orders](https://doi.org/10.64898/2026.07.17.739156).
  bioRxiv preprint.

## Reference-table provenance

Fixed copies are bundled so that a run remains reproducible and does not
depend on a live download:

- Merian elements:
  [`charlottewright/lep_busco_painter`](https://github.com/charlottewright/lep_busco_painter)
- Diptera:
  [`ALGs_syngraph_diptera.tsv`](https://raw.githubusercontent.com/Obscuromics/diptera-ALGs/refs/heads/main/tables/ALGs_syngraph_diptera.tsv)
- Brachycera:
  [`ALGs_syngraph_brachycera.tsv`](https://raw.githubusercontent.com/Obscuromics/diptera-ALGs/refs/heads/main/tables/ALGs_syngraph_brachycera.tsv)
- Coleoptera:
  [`TableS2.csv`](https://raw.githubusercontent.com/Obscuromics/coleoptera-ALGs/refs/heads/main/tables/TableS2.csv)

SHA-256 values for the bundled snapshots are recorded in
`REFERENCE_CHECKSUMS.sha256`.

The Coleoptera source table contains 3,729
`coleoptera_odb12` BUSCOs: 3,603 assigned to C1–C7 or CX and 126 reported as
`unassigned`. The latter are deliberately excluded from the reference mapping.

## Installation

```bash
cd /Users/kh18/code/busco-alg-painter
python3 -m pip install -e .
```

Installed commands:

```text
busco-alg-painter
bap
diptera-busco-painter
dbp
merian-busco-painter
mbp
buscopainter
plot-buscopainter
```

The Diptera and Merian command names are compatibility entry points. New
workflows should use `busco-alg-painter` or `bap`.

## Quick start

Run the complete mapping and plotting workflow:

```bash
busco-alg-painter run \
  --query-table full_table.tsv \
  --profile coleoptera \
  --prefix output/
```

The same command works with `merian`, `diptera` or `brachycera`.

`--profile auto` reads the BUSCO dataset header:

- `lepidoptera_odb10` selects `merian`
- `coleoptera_odb12` selects `coleoptera`
- `diptera_odb12` selects `diptera`, or `brachycera` when an accession, taxid
  or taxonomy string establishes Brachycera membership

For example:

```bash
busco-alg-painter run \
  --query-table full_table.tsv \
  --profile auto \
  --accession GCA_965649395.1 \
  --prefix output/ \
  --write-summary
```

Without an accession, provide a local FASTA index for true scaffold lengths:

```bash
busco-alg-painter run \
  --query-table full_table.tsv \
  --profile merian \
  --lengths assembly.fa.fai \
  --assembly-mode draft \
  --prefix output/
```

When `--accession` is supplied, chromosome lengths and public GenBank
accessions are obtained from the NCBI Datasets API. Unlocalized scaffolds are
added to the length of their parent chromosome.

## BUSCO compatibility checks

BUSCO identifiers are dataset-specific. A profile/input mismatch is therefore
an error rather than a warning:

```text
BUSCO table reports coleoptera_odb12, but profile 'merian'
expects lepidoptera_odb10
```

`--allow-lineage-mismatch` is available for a reference table known to be
compatible with a differently named BUSCO dataset.

## Outputs

`run` produces:

- `all_location.tsv`
- `chrom_lengths.tsv` when `--accession` is used
- `summary.tsv` when `--write-summary` is used
- PNG and SVG plots

The standardized location columns are:

```text
buscoID    query_chr    position    assigned_alg    status
```

The plotter also accepts older Merian tables containing `assigned_chr`.

The lower-level commands are:

```bash
busco-alg-painter paint --help
busco-alg-painter plot --help
```

## Plot controls

Chromosome labels can either include every ALG represented by at least
`--label-threshold` BUSCOs or use the dominant assignment in genomic windows:

```bash
busco-alg-painter run \
  --query-table full_table.tsv \
  --profile auto \
  --prefix output/ \
  --label-window-mb 10 \
  --label-window-min-buscos 5 \
  --label-window-min-fraction 0.5
```

Merian additionally retains the named palettes `merianbow4`, `merianbow`,
`categorical` and `spectrum`:

```bash
busco-alg-painter plot \
  --file output/all_location.tsv \
  --profile merian \
  --palette merianbow4 \
  --prefix output/merian
```

## Batch workflow

`busco_to_algs.sh` processes a ToLID/accession table. It expects BUSCO results
at `${BUSCO_DIR}/${ToLID}/full_table.tsv`.

```bash
DATA_ROOT=/path/to/project_data \
BUSCO_DIR=/path/to/project_data/busco \
OUTPUT_DIR=/path/to/project_data/algs \
ACCESSION_FILE=/path/to/tolid_accessions.tsv \
PROFILE=auto \
bash busco_to_algs.sh
```

The accession table must contain `ToLID` and assembly accession as its first
two tab-separated columns.

## Custom profiles

Copy one of the bundled TOML files and edit its biological metadata, reference
schema, labels and palette. Relative reference-table paths are resolved from
the custom configuration file:

```bash
busco-alg-painter run \
  --config /path/to/my_profile.toml \
  --query-table full_table.tsv \
  --prefix output/
```

`--reference-table` overrides only the table. `--config` is needed when the
labels, palette or table schema also differ.

## Code provenance

This repository consolidates `merian-busco-painter` and
`diptera-busco-painter`. The Merian implementation was adapted from
[`charlottewright/lep_busco_painter`](https://github.com/charlottewright/lep_busco_painter);
its MIT licence and attribution are retained in `LICENSE`.

## Development

Run the test suite without installing:

```bash
PYTHONPATH=src MPLCONFIGDIR=/tmp/busco-alg-painter-mpl \
  python3 -m unittest discover -s tests -v
```
