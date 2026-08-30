#!/usr/bin/env bash

# LSF worker for one fresh diptera_odb12 BUSCO run.

set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 TOLID FASTA DESTINATION LINEAGE CPU" >&2
  exit 2
fi

TOLID="$1"
FASTA="$2"
DEST="$3"
LINEAGE="$4"
CPU="$5"

[[ "$TOLID" =~ ^[A-Za-z0-9._-]+$ ]] || {
  echo "Invalid ToLID: $TOLID" >&2
  exit 2
}
[[ -s "$FASTA" ]] || {
  echo "Input FASTA not found or empty: $FASTA" >&2
  exit 2
}
[[ -d "$LINEAGE" ]] || {
  echo "BUSCO lineage not found: $LINEAGE" >&2
  exit 2
}
[[ "$CPU" =~ ^[1-9][0-9]*$ ]] || {
  echo "Invalid CPU count: $CPU" >&2
  exit 2
}

LINEAGE_NAME="$(basename "$LINEAGE")"
[[ "$LINEAGE_NAME" == "diptera_odb12" ]] || {
  echo "Expected diptera_odb12 lineage, found: $LINEAGE_NAME" >&2
  exit 2
}

mkdir -p -- "$DEST"

if ! type module >/dev/null 2>&1 && [[ -r /etc/profile ]]; then
  # LSF may start a non-login shell; initialise Environment Modules if needed.
  set +u
  # shellcheck disable=SC1091
  source /etc/profile
  set -u
fi
type module >/dev/null 2>&1 || {
  echo "The module command is not available" >&2
  exit 1
}
module load busco/5.8.2--pyhdfd78af_0
command -v busco >/dev/null 2>&1 || {
  echo "busco is not available after loading the module" >&2
  exit 1
}

OUT_NAME="BUSCO_${TOLID}"
if [[ "${LSB_JOBID:-}" =~ ^[0-9]+$ ]]; then
  JOB_TAG="$LSB_JOBID"
else
  JOB_TAG="manual_$$"
fi
BUSCO_TMP_BASE="${BUSCO_TMP_BASE:-/nfs/treeoflife-01/teams/grit/users/kh18/tmp}"
TMPDIR="${BUSCO_TMP_BASE%/}/busco_${JOB_TAG}_${TOLID}"
OUT_DIR="${TMPDIR}/${OUT_NAME}"
export TMPDIR

mkdir -p -- "$TMPDIR"

export OMP_NUM_THREADS="$CPU"
export OPENBLAS_NUM_THREADS="$CPU"
export MKL_NUM_THREADS="$CPU"
export NUMEXPR_NUM_THREADS="$CPU"

busco_table_lineage() {
  local table="$1"

  awk '
    /^# The lineage dataset is:/ {
      sub(/^# The lineage dataset is:[[:space:]]*/, "")
      print $1
      exit
    }
  ' "$table"
}

finish() {
  local rc=$?
  local run_dir="${OUT_DIR}/run_${LINEAGE_NAME}"
  local full_table=""
  local observed_lineage=""
  local published=0
  local submitted_job_id=""
  local summary

  trap - EXIT
  set +e

  if [[ -s "${run_dir}/full_table.tsv" ]]; then
    full_table="${run_dir}/full_table.tsv"
  else
    full_table="$(find "$OUT_DIR" -maxdepth 4 -type f -name full_table.tsv \
      -print 2>/dev/null | sort | head -n 1)"
  fi

  if [[ $rc -eq 0 ]]; then
    if [[ -n "$full_table" && -s "$full_table" ]]; then
      observed_lineage="$(busco_table_lineage "$full_table")"
      if [[ "$observed_lineage" == "$LINEAGE_NAME" ]]; then
        if cp -f -- "$full_table" "${DEST}/full_table.tsv.tmp" &&
           mv -f -- "${DEST}/full_table.tsv.tmp" "${DEST}/full_table.tsv"; then
          published=1
        else
          echo "Could not publish BUSCO table to ${DEST}" >&2
          rc=1
        fi
      else
        echo "BUSCO table lineage mismatch: expected ${LINEAGE_NAME}, found ${observed_lineage:-unknown}" >&2
        rc=1
      fi
    else
      echo "BUSCO completed without producing full_table.tsv" >&2
      rc=1
    fi
  fi

  for summary in "${OUT_DIR}"/short_summary.*; do
    [[ -f "$summary" ]] && cp -f -- "$summary" "$DEST/"
  done

  if [[ $rc -eq 0 && $published -ne 1 ]]; then
    rc=1
  fi

  rm -f -- "${DEST}/full_table.tsv.tmp"
  if [[ $rc -eq 0 ]]; then
    {
      printf 'completed_at\t%s\n' "$(date +'%Y-%m-%dT%H:%M:%S%z')"
      printf 'job_id\t%s\n' "$JOB_TAG"
      printf 'input\t%s\n' "$FASTA"
      printf 'lineage\t%s\n' "$LINEAGE_NAME"
    } > "${DEST}/BUSCO_COMPLETE"
    rm -f -- "${DEST}/BUSCO_FAILED"
  else
    {
      printf 'failed_at\t%s\n' "$(date +'%Y-%m-%dT%H:%M:%S%z')"
      printf 'exit_status\t%s\n' "$rc"
      printf 'job_id\t%s\n' "$JOB_TAG"
      printf 'input\t%s\n' "$FASTA"
      printf 'lineage\t%s\n' "$LINEAGE_NAME"
    } > "${DEST}/BUSCO_FAILED"
  fi

  if [[ -s "${DEST}/BUSCO_SUBMITTED" ]]; then
    submitted_job_id="$(awk -F '\t' '$1 == "job_id" {print $2; exit}' \
      "${DEST}/BUSCO_SUBMITTED")"
    if [[ "$submitted_job_id" == "$JOB_TAG" ]]; then
      rm -f -- "${DEST}/BUSCO_SUBMITTED"
    fi
  fi
  rm -rf -- "$TMPDIR"
  exit "$rc"
}
trap finish EXIT

{
  printf 'started_at\t%s\n' "$(date +'%Y-%m-%dT%H:%M:%S%z')"
  printf 'job_id\t%s\n' "$JOB_TAG"
  printf 'input\t%s\n' "$FASTA"
  printf 'lineage\t%s\n' "$LINEAGE_NAME"
} > "${DEST}/BUSCO_SUBMITTED.tmp.${JOB_TAG}"
mv -f -- \
  "${DEST}/BUSCO_SUBMITTED.tmp.${JOB_TAG}" \
  "${DEST}/BUSCO_SUBMITTED"
rm -f -- \
  "${DEST}/full_table.tsv.tmp" \
  "${DEST}/BUSCO_COMPLETE" \
  "${DEST}/BUSCO_FAILED"
printf '%s\n' "$FASTA" > "${DEST}/INPUT_SOURCE.txt"
printf '%s\n' "$LINEAGE_NAME" > "${DEST}/BUSCO_LINEAGE.txt"
printf '%s\n' "$JOB_TAG" > "${DEST}/LSF_JOB_ID.txt"
rm -rf -- "$OUT_DIR"

INPUT_FASTA="$FASTA"
if [[ "$FASTA" == *.gz ]]; then
  INPUT_FASTA="${TMPDIR}/${TOLID}.primary.curated.fa"
  gzip -cd -- "$FASTA" > "$INPUT_FASTA"
  [[ -s "$INPUT_FASTA" ]] || {
    echo "Failed to decompress input FASTA: $FASTA" >&2
    exit 1
  }
fi

busco \
  -i "$INPUT_FASTA" \
  -m genome \
  --miniprot \
  -l "$LINEAGE" \
  --out_path "$TMPDIR" \
  -o "$OUT_NAME" \
  --cpu "$CPU" \
  -f
