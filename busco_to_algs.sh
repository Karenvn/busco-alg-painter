#!/bin/bash

# Batch BUSCO ALG painting for ToLID/accession tables.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="${DATA_ROOT:-.}"
BUSCO_DIR="${BUSCO_DIR:-${DATA_ROOT}/busco}"
OUTPUT_DIR="${OUTPUT_DIR:-${DATA_ROOT}/algs}"
PROFILE="${PROFILE:-${LINEAGE:-auto}}"
ALG_REF="${ALG_REF:-}"
CUSTOM_CONFIG="${CUSTOM_CONFIG:-}"
TOLID_FILE="${TOLID_FILE:-}"
LABEL_WINDOW_MB="${LABEL_WINDOW_MB:-0}"
LABEL_WINDOW_MIN_BUSCOS="${LABEL_WINDOW_MIN_BUSCOS:-5}"
LABEL_WINDOW_MIN_FRACTION="${LABEL_WINDOW_MIN_FRACTION:-0.5}"

ACCESSION_TABLE_CANDIDATES=(
  "tolid_accessions.tsv"
  "tolid_accession.tsv"
  "tolids_accession.tsv"
  "tolids_accessions.tsv"
  "${SCRIPT_DIR}/tolid_accessions.tsv"
  "${SCRIPT_DIR}/tolid_accession.tsv"
  "${SCRIPT_DIR}/tolids_accession.tsv"
  "${SCRIPT_DIR}/tolids_accessions.tsv"
)

if command -v busco-alg-painter >/dev/null 2>&1; then
  RUN_CMD=(busco-alg-painter run)
else
  export PYTHONPATH="${SCRIPT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}"
  RUN_CMD=(python3 -m busco_alg_painter run)
fi

resolve_accession_file() {
  if [[ -n "${ACCESSION_FILE:-}" ]]; then
    printf '%s\n' "$ACCESSION_FILE"
    return 0
  fi

  local candidate
  for candidate in "${ACCESSION_TABLE_CANDIDATES[@]}"; do
    if [[ -f "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

ACCESSION_FILE="$(resolve_accession_file)" || {
  echo "ERROR: no accession table found; set ACCESSION_FILE"
  exit 1
}

if [[ -n "$TOLID_FILE" && ! -f "$TOLID_FILE" ]]; then
  echo "ERROR: TOLID_FILE not found: $TOLID_FILE"
  exit 1
fi
if [[ -n "$ALG_REF" && ! -f "$ALG_REF" ]]; then
  echo "ERROR: ALG_REF not found: $ALG_REF"
  exit 1
fi
if [[ -n "$CUSTOM_CONFIG" && ! -f "$CUSTOM_CONFIG" ]]; then
  echo "ERROR: CUSTOM_CONFIG not found: $CUSTOM_CONFIG"
  exit 1
fi

TOLID_SOURCE="${TOLID_FILE:-$ACCESSION_FILE}"
mkdir -p "$OUTPUT_DIR" || exit 1

get_accession() {
  local tolid="$1"
  awk -F '\t' -v tolid="$tolid" '$1 == tolid {print $2; exit}' "$ACCESSION_FILE"
}

is_header_tolid() {
  local lower
  lower="$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')"
  [[ "$lower" == "tolid" || "$lower" == "tol_id" || "$lower" == "sample" ]]
}

process_tolid() {
  local tolid="$1"
  local accession busco_input tolid_output
  local -a args

  accession="$(get_accession "$tolid")"
  busco_input="${BUSCO_DIR}/${tolid}/full_table.tsv"
  tolid_output="${OUTPUT_DIR}/${tolid}"

  if [[ -z "$accession" ]]; then
    echo "WARN: no accession for $tolid; skipping"
    return 2
  fi
  if [[ ! -f "$busco_input" ]]; then
    echo "WARN: missing $busco_input; skipping"
    return 2
  fi

  mkdir -p "$tolid_output" || return 1
  args=(
    --query-table "$busco_input"
    --prefix "${tolid_output}/"
    --profile "$PROFILE"
    --accession "$accession"
    --write-summary
    --label-window-mb "$LABEL_WINDOW_MB"
    --label-window-min-buscos "$LABEL_WINDOW_MIN_BUSCOS"
    --label-window-min-fraction "$LABEL_WINDOW_MIN_FRACTION"
  )
  if [[ -n "$ALG_REF" ]]; then
    args+=(--reference-table "$ALG_REF")
  fi
  if [[ -n "$CUSTOM_CONFIG" ]]; then
    args+=(--config "$CUSTOM_CONFIG")
  fi

  echo "Processing $tolid ($accession)"
  "${RUN_CMD[@]}" "${args[@]}"
}

total=0
success=0
failed=0
skipped=0

while IFS= read -r line || [[ -n "$line" ]]; do
  tolid="${line%%$'\t'*}"
  [[ -z "$tolid" || "$tolid" =~ ^# ]] && continue
  is_header_tolid "$tolid" && continue

  total=$((total + 1))
  process_tolid "$tolid"
  status=$?
  if [[ $status -eq 0 ]]; then
    success=$((success + 1))
  elif [[ $status -eq 2 ]]; then
    skipped=$((skipped + 1))
  else
    failed=$((failed + 1))
  fi
done < "$TOLID_SOURCE"

echo "Batch complete: total=$total success=$success skipped=$skipped failed=$failed"
echo "Output: $OUTPUT_DIR"
[[ $failed -eq 0 ]]
