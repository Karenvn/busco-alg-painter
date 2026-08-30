#!/usr/bin/env bash

# Resolve curated primary assemblies with speciesops and submit fresh
# diptera_odb12 BUSCO jobs.
#
# Usage:
#   run_busco_diptera_batch.sh [--force] [tolids.txt]
#
# Completed diptera_odb12 ToLIDs are skipped unless --force is used.

set -euo pipefail

FORCE=0
LIST=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --force)
      FORCE=1
      ;;
    -h|--help)
      sed -n '3,9p' "$0"
      exit 0
      ;;
    -*)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
    *)
      if [[ -n "$LIST" ]]; then
        echo "Only one ToLID list may be supplied" >&2
        exit 2
      fi
      LIST="$1"
      ;;
  esac
  shift
done

LIST="${LIST:-/data/tol/users/kh18/nfs/tolids.txt}"
OUTPUT_BASE="${OUTPUT_BASE:-/data/tol/users/kh18/nfs/server_data/busco}"
LINEAGE="${LINEAGE:-/data/tol/resources/busco/latest/lineages/diptera_odb12}"
CPU="${CPU:-16}"
MEMORY_MB="${MEMORY_MB:-65536}"
QUEUE="${QUEUE:-normal}"
WALL_MINUTES="${WALL_MINUTES:-720}"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)"
WORKER="${SCRIPT_DIR}/run_busco_diptera_job.sh"
LINEAGE_NAME="$(basename "$LINEAGE")"

[[ -s "$LIST" ]] || {
  echo "ToLID list not found or empty: $LIST" >&2
  exit 2
}
[[ -x "$WORKER" ]] || {
  echo "BUSCO worker is missing or not executable: $WORKER" >&2
  exit 2
}
[[ -d "$LINEAGE" ]] || {
  echo "BUSCO lineage not found: $LINEAGE" >&2
  exit 2
}
[[ "$LINEAGE_NAME" == "diptera_odb12" ]] || {
  echo "Expected diptera_odb12 lineage, found: $LINEAGE_NAME" >&2
  exit 2
}
[[ "$CPU" =~ ^[1-9][0-9]*$ ]] || {
  echo "Invalid CPU count: $CPU" >&2
  exit 2
}
[[ "$MEMORY_MB" =~ ^[1-9][0-9]*$ ]] || {
  echo "Invalid memory request: $MEMORY_MB" >&2
  exit 2
}
[[ "$WALL_MINUTES" =~ ^[1-9][0-9]*$ ]] || {
  echo "Invalid wall time: $WALL_MINUTES" >&2
  exit 2
}

if { ! type module >/dev/null 2>&1 ||
     ! command -v bsub >/dev/null 2>&1 ||
     ! command -v bjobs >/dev/null 2>&1; } &&
   [[ -r /etc/profile ]]; then
  # Non-login SSH commands do not inherit the module/LSF setup on tol22.
  set +u
  # shellcheck disable=SC1091
  source /etc/profile
  set -u
fi

type module >/dev/null 2>&1 || {
  echo "The module command is not available" >&2
  exit 1
}
module load speciesops || {
  echo "Could not load the speciesops module" >&2
  exit 1
}
command -v speciesops >/dev/null 2>&1 || {
  echo "speciesops is not available after loading the module" >&2
  exit 1
}
command -v bsub >/dev/null 2>&1 || {
  echo "bsub is not available" >&2
  exit 1
}
command -v bjobs >/dev/null 2>&1 || {
  echo "bjobs is not available" >&2
  exit 1
}

mkdir -p -- "$OUTPUT_BASE"

RUNSTAMP="$(date +'%Y%m%d_%H%M%S')"
MANIFEST="${OUTPUT_BASE}/submission_manifest_${RUNSTAMP}.tsv"
JOBS="${OUTPUT_BASE}/submitted_jobs_${RUNSTAMP}.tsv"
SKIPPED="${OUTPUT_BASE}/skipped_tolids_${RUNSTAMP}.tsv"
ACCESSIONS="${OUTPUT_BASE}/tolids_accessions_${RUNSTAMP}.tsv"

printf 'tolid\trelease\tassembly_accession\tspecies_path\tcurated_dir\tinput_fasta\tdestination\n' > "$MANIFEST"
printf 'tolid\tjob_id\tbsub_response\n' > "$JOBS"
printf 'tolid\treason\n' > "$SKIPPED"
printf 'ToLID\tassembly_accession\n' > "$ACCESSIONS"

log() {
  printf '%s %s\n' "$(date +'%F %T')" "$*" >&2
}

warn() {
  printf '%s WARN: %s\n' "$(date +'%F %T')" "$*" >&2
}

trim_both() {
  sed -e 's/^[[:space:]][[:space:]]*//' \
      -e 's/[[:space:]][[:space:]]*$//' \
      -e 's/\r$//'
}

record_skip() {
  local tolid="$1"
  local reason="$2"
  warn "${tolid}: ${reason}"
  printf '%s\t%s\n' "$tolid" "$reason" >> "$SKIPPED"
}

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

marker_value() {
  local marker="$1"
  local key="$2"

  awk -F '\t' -v key="$key" '$1 == key {print $2; exit}' "$marker"
}

publish_accessions() {
  local canonical="${OUTPUT_BASE}/tolids_accessions.tsv"
  local temporary="${canonical}.tmp.$$"

  cp -f -- "$ACCESSIONS" "$temporary"
  mv -f -- "$temporary" "$canonical"
}

publish_tolid_list() {
  local canonical="${OUTPUT_BASE}/tolids.txt"
  local temporary="${canonical}.tmp.$$"

  cp -f -- "$LIST" "$temporary"
  mv -f -- "$temporary" "$canonical"
}

write_submitted_marker() {
  local destination="$1"
  local job_id="$2"
  local input_fasta="$3"
  local marker="${destination}/BUSCO_SUBMITTED"
  local temporary="${marker}.tmp.$$"

  {
    printf 'submitted_at\t%s\n' "$(date +'%Y-%m-%dT%H:%M:%S%z')"
    printf 'job_id\t%s\n' "$job_id"
    printf 'input\t%s\n' "$input_fasta"
    printf 'lineage\t%s\n' "$LINEAGE_NAME"
  } > "$temporary"
  mv -f -- "$temporary" "$marker"
}

lsf_job_is_active() {
  local job_id="$1"
  local status=""

  status="$(bjobs -a -noheader -o stat "$job_id" 2>/dev/null |
    awk 'NR == 1 {print $1}' || true)"
  case "$status" in
    PEND|RUN|WAIT|PSUSP|USUSP|SSUSP|PROV|UNKWN|ZOMBI)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

get_species_path() {
  local tolid="$1"
  local out=""
  local path=""
  local attempt

  for attempt in 1 2 3; do
    out="$(LC_ALL=C speciesops getdir --tolid "$tolid" 2>&1 || true)"
    path="$(printf '%s\n' "$out" |
      sed -n 's/^[[:space:]]*Species[[:space:]][[:space:]]*directory[[:space:]][[:space:]]*path:[[:space:]]*\(\/[^[:space:]][^[:space:]]*\)[[:space:]]*$/\1/p' |
      head -n 1)"
    if [[ -z "$path" ]]; then
      path="$(printf '%s\n' "$out" |
        grep -Eo '^/[^[:space:]]+' |
        head -n 1 || true)"
    fi
    if [[ -n "$path" && -d "$path" ]]; then
      printf '%s\n' "$path"
      return 0
    fi
    sleep 1
  done

  printf '%s\n' "$out" > "${OUTPUT_BASE}/speciesops_${tolid}_${RUNSTAMP}.txt"
  return 1
}

latest_primary_release_name() {
  local species_path="$1"
  local tolid="$2"
  local release_dir="${species_path%/}/assembly/release"
  local hit=""

  [[ -d "$release_dir" ]] || return 1

  hit="$(find "$release_dir" -maxdepth 1 \( -type l -o -type d \) \
    \( -name "${tolid}.hap1.[0-9]*" -o -name "${tolid}.[0-9]*" \) \
    ! -name '*_alternate_haplotype*' \
    ! -name "${tolid}.hap2.*" \
    -print 2>/dev/null |
    sort -V |
    tail -n 1)"

  [[ -n "$hit" ]] || return 1
  basename "$hit"
}

release_accession() {
  local species_path="$1"
  local release_name="$2"
  local release_link="${species_path%/}/assembly/release/${release_name}"
  local target=""
  local accession=""

  [[ -e "$release_link" || -L "$release_link" ]] || return 1

  target="$(readlink "$release_link" 2>/dev/null || true)"
  accession="$(basename "$target")"
  if [[ "$accession" =~ ^GC[AF]_[0-9]+\.[0-9]+$ ]]; then
    printf '%s\n' "$accession"
    return 0
  fi

  target="$(readlink -f "$release_link" 2>/dev/null || true)"
  accession="$(basename "$target")"
  if [[ "$accession" =~ ^GC[AF]_[0-9]+\.[0-9]+$ ]]; then
    printf '%s\n' "$accession"
    return 0
  fi

  return 1
}

curated_dir_matching_release() {
  local species_path="$1"
  local tolid="$2"
  local release_name="$3"
  local curated_root="${species_path%/}/assembly/curated"
  local candidate=""
  local base=""
  local version=""

  [[ -d "$curated_root" ]] || return 1

  candidate="${curated_root}/${release_name}"
  if [[ -d "$candidate" ]]; then
    printf '%s\n' "$candidate"
    return 0
  fi

  [[ "$release_name" =~ \.([0-9]+)$ ]] || return 1
  version="${BASH_REMATCH[1]}"

  while IFS= read -r candidate; do
    [[ -d "$candidate" ]] || continue
    base="$(basename "$candidate")"
    [[ "$base" == "${tolid}.hap2."* ]] && continue
    [[ "$base" == *"_alternate_haplotype" ]] && continue
    [[ "$base" =~ \.([0-9]+)$ ]] || continue
    [[ "${BASH_REMATCH[1]}" == "$version" ]] || continue
    printf '%s\n' "$candidate"
    return 0
  done < <(
    compgen -G "${curated_root}/${tolid}.*" 2>/dev/null |
      sort -V
  )

  return 1
}

latest_curated_dir() {
  local species_path="$1"
  local tolid="$2"
  local candidate=""
  local base=""
  local hit=""

  while IFS= read -r candidate; do
    [[ -d "$candidate" ]] || continue
    base="$(basename "$candidate")"
    [[ "$base" == "${tolid}.hap2."* ]] && continue
    [[ "$base" == *"_alternate_haplotype" ]] && continue
    hit="$candidate"
  done < <(
    compgen -G "${species_path%/}/assembly/curated/${tolid}.*" 2>/dev/null |
      sort -V
  )

  [[ -n "$hit" ]] || return 1
  printf '%s\n' "$hit"
}

select_curated_dir() {
  local species_path="$1"
  local tolid="$2"
  local release_name=""
  local curated_dir=""

  release_name="$(latest_primary_release_name "$species_path" "$tolid" || true)"
  if [[ -n "$release_name" ]]; then
    curated_dir="$(
      curated_dir_matching_release "$species_path" "$tolid" "$release_name" ||
        true
    )"
    if [[ -n "$curated_dir" ]]; then
      printf '%s\n' "$curated_dir"
      return 0
    fi
    return 1
  fi

  latest_curated_dir "$species_path" "$tolid"
}

choose_primary_curated_fasta() {
  local curated_dir="$1"
  local tolid="$2"
  local base=""
  local version=""
  local pattern=""
  local hit=""
  local -a patterns=()

  base="$(basename "$curated_dir")"
  if [[ "$base" =~ \.([0-9]+)$ ]]; then
    version="${BASH_REMATCH[1]}"
  fi

  if [[ -n "$version" ]]; then
    patterns+=(
      "${curated_dir}/${tolid}.hap1.${version}.primary.curated.fa"
      "${curated_dir}/${tolid}.hap1.${version}.primary.curated.fa.gz"
      "${curated_dir}/${tolid}.${version}.primary.curated.fa"
      "${curated_dir}/${tolid}.${version}.primary.curated.fa.gz"
    )
  fi
  patterns+=(
    "${curated_dir}/${base}.primary.curated.fa"
    "${curated_dir}/${base}.primary.curated.fa.gz"
    "${curated_dir}/${tolid}.hap1.*.primary.curated.fa"
    "${curated_dir}/${tolid}.hap1.*.primary.curated.fa.gz"
    "${curated_dir}/${tolid}.[0-9]*.primary.curated.fa"
    "${curated_dir}/${tolid}.[0-9]*.primary.curated.fa.gz"
  )

  for pattern in "${patterns[@]}"; do
    hit="$(compgen -G "$pattern" 2>/dev/null |
      grep -v -E 'hap2|all_haplotigs|alternate_haplotype|pre_curation' |
      sort -V |
      tail -n 1 || true)"
    if [[ -n "$hit" && -s "$hit" ]]; then
      printf '%s\n' "$hit"
      return 0
    fi
  done

  return 1
}

seen=0
resolved=0
submitted=0
skipped=0
failed=0

while IFS= read -r raw || [[ -n "$raw" ]]; do
  TOLID="$(printf '%s' "$raw" | trim_both)"
  [[ -z "$TOLID" || "$TOLID" =~ ^# ]] && continue
  seen=$((seen + 1))

  if [[ ! "$TOLID" =~ ^[A-Za-z0-9._-]+$ ]]; then
    record_skip "$TOLID" "invalid ToLID"
    skipped=$((skipped + 1))
    continue
  fi

  SPECIES_PATH="$(get_species_path "$TOLID" || true)"
  if [[ -z "$SPECIES_PATH" || ! -d "$SPECIES_PATH" ]]; then
    record_skip "$TOLID" "speciesops directory could not be resolved"
    skipped=$((skipped + 1))
    continue
  fi

  RELEASE_NAME="$(latest_primary_release_name "$SPECIES_PATH" "$TOLID" || true)"
  if [[ -z "$RELEASE_NAME" ]]; then
    record_skip "$TOLID" "current primary release could not be resolved"
    skipped=$((skipped + 1))
    continue
  fi

  ACCESSION="$(release_accession "$SPECIES_PATH" "$RELEASE_NAME" || true)"
  if [[ -z "$ACCESSION" ]]; then
    record_skip "$TOLID" "assembly accession could not be resolved from ${RELEASE_NAME}"
    skipped=$((skipped + 1))
    continue
  fi

  CURATED_DIR="$(select_curated_dir "$SPECIES_PATH" "$TOLID" || true)"
  if [[ -z "$CURATED_DIR" || ! -d "$CURATED_DIR" ]]; then
    record_skip "$TOLID" "current curated assembly directory not found"
    skipped=$((skipped + 1))
    continue
  fi

  FASTA="$(choose_primary_curated_fasta "$CURATED_DIR" "$TOLID" || true)"
  if [[ -z "$FASTA" || ! -s "$FASTA" ]]; then
    record_skip "$TOLID" "primary curated FASTA not found in ${CURATED_DIR}"
    skipped=$((skipped + 1))
    continue
  fi

  DEST="${OUTPUT_BASE}/${TOLID}"
  mkdir -p -- "$DEST"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$TOLID" "$RELEASE_NAME" "$ACCESSION" "$SPECIES_PATH" \
    "$CURATED_DIR" "$FASTA" "$DEST" >> "$MANIFEST"
  printf '%s\t%s\n' "$TOLID" "$ACCESSION" >> "$ACCESSIONS"
  resolved=$((resolved + 1))

  SUBMITTED_MARKER="${DEST}/BUSCO_SUBMITTED"
  if [[ -s "$SUBMITTED_MARKER" ]]; then
    ACTIVE_JOB_ID="$(marker_value "$SUBMITTED_MARKER" job_id || true)"
    if [[ "$ACTIVE_JOB_ID" =~ ^[0-9]+$ ]] &&
       lsf_job_is_active "$ACTIVE_JOB_ID"; then
      record_skip "$TOLID" "BUSCO job ${ACTIVE_JOB_ID} is already queued or running"
      skipped=$((skipped + 1))
      continue
    fi
    warn "${TOLID}: ignoring stale BUSCO_SUBMITTED marker"
    if ! rm -f -- "$SUBMITTED_MARKER"; then
      warn "${TOLID}: stale BUSCO_SUBMITTED marker could not be removed"
    fi
  fi

  if [[ $FORCE -eq 0 && -s "${DEST}/full_table.tsv" ]]; then
    EXISTING_LINEAGE="$(busco_table_lineage "${DEST}/full_table.tsv" || true)"
    COMPLETE_MARKER="${DEST}/BUSCO_COMPLETE"
    COMPLETE_INPUT=""
    COMPLETE_LINEAGE=""
    if [[ -s "$COMPLETE_MARKER" ]]; then
      COMPLETE_INPUT="$(marker_value "$COMPLETE_MARKER" input || true)"
      COMPLETE_LINEAGE="$(marker_value "$COMPLETE_MARKER" lineage || true)"
    fi

    if [[ -s "${DEST}/BUSCO_FAILED" ]]; then
      warn "${TOLID}: retrying because the previous BUSCO attempt failed"
    elif [[ "$EXISTING_LINEAGE" != "$LINEAGE_NAME" ]]; then
      if [[ -n "$EXISTING_LINEAGE" ]]; then
        warn "${TOLID}: replacing existing ${EXISTING_LINEAGE} output with ${LINEAGE_NAME}"
      else
        warn "${TOLID}: existing BUSCO table has no readable lineage; replacing it"
      fi
    elif [[ -z "$COMPLETE_INPUT" || -z "$COMPLETE_LINEAGE" ]]; then
      warn "${TOLID}: existing BUSCO table has no complete-run metadata; replacing it"
    elif [[ "$COMPLETE_LINEAGE" != "$LINEAGE_NAME" ]]; then
      warn "${TOLID}: completion marker lineage is ${COMPLETE_LINEAGE}; replacing it"
    elif [[ "$COMPLETE_INPUT" != "$FASTA" ]]; then
      warn "${TOLID}: current curated FASTA differs from the completed BUSCO input; replacing it"
    else
      record_skip "$TOLID" "complete ${LINEAGE_NAME} output already exists for the current FASTA"
      skipped=$((skipped + 1))
      continue
    fi
  fi

  BSUB_ARGS=(
    -J "busco_${TOLID}"
    -q "$QUEUE"
    -n "$CPU"
    -W "$WALL_MINUTES"
    -M "$MEMORY_MB"
    -R "select[mem>${MEMORY_MB}] rusage[mem=${MEMORY_MB}] span[hosts=1]"
    -cwd "$DEST"
    -o "${DEST}/BUSCO_${TOLID}.%J.out"
    -e "${DEST}/BUSCO_${TOLID}.%J.err"
    bash "$WORKER" "$TOLID" "$FASTA" "$DEST" "$LINEAGE" "$CPU"
  )

  if RESPONSE="$(bsub "${BSUB_ARGS[@]}" 2>&1)"; then
    JOB_ID="$(printf '%s\n' "$RESPONSE" |
      sed -n 's/^Job <\([0-9][0-9]*\)>.*/\1/p' |
      head -n 1)"
    printf '%s\t%s\t%s\n' "$TOLID" "${JOB_ID:-unknown}" "$RESPONSE" >> "$JOBS"
    if ! write_submitted_marker "$DEST" "${JOB_ID:-unknown}" "$FASTA"; then
      warn "${TOLID}: job submitted, but BUSCO_SUBMITTED marker could not be written"
    fi
    log "${TOLID}: ${RESPONSE}"
    submitted=$((submitted + 1))
  else
    RC=$?
    record_skip "$TOLID" "bsub failed: ${RESPONSE//$'\t'/ }"
    failed=$((failed + 1))
  fi
done < "$LIST"

publish_tolid_list
publish_accessions

log "seen=${seen} resolved=${resolved} submitted=${submitted} skipped=${skipped} failed=${failed}"
log "manifest: ${MANIFEST}"
log "jobs: ${JOBS}"
log "skipped: ${SKIPPED}"
log "accessions: ${ACCESSIONS}"

[[ $failed -eq 0 ]]
