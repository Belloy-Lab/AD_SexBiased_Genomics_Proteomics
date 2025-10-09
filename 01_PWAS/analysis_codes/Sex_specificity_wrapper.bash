#!/usr/bin/env bash
set -euo pipefail

# ----- Paths & resources (edit as needed) -----
R_SCRIPT="sex_specific_genes.R"

BRAIN_MAP="SupportingFiles/Hg38_Biomart07082025.txt"
CSF_MAP="SupportingFiles/CSF_SS_library.txt"

LOG_DIR="${LOG_DIR:-$PWD/logs}"
mkdir -p "$LOG_DIR"

QUEUE="${QUEUE:-subscription}"
N_CORES="${N_CORES:-3}"
MEM_MB="${MEM_MB:-100000}"  # 100 GB
RSTR="span[hosts=1] rusage[mem=${MEM_MB}]"

# Optional LSF accounting (customize or export before running)
LSF_PROJECT="${LSF_PROJECT:-compute-belloy}"
LSF_TIER="${LSF_TIER:-t1}"
JOB_GROUP="${JOB_GROUP:-"/$USER/${LSF_PROJECT}"}"
ACCOUNT_GROUP="${ACCOUNT_GROUP:-"${LSF_PROJECT}-${LSF_TIER}"}"
SERVICE_CLASS="${SERVICE_CLASS:-"${LSF_PROJECT}-${LSF_TIER}"}"

# Build common bsub flags
BSUB_BASE=( -q "$QUEUE" -n "$N_CORES" -M "$MEM_MB" -R "$RSTR" )
BSUB_GRP=( -g "$JOB_GROUP" -G "$ACCOUNT_GROUP" -sla "$SERVICE_CLASS" )

# ----- Job matrix -----
# Cohort -> work_dir and tag used in prefix
declare -A WORKDIR=(
  [EU_all]="/PWAS/EU_all"
  [EU_noUKB]="/PWAS/EU_noUKB"
  [AFR]="/PWAS/AFR"
)
declare -A TAG=(
  [EU_all]="EUall"
  [EU_noUKB]="EUnoUKB"
  [AFR]="AFR"
)

TISSUES=(Brain CSF)
SEXES=(Male Female)

# ----- Submit all combinations -----
for cohort in EU_all EU_noUKB AFR; do
  for tissue in "${TISSUES[@]}"; do
    for sex in "${SEXES[@]}"; do
      prefix="${TAG[$cohort]}_${sex}_${tissue}"
      work_dir="${WORKDIR[$cohort]}"

      if [[ "$tissue" == "Brain" ]]; then
        MAP_FLAG="--brain_map"
        MAP_PATH="$BRAIN_MAP"
      else
        MAP_FLAG="--csf_map"
        MAP_PATH="$CSF_MAP"
      fi

      JOB_NAME="SSG_V6_${cohort}_${tissue}_${sex}"

      bsub \
        "${BSUB_GRP[@]}" \
        "${BSUB_BASE[@]}" \
        -J "$JOB_NAME" \
        -o "${LOG_DIR}/${JOB_NAME}.%J.out" \
        -e "${LOG_DIR}/${JOB_NAME}.%J.err" \
        Rscript "$R_SCRIPT" \
          --work_dir "$work_dir" \
          --tissue "$tissue" \
          --sex "$sex" \
          --prefix "$prefix" \
          "$MAP_FLAG" "$MAP_PATH"

    done
  done
done
