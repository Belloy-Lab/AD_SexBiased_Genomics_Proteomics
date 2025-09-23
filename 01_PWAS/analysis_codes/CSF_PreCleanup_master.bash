#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# CSF Pre-Cleanup job submitter (LSF bsub)
# One row in CSV = one job submission
# -------------------------------
# NOTE: Update the script path below if the location changes
LOG_ROOT="PWAS"

CSV_FILE="${1:-}"
DRY_RUN="${2:-}"

if [[ -z "${CSV_FILE}" || ! -f "${CSV_FILE}" ]]; then
  echo "Usage: $0 CSF_PreCleanup_input.csv [--dry-run]"
  exit 1
fi

if [[ "${DRY_RUN:-}" == "--dry-run" ]]; then
  echo ">>> DRY-RUN MODE: commands will be printed but not submitted"
fi

# CSV columns:
# cohort_dir,cohort_tag,sex,in_dir,in_GWAS,out_GWAS,Nsize_GWAS,out_dir
tail -n +2 "${CSV_FILE}" | while IFS=',' read -r cohort_dir cohort_tag sex in_dir in_GWAS out_GWAS Nsize_GWAS out_dir; do
  cohort_dir="$(echo "$cohort_dir" | xargs)"
  cohort_tag="$(echo "$cohort_tag" | xargs)"
  sex="$(echo "$sex" | xargs)"
  in_dir="$(echo "$in_dir" | xargs)"
  in_GWAS="$(echo "$in_GWAS" | xargs)"
  out_GWAS="$(echo "$out_GWAS" | xargs)"
  Nsize_GWAS="$(echo "$Nsize_GWAS" | xargs)"
  out_dir="$(echo "$out_dir" | xargs)"

  # Expand $USER if present in CSV
  in_dir="${in_dir//\$USER/$USER}"
  out_dir="${out_dir//\$USER/$USER}"

  # Logs
  LOG_DIR="${LOG_ROOT}/${cohort_dir}/logs"
  mkdir -p "${LOG_DIR}"

  # Names: Job uses capitalized sex; logs use lowercase sex
  SexCap="${sex^}"
  sex_lc="${sex,,}"
  JOB_NAME="CSF_PreClean_${SexCap}_${cohort_tag}"
  LOG_PREFIX="CSF_PreClean_${sex_lc}_${cohort_tag}"

  CMD=(bsub
    -g /$USER/compute-belloy
    -J "${JOB_NAME}"
    -n 1
    -N
    -o "${LOG_DIR}/${LOG_PREFIX}.%J.out"
    -e "${LOG_DIR}/${LOG_PREFIX}.%J.err"
    -q subscription
    -R "rusage[mem=40GB] span[hosts=1]"
    -G compute-belloy-t1
    -sla compute-belloy-t1
    -a docker(dmr07083/fusion-project:4.3.2)
    # NOTE: Update the script path below if the location changes
    bash PWAS/analysis_codes/CSF_PreCleanup.bash
      --in_dir "${in_dir}"
      --in_GWAS "${in_GWAS}"
      --out_GWAS "${out_GWAS}"
      --Nsize_GWAS "${Nsize_GWAS}"
      --out_dir "${out_dir}"
  )

  echo "Submitting: ${JOB_NAME}"
  printf '  %q ' "${CMD[@]}"; echo

  if [[ "${DRY_RUN:-}" != "--dry-run" ]]; then
    "${CMD[@]}"
    sleep 0.2
  fi
done
