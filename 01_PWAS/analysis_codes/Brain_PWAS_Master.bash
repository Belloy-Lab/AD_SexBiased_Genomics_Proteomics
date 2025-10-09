#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# Master PWAS job submitter (LSF bsub)
# One row in CSV = one job submission
# -------------------------------
# Logs root (jobs will auto-organize per cohort here)
LOG_ROOT="/storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS"

# ---- ARGUMENTS ----
CSV_FILE="${1:-}"
DRY_RUN="${2:-}"

if [[ -z "${CSV_FILE}" || ! -f "${CSV_FILE}" ]]; then
  echo "Usage: $0 <pwas_jobs.csv> [--dry-run]"
  exit 1
fi

if [[ "${DRY_RUN:-}" == "--dry-run" ]]; then
  echo ">>> DRY-RUN MODE: commands will be printed but not submitted"
fi

# ---- SUBMIT LOOP ----
# CSV columns:
# cohort,ss_dir,ss,sex,out_dir,out_file
# Skip header, robust IFS, allow spaces in paths if quoted in CSV.
tail -n +2 "${CSV_FILE}" | while IFS=',' read -r cohort ss_dir ss sex out_dir out_file; do
  # Trim potential whitespace
  cohort="$(echo "$cohort" | xargs)"
  ss_dir="$(echo "$ss_dir" | xargs)"
  ss="$(echo "$ss" | xargs)"
  sex="$(echo "$sex" | xargs)"
  out_dir="$(echo "$out_dir" | xargs)"
  out_file="$(echo "$out_file" | xargs)"

  # Expand $USER in CSV paths (necessary so CSV with $USER works)
  ss_dir="${ss_dir//\$USER/$USER}"
  out_dir="${out_dir//\$USER/$USER}"

  # Per-job logs
  LOG_DIR="${LOG_ROOT}/${cohort}/logs"
  mkdir -p "${LOG_DIR}"

  # Job name
  JOB_NAME="PWAS_${cohort}_${sex}_Brain"

  # Output existence check (skip if already completed output exists)
  # Adjust this condition to your pipeline's final output(s) if needed:
  if [[ -f "${out_dir}/${out_file}.fusion.z" || -f "${out_dir}/${out_file}.done" ]]; then
    echo "SKIP: ${JOB_NAME} seems done (found output in ${out_dir})."
    continue
  fi

  # Build the bsub command
  CMD=(bsub
    -q subscription \
    -J "${JOB_NAME}" \
    -n 3 \
    -o "${LOG_DIR}/${JOB_NAME}.%J.out" \
    -e "${LOG_DIR}/${JOB_NAME}.%J.err" \
    -R "span[hosts=1] rusage[mem=100000]" \
    -a docker(satheshsiva27/multiomics-toolkit:0.1)
    # NOTE: Update the script path below if the location changes
    bash PWAS/analysis_codes/Brain_pwas_wrapper.bash
      --ss_dir "${ss_dir}"
      --ss "${ss}"
      --ss_sex "${sex}"
      --pw_m_dir PWAS/Wingoetal2023_protein_weights/Brain/Males # <-- FIX the male brin protein weights path
      --pw_m train_weights.pos # <-- fix the file name
      --pw_f_dir PWAS/Wingoetal2023_protein_weights/Brain/Females   # <-- FIX the female brin protein weights path
      --pw_f train_weights.pos # <-- fix the file name
      --pw_b_dir PWAS/Wingoetal2023_protein_weights/Brain/joint # <-- FIX the joint (both) brin protein weights path
      --pw_b train_weights.pos # <-- fix the file name
      --ref_ld_chr 1000G_Grch37/1KG_EUR/file_name_by_chr(1-22) # <-- FIX 1000G LD ref panel path
      --out_dir "${out_dir}"
      --out_file "${out_file}"
  )

  # Print and optionally submit
  echo "Submitting: ${JOB_NAME}"
  printf '  %q ' "${CMD[@]}"; echo

  if [[ "${DRY_RUN:-}" != "--dry-run" ]]; then
    "${CMD[@]}"
    # tiny jitter to avoid log stamp collisions in a large burst
    sleep 0.2
  fi
done
