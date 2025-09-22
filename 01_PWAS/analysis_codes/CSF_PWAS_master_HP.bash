#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# Master CSF PWAS job submitter (LSF bsub)
# One row in CSV = one job submission
# -------------------------------
# Logs root (jobs will auto-organize per cohort here)
LOG_ROOT="/storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS"

# ---- ARGUMENTS ----
CSV_FILE="${1:-}"
DRY_RUN="${2:-}"

if [[ -z "${CSV_FILE}" || ! -f "${CSV_FILE}" ]]; then
  echo "Usage: $0 <csf_pwas_jobs.csv> [--dry-run]"
  exit 1
fi

if [[ "${DRY_RUN:-}" == "--dry-run" ]]; then
  echo ">>> DRY-RUN MODE: commands will be printed but not submitted"
fi

# ---- SUBMIT LOOP ----
# CSV columns:
# cohort,ss_dir,ssx,ssn,sex,out_dir,out_file
# Skip header, robust IFS, allow spaces in paths if quoted in CSV.
tail -n +2 "${CSV_FILE}" | while IFS=',' read -r cohort ss_dir ssx ssn sex out_dir out_file; do
  # Trim potential whitespace
  cohort="$(echo "$cohort" | xargs)"
  ss_dir="$(echo "$ss_dir" | xargs)"
  ssx="$(echo "$ssx" | xargs)"
  ssn="$(echo "$ssn" | xargs)"
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
  JOB_NAME="PWAS_${cohort}_${sex}_CSF"

  # Output existence check (skip if already completed output exists)
  # Adjust this condition to your pipeline's final output(s) if needed:
  if [[ -f "${out_dir}/${out_file}.fusion.z" || -f "${out_dir}/${out_file}.done" ]]; then
    echo "SKIP: ${JOB_NAME} seems done (found output in ${out_dir})."
    continue
  fi

  # Build the bsub command
  CMD=(bsub
    -g /$USER/compute-belloy
    -q subscription
    -J "${JOB_NAME}"
    -n 3
    -o "${LOG_DIR}/${JOB_NAME}.%J.out"
    -e "${LOG_DIR}/${JOB_NAME}.%J.err"
    -R "span[hosts=1] rusage[mem=150000]"
    -M 150000
    -G compute-belloy-t1
    -sla compute-belloy-t1
    -a docker(dmr07083/fusion-project:4.3.2)
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/CSF_pwas.bash
      --ss_dir "${ss_dir}"
      --ssx "${ssx}"
      --ssn "${ssn}"
      --ss_sex "${sex}"
      --pw_m_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt
      --pw_m NGI_CSF_male_cis_HP-weights.pos
      --pw_f_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt
      --pw_f NGI_CSF_female_cis_HP-weights.pos
      --pw_b_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_pQTL_files/PWAS_Weights
      --pw_b NGI_CSF_cis_HP-weights.pos
      --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/LD/LD/CSF_proteomics_3107_samples_maf_0032_geno_1_chr_
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
