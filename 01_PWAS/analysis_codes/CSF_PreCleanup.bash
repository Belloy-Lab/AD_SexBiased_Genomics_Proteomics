#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# Master pipeline: GWAS → liftover/refpanel (R) → munge (bash) → clean (R)
# Uses GNU getopt (supports long options)
# -------------------------------

# ----- parse long options -----
OPTS=$(
  getopt -o '' \
    -l in_dir:,in_GWAS:,out_dir:, \
    -l Nsize_GWAS:,out_GWAS:,help \
    -n "$0" -- "$@"
)
eval set -- "$OPTS"

IN_DIR=""; IN_GWAS=""; OUT_DIR=""; OUT_GWAS=""; NSIZE_GWAS="";

while true; do
  case "$1" in
    --in_dir)        IN_DIR="$2"; shift 2 ;;
    --in_GWAS)       IN_GWAS="$2"; shift 2 ;;
    --out_GWAS)      OUT_GWAS="$2"; shift 2 ;;
    --Nsize_GWAS)    NSIZE_GWAS="$2"; shift 2 ;;
    --out_dir)       OUT_DIR="$2"; shift 2 ;;
    --help) cat <<'EOF'

Usage:
  CSF_PreCleanup_Master.sh \
    --in_dir    /path/to/GWAS_summary_stats_dir \
    --in_GWAS   <gwas_input_filename> \
    --out_GWAS  <output_stem> \
    --Nsize_GWAS <integer_sample_size> \
    --out_dir   /path/to/output_dir

Required arguments:
  --in_dir       Directory containing the GWAS summary stats file (sex-specific or unisex).
                 The script will read: <in_dir>/<in_GWAS>

  --in_GWAS      GWAS summary stats filename (no path). Must contain columns:
                 CHR, BP, ALLELE1 (A1), ALLELE0 (A2), SNP (optional), P, A1FREQ (for munge).

  --out_GWAS     Output stem (no extension). Used to name all downstream outputs:
                 <out_dir>/<out_GWAS>.CSFsexstrat.* and <out_dir>/<out_GWAS>.CSFnonsex.* 

--Nsize_GWAS   Integer sample size used by LDSC munge_sumstats.py (--N).

  --out_dir      Output directory; will be created if missing.

What this script does:
  Step 1 (R): Intersects <in_GWAS> with two CSF LD reference panels, producing:
      <out_dir>/<out_GWAS>.CSFsexstrat.txt
      <out_dir>/<out_GWAS>.CSFnonsex.txt
  Step 2 (LDSC): Munges each intersected file, producing:
      <out_dir>/<out_GWAS>.CSFsexstrat.sumstats.gz
      <out_dir>/<out_GWAS>.CSFnonsex.sumstats.gz
      (then unzips to .sumstats)
  Step 3 (R): Final cleanup of each .sumstats file in <out_dir>.

Example:
  CSF_PreCleanup_Master.sh \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update \
    --in_GWAS ADGC_ADSP_UKB_FinnGen_Females_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var \
    --out_GWAS AD_Females_cc_rb.gen090.noAPOE.shared_var \
    --Nsize_GWAS 636154 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/CSF

Notes:
  - Input paths are resolved as: <in_dir>/<in_GWAS>
  - Outputs are written to: <out_dir>
  - This script assumes a working 'ldsc' conda environment is available and activated in Step 2.
EOF
      exit 0 ;;
    --) shift; break ;;
    *) echo "Internal getopt error"; exit 2 ;;
  esac
done

# ----- validate required values are provided -----
[[ -z "$IN_DIR"  || -z "$IN_GWAS" || -z "$OUT_DIR" ]] && {
  echo "ERROR: --in_dir, --in_GWAS, --out_dir are required"; exit 2; }
[[ -z "$OUT_GWAS" || -z "$NSIZE_GWAS" ]] && {
  echo "ERROR: --out_GWAS, --Nsize_GWAS are required"; exit 2; }

# ----- existence/type checks -----
[[ -d "$IN_DIR" ]] || { echo "ERROR: --in_dir not found: $IN_DIR"; exit 3; }
[[ -s "$IN_DIR/$IN_GWAS" ]] || { echo "ERROR: GWAS input not found: $IN_DIR/$IN_GWAS"; exit 3; }
[[ "$NSIZE_GWAS" =~ ^[0-9]+$ ]] || { echo "ERROR: --Nsize_GWAS must be integer"; exit 3; }

# helper: time block
timeit() {
  local label="$1"; shift
  echo "==> [$label] started: $(date)"
  local t0
  t0=$(date +%s)
  eval "$@"
  local t1
  t1=$(date +%s)
  echo "==> [$label] finished: $(date) | elapsed: $((t1-t0))s)"
}

# ==============================
# STEP 1: Liftover + EUR LD intersections (R)
# ==============================
mkdir -p "$OUT_DIR"

# GWAS summary stat sex-strat ref panel intersect.
Rscript /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/LDintersect_CSF.R \
  --in_dir "${IN_DIR}/" \
  --gwas_file "$IN_GWAS" \
  --out_dir "${OUT_DIR}/" \
  --out_file "${OUT_GWAS}.CSFsexstrat.txt" \
  --LD_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/LD/LD/ \
  --LD_file chr1-22.txt

Rscript /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/LDintersect_CSF.R \
  --in_dir "${IN_DIR}/" \
  --gwas_file "$IN_GWAS" \
  --out_dir "${OUT_DIR}/" \
  --out_file "${OUT_GWAS}.CSFnonsex.txt" \
  --LD_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_pQTL_files/PWAS_LD_Reference_Files/ \
  --LD_file chr1-22.txt


# ==============================
# STEP 2: Munge / clean summary stats (bash)
# ==============================
# LDSC path
LDSC="/storage1/fs1/belloy/Active/01_References/02_Commands_and_Tutorials/fusion_twas-master/ldsc"
cd "${LDSC}"

# Activate conda env (load conda first if needed)
if [ -f "${HOME}/.bashrc" ]; then
   source "${HOME}/.bashrc"
fi
if command -v conda >/dev/null 2>&1; then
  :
elif [ -f /opt/conda/etc/profile.d/conda.sh ]; then
  source /opt/conda/etc/profile.d/conda.sh
else
  echo "ERROR: conda not found. Make sure conda is available on this node." >&2
  exit 1
fi
conda activate ldsc

# ---- GWAS summary statistics ----
munge_sumstats.py --sumstats ${OUT_DIR}/"${OUT_GWAS}.CSFsexstrat.txt" \
  --keep-maf \
  --maf-min 0.001 \
  --a1 ALLELE1 \
  --a2 ALLELE0 \
  --snp SNP \
  --p P \
  --frq A1FREQ \
  --N $NSIZE_GWAS \
  --out ${OUT_DIR}/"${OUT_GWAS}.sexstrat"
  
munge_sumstats.py --sumstats ${OUT_DIR}/"${OUT_GWAS}.CSFnonsex.txt" \
  --keep-maf \
  --maf-min 0.001 \
  --a1 ALLELE1 \
  --a2 ALLELE0 \
  --snp SNP \
  --p P \
  --frq A1FREQ \
  --N $NSIZE_GWAS \
  --out ${OUT_DIR}/"${OUT_GWAS}.nonsex"

# navigate to working directory
cd "${OUT_DIR}"

gzip -d "${OUT_GWAS}.CSFsexstrat.sumstats.gz" "${OUT_GWAS}.CSFnonsex.sumstats.gz"

# ==============================
# STEP 3: Final R cleanup on munged
# ==============================
Rscript /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/Post_mungestat_cleanup.R \
    --file "${OUT_GWAS}.sexstrat.sumstats" \
    --dir "${OUT_DIR}"

Rscript /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/Post_mungestat_cleanup.R \
    --file "${OUT_GWAS}.nonsex.sumstats" \
    --dir "${OUT_DIR}"
