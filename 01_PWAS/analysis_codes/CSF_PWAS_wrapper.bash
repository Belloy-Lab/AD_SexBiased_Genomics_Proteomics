#!/usr/bin/env bash
set -euo pipefail

# ========== DROP-IN: sumstats / weights / LD-ref / outdir flags ==========
SSX=""; SSN=""; SS_DIR=""; SS_SEX=""
PW_M_DIR=""; PW_M=""
PW_F_DIR=""; PW_F=""
PW_B_DIR=""; PW_B=""
REF_LD_CHR=""
OUT_DIR=""; OUT_FILE=""

usage() { cat <<'USAGE'
PWAS wrapper flags (CSF)

Inputs (all required):
  --ss_dir  PATH   Directory containing summary-stats files
  --ssx     FILE   GWAS sumstats intersected with sex-strat LD panel (filename only)
  --ssn     FILE   GWAS sumstats intersected with non-strat LD panel (filename only)

Sex selector (required; accepts M/F/B):
  --ss_sex  STR    male|female|both  (also M|F|B; case-insensitive)

Protein weights (CSF):
  --pw_m_dir  PATH   Dir with MALE weights (*.wgt.RDat)
  --pw_m      FILE   Male CSF weights file (inside pw_m_dir)
  --pw_f_dir  PATH   Dir with FEMALE weights
  --pw_f      FILE   Female CSF weights file (inside pw_f_dir)
  --pw_b_dir  PATH   Dir with BOTH/combined weights
  --pw_b      FILE   Both CSF weights file (inside pw_b_dir)

LD reference (per-chromosome):
  --ref_ld_chr PATH  Per-chromosome LD ref path (template OK: chr{CHR}, chr{chrom}, or %d)

Outputs:
  --out_dir  PATH   Results root (created if missing)
  --out_file STR    Output prefix/extension (e.g., EUall.cleaned.csf)

Aliases also accepted:
  --sumstats, --sumstats_dir, --sumstats_sex, --sex,
  --proteinweights_male_dir, --proteinweights_female_dir, --proteinweights_both_dir
  --proteinweights_male, --proteinweights_female, --proteinweights_both
USAGE
}

# ---- getopt ----
OPTS=$(
  getopt -o '' \
    -l ssx:,ssn:,ss_dir:,sumstats:,sumstats_dir:, \
    -l ss_sex:,sumstats_sex:,sex:, \
    -l pw_m_dir:,pw_m:,proteinweights_male_dir:,proteinweights_male:, \
    -l pw_f_dir:,pw_f:,proteinweights_female_dir:,proteinweights_female:, \
    -l pw_b_dir:,pw_b:,proteinweights_both_dir:,proteinweights_both:, \
    -l ref_ld_chr:,out_dir:,out_file:,help \
    -n "$0" -- "$@"
) || { echo "Error: failed to parse options"; exit 2; }
eval set -- "$OPTS"

# ---- parse ----
while true; do
  case "$1" in
    --ssx)                                        SSX="$2"; shift 2 ;;
    --ssn)                                        SSN="$2"; shift 2 ;;
    --ss_dir|--sumstats_dir)                      SS_DIR="$2"; shift 2 ;;
    --ss_sex|--sumstats_sex|--sex)                SS_SEX="$2"; shift 2 ;;
    --pw_m_dir|--proteinweights_male_dir)         PW_M_DIR="$2"; shift 2 ;;
    --pw_m|--proteinweights_male)                 PW_M="$2"; shift 2 ;;
    --pw_f_dir|--proteinweights_female_dir)       PW_F_DIR="$2"; shift 2 ;;
    --pw_f|--proteinweights_female)               PW_F="$2"; shift 2 ;;
    --pw_b_dir|--proteinweights_both_dir)         PW_B_DIR="$2"; shift 2 ;;
    --pw_b|--proteinweights_both)                 PW_B="$2"; shift 2 ;;
    --ref_ld_chr)                                 REF_LD_CHR="$2"; shift 2 ;;
    --out_dir)                                    OUT_DIR="$2"; shift 2 ;;
    --out_file)                                   OUT_FILE="$2"; shift 2 ;;
    --help) usage; exit 0 ;;
    --) shift; break ;;
    *) echo "Unknown option: $1"; exit 2 ;;
  esac
done

# ---- normalize & validate sex ----
if [[ -z "$SS_SEX" ]]; then
  echo "Error: --ss_sex is required (male|female|both)."; exit 2
fi
SS_SEX=$(printf '%s' "$SS_SEX" | awk '{print tolower($0)}')
case "$SS_SEX" in
  m|male)   SS_SEX="M" ;;
  f|female) SS_SEX="F" ;;
  b|both)   SS_SEX="B" ;;
  *) echo "Error: --ss_sex must be male|female|both (or M|F|B)."; exit 2 ;;
esac

# ---- sumstats: require dir + both files and check existence ----
[[ "$SSX" == */* ]] && { echo "Error: --ssx should be a filename, not a path"; exit 2; }
[[ "$SSN" == */* ]] && { echo "Error: --ssn should be a filename, not a path"; exit 2; }
if [[ -z "$SS_DIR" || -z "$SSX" || -z "$SSN" ]]; then
  echo "Error: provide --ss_dir, --ssx, and --ssn."; exit 2
fi
SSX_PATH="${SS_DIR%/}/$(basename "$SSX")"
SSN_PATH="${SS_DIR%/}/$(basename "$SSN")"
[[ -f "$SSX_PATH" ]] || { echo "Error: ssx not found: $SSX_PATH"; exit 2; }
[[ -f "$SSN_PATH" ]] || { echo "Error: ssn not found: $SSN_PATH"; exit 2; }

# ---- required outputs ----
[[ -z "$OUT_DIR"  ]] && { echo "Error: --out_dir is required.";  exit 2; }
[[ -z "$OUT_FILE" ]] && { echo "Error: --out_file is required."; exit 2; }
mkdir -p "$OUT_DIR"

# ---- LD ref basic check ----
[[ -z "$REF_LD_CHR" ]] && { echo "Error: --ref_ld_chr is required."; exit 2; }

# ==============================
# Step 1. Create directory system for PWAS (unchanged)
# ==============================
mkdir -p "${OUT_DIR}/CSF"

sex="${SS_SEX,,}"  # lowercase

if [[ "$sex" == "m" || "$sex" == "male" ]]; then
  mkdir -p "${OUT_DIR}/CSF/Male/Primary" \
           "${OUT_DIR}/CSF/Male/Secondary" \
           "${OUT_DIR}/CSF/Male/Opposite"

elif [[ "$sex" == "f" || "$sex" == "female" ]]; then
  mkdir -p "${OUT_DIR}/CSF/Female/Primary" \
           "${OUT_DIR}/CSF/Female/Secondary" \
           "${OUT_DIR}/CSF/Female/Opposite"

elif [[ "$sex" == "b" || "$sex" == "both" ]]; then
  mkdir -p "${OUT_DIR}/CSF/Both"
fi

# ==============================
# Step 2. PWAS association (3 loops in parallel where applicable)
# ==============================

need_file() { [[ -f "$1" ]] || { echo "Error: missing file: $1"; exit 2; }; }
need_dir()  { [[ -d "$1" ]] || { echo "Error: missing directory: $1"; exit 2; }; }

do_assoc() {
  local sumstats="$1" wdir="$2" wfile="$3" ref="$4" chr="$5" outfile="$6"
  if [[ -s "$outfile" ]]; then
    echo "[SKIP] $(basename "$outfile")"
    return 0
  fi
  Rscript /storage1/fs1/belloy/Active/01_References/02_Commands_and_Tutorials/fusion_twas-master/FUSION.assoc_test.R \
    --sumstats    "$sumstats" \
    --weights     "${wdir}/${wfile}" \
    --weights_dir "${wdir}" \
    --ref_ld_chr  "$ref" \
    --chr         "$chr" \
    --out         "$outfile"
}

assoc_loop() {
  # args: label sumstats wdir wfile outdir prefix
  local label="$1" sumstats="$2" wdir="$3" wfile="$4" outdir="$5" prefix="$6"
  need_dir "$wdir"; need_file "${wdir}/${wfile}"
  mkdir -p "$outdir"
  echo "[INFO] Starting loop: $label"
  for chr in {1..22}; do
    do_assoc "$sumstats" "$wdir" "$wfile" "$REF_LD_CHR" "$chr" \
             "${outdir}/${prefix}.${chr}.dat"
  done
  echo "[INFO] Finished loop: $label"
}

merge_arm() {
  # args: outdir prefix merged_out
  local outdir="$1" prefix="$2" merged="$3"
  local files=() chr f
  for chr in {1..22}; do
    f="${outdir}/${prefix}.${chr}.dat"
    [[ -s "$f" ]] && files+=("$f") || echo "[WARN] Missing or empty: $(basename "$f")"
  done
  (( ${#files[@]} > 0 )) || { echo "[ERR] No per-chr files for ${prefix}"; return 1; }
  awk 'FNR==1 && NR!=1{next} {print}' "${files[@]}" > "$merged"
  echo "[INFO] Merged ${#files[@]} files into $(basename "$merged")"
}

# keep external libs single-threaded to match -n 3 in LSF
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

if [[ "$sex" == "m" || "$sex" == "male" ]]; then
  pids=()
  # Primary: SSX vs Male weights
  ( assoc_loop "CSF Male:Primary"   "$SSX_PATH" "$PW_M_DIR" "$PW_M" "${OUT_DIR}/CSF/Male/Primary"   "${OUT_FILE}.CSFMalePW"   ) & pids+=($!)
  # Secondary: SSN vs Both weights
  ( assoc_loop "CSF Male:Secondary" "$SSN_PATH" "$PW_B_DIR" "$PW_B" "${OUT_DIR}/CSF/Male/Secondary" "${OUT_FILE}.CSFBothPW"   ) & pids+=($!)
  # Opposite: SSX vs Female weights
  ( assoc_loop "CSF Male:Opposite"  "$SSX_PATH" "$PW_F_DIR" "$PW_F" "${OUT_DIR}/CSF/Male/Opposite"  "${OUT_FILE}.CSFFemalePW" ) & pids+=($!)

  fails=0; for pid in "${pids[@]}"; do wait "$pid" || fails=$((fails+1)); done
  (( fails == 0 )) || { echo "[ERR] $fails loop(s) failed"; exit 1; }

  merge_arm "${OUT_DIR}/CSF/Male/Primary"   "${OUT_FILE}.CSFMalePW"   "${OUT_DIR}/CSF/Male/Primary/${OUT_FILE}.CSFMalePW.allchr.dat"
  merge_arm "${OUT_DIR}/CSF/Male/Secondary" "${OUT_FILE}.CSFBothPW"   "${OUT_DIR}/CSF/Male/Secondary/${OUT_FILE}.CSFBothPW.allchr.dat"
  merge_arm "${OUT_DIR}/CSF/Male/Opposite"  "${OUT_FILE}.CSFFemalePW" "${OUT_DIR}/CSF/Male/Opposite/${OUT_FILE}.CSFFemalePW.allchr.dat"

elif [[ "$sex" == "f" || "$sex" == "female" ]]; then
  pids=()
  # Primary: SSX vs Female weights
  ( assoc_loop "CSF Female:Primary"   "$SSX_PATH" "$PW_F_DIR" "$PW_F" "${OUT_DIR}/CSF/Female/Primary"   "${OUT_FILE}.CSFFemalePW" ) & pids+=($!)
  # Secondary: SSN vs Both weights
  ( assoc_loop "CSF Female:Secondary" "$SSN_PATH" "$PW_B_DIR" "$PW_B" "${OUT_DIR}/CSF/Female/Secondary" "${OUT_FILE}.CSFBothPW"   ) & pids+=($!)
  # Opposite: SSX vs Male weights
  ( assoc_loop "CSF Female:Opposite"  "$SSX_PATH" "$PW_M_DIR" "$PW_M" "${OUT_DIR}/CSF/Female/Opposite"  "${OUT_FILE}.CSFMalePW"   ) & pids+=($!)

  fails=0; for pid in "${pids[@]}"; do wait "$pid" || fails=$((fails+1)); done
  (( fails == 0 )) || { echo "[ERR] $fails loop(s) failed"; exit 1; }

  merge_arm "${OUT_DIR}/CSF/Female/Primary"   "${OUT_FILE}.CSFFemalePW" "${OUT_DIR}/CSF/Female/Primary/${OUT_FILE}.CSFFemalePW.allchr.dat"
  merge_arm "${OUT_DIR}/CSF/Female/Secondary" "${OUT_FILE}.CSFBothPW"   "${OUT_DIR}/CSF/Female/Secondary/${OUT_FILE}.CSFBothPW.allchr.dat"
  merge_arm "${OUT_DIR}/CSF/Female/Opposite"  "${OUT_FILE}.CSFMalePW"   "${OUT_DIR}/CSF/Female/Opposite/${OUT_FILE}.CSFMalePW.allchr.dat"

elif [[ "$sex" == "b" || "$sex" == "both" ]]; then
  # Both analysis: SSN vs Both weights
  assoc_loop "CSF Both:Primary" "$SSN_PATH" "$PW_B_DIR" "$PW_B" "${OUT_DIR}/CSF/Both" "${OUT_FILE}.CSFBothPW"
  merge_arm  "${OUT_DIR}/CSF/Both" "${OUT_FILE}.CSFBothPW" "${OUT_DIR}/CSF/Both/${OUT_FILE}.CSFBothPW.allchr.dat"
fi
