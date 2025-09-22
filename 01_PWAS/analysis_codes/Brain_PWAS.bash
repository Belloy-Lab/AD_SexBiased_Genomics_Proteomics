#!/usr/bin/env bash
set -euo pipefail

# ========== DROP-IN: sumstats / weights / LD-ref / outdir flags ==========
SS=""; SS_DIR=""; SS_SEX=""
PW_M_DIR=""; PW_M=""
PW_F_DIR=""; PW_F=""
PW_B_DIR=""; PW_B=""
REF_LD_CHR=""
OUT_DIR=""; OUT_FILE=""

usage() { cat <<'USAGE'
PWAS wrapper flags

Inputs (both required):
  --ss_dir        PATH   Directory containing the summary-stats file
  --ss            FILE   Summary-stats filename (no path)

Sex selector (required; accepts M/F/B):
  --ss_sex        STR    male|female|both  (also accepts M|F|B; case-insensitive)

Protein weights (each is a single file paired with its directory):
  --pw_m_dir      PATH   Directory with MALE weights (*.wgt.RDat)
  --pw_m          FILE   Male protein weights file from Brain (inside pw_m_dir)
  --pw_f_dir      PATH   Directory with FEMALE weights
  --pw_f          FILE   Female protein weights file from Brain (inside pw_f_dir)
  --pw_b_dir      PATH   Directory with BOTH/combined weights
  --pw_b          FILE   Both male & female combined protein weights file (inside pw_b_dir)

LD reference (per-chromosome):
  --ref_ld_chr    PATH   Per-chromosome LD ref file with full path.
                         Templates like "...chr{CHR}.ld" / "...chr{chrom}.ld" / "%d" are allowed.

Outputs:
  --out_dir       PATH   Results root directory (created if missing)
  --out_file      STR    Output file name extension (e.g., EUall.cleaned.male)

Aliases also accepted:
  --sumstats, --sumstats_dir, --sumstats_sex, --sex,
  --proteinweights_male_dir, --proteinweights_female_dir, --proteinweights_both_dir
  --proteinweights_male, --proteinweights_female, --proteinweights_both
USAGE
}

# ---- getopt ----
OPTS=$(
  getopt -o '' \
    -l ss:,ss_dir:,sumstats:,sumstats_dir:, \
    -l ss_sex:,sumstats_sex:,sex:, \
    -l pw_m_dir:,pw_m:,proteinweights_male_dir:, \
    -l proteinweights_male:,proteinweights_female:, \
    -l pw_f_dir:,pw_f:,proteinweights_female_dir:, \
    -l pw_b_dir:,pw_b:,proteinweights_both_dir:, \
    -l proteinweights_both:,ref_ld_chr:,out_dir:,out_file:,help \
    -n "$0" -- "$@"
) || { echo "Error: failed to parse options"; exit 2; }
eval set -- "$OPTS"

# ---- parse ----
while true; do
  case "$1" in
    --ss|--sumstats)                           SS="$2"; shift 2 ;;
    --ss_dir|--sumstats_dir)                   SS_DIR="$2"; shift 2 ;;
    --ss_sex|--sumstats_sex|--sex)             SS_SEX="$2"; shift 2 ;;
    --pw_m_dir|--proteinweights_male_dir)      PW_M_DIR="$2"; shift 2 ;;
    --pw_m|--proteinweights_male)              PW_M="$2"; shift 2 ;;
    --pw_f_dir|--proteinweights_female_dir)    PW_F_DIR="$2"; shift 2 ;;
    --pw_f|--proteinweights_female)            PW_F="$2"; shift 2 ;;
    --pw_b_dir|--proteinweights_both_dir)      PW_B_DIR="$2"; shift 2 ;;
    --pw_b|--proteinweights_both)              PW_B="$2"; shift 2 ;;
    --ref_ld_chr)                              REF_LD_CHR="$2"; shift 2 ;;
    --out_dir)                                 OUT_DIR="$2"; shift 2 ;;
    --out_file)                                OUT_FILE="$2"; shift 2 ;;
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

# ---- sumstats: require dir + file and check existence ----
[[ "$SS" == */* ]] && { echo "Error: --ss should be a filename, not a path"; exit 2; }
if [[ -z "$SS_DIR" || -z "$SS" ]]; then
  echo "Error: provide both --ss_dir (directory) and --ss (filename)."; exit 2
fi
SUMSTATS="${SS_DIR%/}/$(basename "$SS")"
[[ -f "$SUMSTATS" ]] || { echo "Error: sumstats not found: $SUMSTATS"; exit 2; }

# ---- required outputs ----
[[ -z "$OUT_DIR"  ]] && { echo "Error: --out_dir is required.";  exit 2; }
[[ -z "$OUT_FILE" ]] && { echo "Error: --out_file is required."; exit 2; }
mkdir -p "$OUT_DIR"

# ---- LD ref: allow concrete file or template ----
if [[ -z "$REF_LD_CHR" ]]; then
  echo "Error: --ref_ld_chr is required (per-chromosome path or template)."; exit 2
fi
if [[ "$REF_LD_CHR" == *"{CHR}"* || "$REF_LD_CHR" == *"{chr}"* || "$REF_LD_CHR" == *"{chrom}"* || "$REF_LD_CHR" == *"%d"* ]]; then
  echo "[INFO] LD-ref template detected: $REF_LD_CHR"
elif [[ ! -e "$REF_LD_CHR" ]]; then
  echo "Warning: --ref_ld_chr '$REF_LD_CHR' not found now; assuming your pipeline expands per chromosome."
fi

# ==============================
# Step 1. Create directory system for PWAS
# ==============================
mkdir -p "${OUT_DIR}/Brain"

sex="${SS_SEX,,}"  # lowercase

if [[ "$sex" == "m" || "$sex" == "male" ]]; then
  mkdir -p "${OUT_DIR}/Brain/Male/Primary" \
           "${OUT_DIR}/Brain/Male/Secondary" \
           "${OUT_DIR}/Brain/Male/Opposite"

elif [[ "$sex" == "f" || "$sex" == "female" ]]; then
  mkdir -p "${OUT_DIR}/Brain/Female/Primary" \
           "${OUT_DIR}/Brain/Female/Secondary" \
           "${OUT_DIR}/Brain/Female/Opposite"

elif [[ "$sex" == "b" || "$sex" == "both" ]]; then
  mkdir -p "${OUT_DIR}/Brain/Both"
fi

# ==============================
# Step 2. PWAS association (3 loops in parallel where applicable)
# ==============================

# quick preflight & helpers
need_file() { [[ -f "$1" ]] || { echo "Error: missing file: $1"; exit 2; }; }
need_dir()  { [[ -d "$1" ]] || { echo "Error: missing directory: $1"; exit 2; }; }

do_assoc() {
  local sumstats="$1" wdir="$2" wfile="$3" ref="$4" chr="$5" outfile="$6"
  # skip if exists (resumable)
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
  # args: label wdir wfile outdir prefix
  local label="$1" wdir="$2" wfile="$3" outdir="$4" prefix="$5"
  need_dir "$wdir"; need_file "${wdir}/${wfile}"
  mkdir -p "$outdir"
  echo "[INFO] Starting loop: $label"
  for chr in {1..22}; do
    do_assoc "$SUMSTATS" "$wdir" "$wfile" "$REF_LD_CHR" "$chr" \
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
    if [[ -s "$f" ]]; then
      files+=("$f")
    else
      echo "[WARN] Missing or empty: $(basename "$f")"
    fi
  done
  (( ${#files[@]} > 0 )) || { echo "[ERR] No per-chr files found for ${prefix}"; return 1; }
  # keep header from first file, drop from subsequent
  awk 'FNR==1 && NR!=1{next} {print}' "${files[@]}" > "$merged"
  echo "[INFO] Merged ${#files[@]} files into $(basename "$merged")"
}

# keep external libs single-threaded to match -n 3 in LSF
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

if [[ "$sex" == "m" || "$sex" == "male" ]]; then
  pids=()
  ( assoc_loop "Male:Primary"   "$PW_M_DIR" "$PW_M" "${OUT_DIR}/Brain/Male/Primary"   "${OUT_FILE}.BrainMalePW"   ) & pids+=($!)
  ( assoc_loop "Male:Secondary" "$PW_B_DIR" "$PW_B" "${OUT_DIR}/Brain/Male/Secondary" "${OUT_FILE}.BrainBothPW"   ) & pids+=($!)
  ( assoc_loop "Male:Opposite"  "$PW_F_DIR" "$PW_F" "${OUT_DIR}/Brain/Male/Opposite"  "${OUT_FILE}.BrainFemalePW" ) & pids+=($!)

  fails=0
  for pid in "${pids[@]}"; do
    wait "$pid" || fails=$((fails+1))
  done
  (( fails == 0 )) || { echo "[ERR] $fails loop(s) failed"; exit 1; }

  # post-merge per arm
  merge_arm "${OUT_DIR}/Brain/Male/Primary"   "${OUT_FILE}.BrainMalePW"   "${OUT_DIR}/Brain/Male/Primary/${OUT_FILE}.BrainMalePW.allchr.dat"
  merge_arm "${OUT_DIR}/Brain/Male/Secondary" "${OUT_FILE}.BrainBothPW"   "${OUT_DIR}/Brain/Male/Secondary/${OUT_FILE}.BrainBothPW.allchr.dat"
  merge_arm "${OUT_DIR}/Brain/Male/Opposite"  "${OUT_FILE}.BrainFemalePW" "${OUT_DIR}/Brain/Male/Opposite/${OUT_FILE}.BrainFemalePW.allchr.dat"

elif [[ "$sex" == "f" || "$sex" == "female" ]]; then
  pids=()
  ( assoc_loop "Female:Primary"   "$PW_F_DIR" "$PW_F" "${OUT_DIR}/Brain/Female/Primary"   "${OUT_FILE}.BrainFemalePW" ) & pids+=($!)
  ( assoc_loop "Female:Secondary" "$PW_B_DIR" "$PW_B" "${OUT_DIR}/Brain/Female/Secondary" "${OUT_FILE}.BrainBothPW"   ) & pids+=($!)
  ( assoc_loop "Female:Opposite"  "$PW_M_DIR" "$PW_M" "${OUT_DIR}/Brain/Female/Opposite"  "${OUT_FILE}.BrainMalePW"   ) & pids+=($!)

  fails=0
  for pid in "${pids[@]}"; do
    wait "$pid" || fails=$((fails+1))
  done
  (( fails == 0 )) || { echo "[ERR] $fails loop(s) failed"; exit 1; }

  # post-merge per arm
  merge_arm "${OUT_DIR}/Brain/Female/Primary"   "${OUT_FILE}.BrainFemalePW" "${OUT_DIR}/Brain/Female/Primary/${OUT_FILE}.BrainFemalePW.allchr.dat"
  merge_arm "${OUT_DIR}/Brain/Female/Secondary" "${OUT_FILE}.BrainBothPW"   "${OUT_DIR}/Brain/Female/Secondary/${OUT_FILE}.BrainBothPW.allchr.dat"
  merge_arm "${OUT_DIR}/Brain/Female/Opposite"  "${OUT_FILE}.BrainMalePW"   "${OUT_DIR}/Brain/Female/Opposite/${OUT_FILE}.BrainMalePW.allchr.dat"

elif [[ "$sex" == "b" || "$sex" == "both" ]]; then
  assoc_loop "Both:Primary" "$PW_B_DIR" "$PW_B" "${OUT_DIR}/Brain/Both" "${OUT_FILE}.BrainBothPW"
  # post-merge
  merge_arm "${OUT_DIR}/Brain/Both" "${OUT_FILE}.BrainBothPW" "${OUT_DIR}/Brain/Both/${OUT_FILE}.BrainBothPW.allchr.dat"
fi
