# SMR Analysis: Execution Guide

This document outlines the computational steps and corresponding job submission commands to perform Summary-based Mendelian Randomization (SMR) analysis using sex-stratified GWAS and pQTL data from CSF and brain tissue.

## Environment Setup (Required Before Running Any Script)

Before launching jobs or interactive sessions, set the following environment variables:

```bash
export LSF_DOCKER_VOLUMES="/storage1/fs1/belloy/Active:/storage1/fs1/belloy/Active \
/storage2/fs1/belloy2/Active:/storage2/fs1/belloy2/Active /scratch1/fs1/belloy:/scratch1/fs1/belloy $HOME:$HOME"

export CONDA_ENVS_DIRS="/storage1/fs1/belloy/Active/conda/envs/"
export CONDA_PKGS_DIRS="/storage1/fs1/belloy/Active/conda/pkgs/"
export PATH="/opt/conda/bin:$PATH"
export LSF_DOCKER_ENTRYPOINT=/bin/bash
```

### Launch interactive session:

```bash
bsub -Is -G compute-belloy-t1 -q subscription -R 'rusage[mem=40GB]' -a 'docker(dmr07083/fusion-project:4.3.2)' /bin/bash
```

## Project Directory Structure for Analysis
### Set Working Directory
Navigate to your project directory (choose the appropriate path based on your system). Please replace '$USER' with your actual username or preferred folder name. Also, make sure to replace '$USER' with your actual username or preferred folder name in all R and Bash scripts located in the analysis_codes directory.

```bash
cd /storage2/fs1/belloy2/Active/05_Projects/$USER/
mkdir -p SMR/logs
```

### Set Up Code Directory
Navigate to 04_Code directory from storage1 or storage2 and set up project-specific code folders:
```bash
cd /storage2/fs1/belloy2/Active/04_Code/$USER/
mkdir SMR
```
Please download the analysis_codes folder from this repository and copy its entire contents into the 04_Code/$USER/SMR/ directory on the server.
If any folder paths differ from the default setup, make sure to update them accordingly in all R and Bash scripts.

---

## 🧹 Step 1: Clean GWAS Summary Statistics

```bash
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/SMR/Clean_GWAS.R \
  --gwasdir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update \
  --female_gwas ADGC_ADSP_UKB_FinnGen_Females_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var \
  --male_gwas ADGC_ADSP_UKB_FinnGen_Males_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var \
  --outdir /storage2/fs1/belloy2/Active/$USER/SMR \
  --female_ma AD_female.ma \
  --male_ma AD_male.ma
```

---

## 🧼 Step 2: Clean CSF and Brain pQTL Files

### 2A. Clean CSF pQTL files:

```bash
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/SMR/Clean_pQTL_CSF.R \
  --pqtl_female_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/other/pQTL-female-results \
  --pqtl_male_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/other/pQTL-male-results \
  --pqtl_both_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_pQTL_files/Other_Files/02_no_na_no_palindromics \
  --gwasdir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update \
  --female_gwas ADGC_ADSP_UKB_FinnGen_Females_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var \
  --male_gwas ADGC_ADSP_UKB_FinnGen_Males_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var \
  --gene /storage2/fs1/belloy2/Active/05_Projects/$USER/SMR/candidate_gene.txt \
  --output_dir /storage2/fs1/belloy2/Active/$USER/SMR
```

### 2B. Clean Brain pQTL files:

```bash
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/SMR/Clean_pQTL_Brain.R \
  --pqtl_dir /storage1/fs1/belloy/Active/05_Projects/noahc/opt_QTL/liftover/lifted_data \
  --pqtl_female wingo_female_pQTL_appended_lifted_to38.txt \
  --pqtl_male wingo_male_pQTL_appended_lifted_to38.txt \
  --pqtl_both wingo_both_pQTL_appended_lifted_to38.txt \
  --gwasdir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update \
  --female_gwas ADGC_ADSP_UKB_FinnGen_Females_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var \
  --male_gwas ADGC_ADSP_UKB_FinnGen_Males_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var \
  --gene /storage2/fs1/belloy2/Active/05_Projects/$USER/SMR/candidate_gene.txt \
  --output_dir /storage2/fs1/belloy2/Active/$USER/SMR
```

---

## 📦 Step 3: Combine flists & Create `.besd` Format

```bash
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/SMR/Make_besd_file.R \
  --flist_dir /storage2/fs1/belloy2/Active/$USER/SMR \
  --flist_CSF CSF_my.flist \
  --flist_Brain Brain_my.flist \
  --flist_combined CSF_Brain_my.flist \
  --output_dir /storage2/fs1/belloy2/Active/$USER/SMR
```

---

## 🚀 Step 4: Run SMR Analysis

```bash
bsub -g /$USER/compute-belloy -J SMR -n 10 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/SMR/logs/SMR.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/SMR/EU_all/logs/SMR.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' \
  -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/$USER/SMR/SMR_submit.sh
```

---

## 🧶 Step 5: Post-process SMR Results

```bash
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/SMR/Process_SMR_results.R \
  --gene /storage2/fs1/belloy2/Active/05_Projects/$USER/SMR/candidate_gene.txt \
  --smr_res_female smr_res_female.txt \
  --smr_res_male smr_res_male.txt \
  --work_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/SMR \
  --output SMR_res_filtered.txt
```

---

## 📁 Output Summary

* `*.ma` — Cleaned GWAS files
* `*.esd` — QTL summary files per gene
* `my.flist` — Index for ESDs
* `.besd` — Binary QTL format for SMR
* `smr_res_*.txt` — SMR output
* `SMR_res_filtered.txt` — Final filtered results with FDR correction

---
