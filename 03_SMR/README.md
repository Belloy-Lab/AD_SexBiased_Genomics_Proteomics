# Summary-based Mendelian Randomization (SMR) Analysis Pipeline

This repository documents a complete pipeline to perform SMR analysis integrating sex-stratified GWAS and pQTL datasets from cerebrospinal fluid (CSF) and brain tissue to identify potential causal genes for Alzheimer's Disease (AD).

---

## 🔬 Overview

SMR analysis is used to test for causal relationships between gene expression (or protein abundance) and complex traits using summary-level data from GWAS and QTL studies. This pipeline focuses on:

* Sex-stratified SMR for female and male AD GWAS datasets
* Integration with protein QTLs (pQTLs) derived from CSF and brain tissues
* Generation of `.besd` format for QTL datasets

Tool Reference: [SMR & HEIDI Software](https://yanglab.westlake.edu.cn/software/smr/#SMR&HEIDIanalysis)

---

## 📁 Input Files

### GWAS Files

* Female GWAS summary: `ADGC_ADSP_UKB_FinnGen_Females_case_control_...`
* Male GWAS summary: `ADGC_ADSP_UKB_FinnGen_Males_case_control_...`

### pQTL Files

* CSF pQTL: separate female, male, and both-sex result files
* Brain pQTL: Wingo et al. pQTL files for female, male, and both

### Gene List

* Candidate gene list with the following required columns:

  * `gene_id`
  * `updated_Gene`
  * `CHR`
  * `Start_hg38`, `End_hg38`
  * `Tissue` (CSF or Brain)

### Reference Genotype

* Binary PLINK reference panel, e.g.:
  `/storage1/fs1/belloy/.../TOPMed_merged`

---

## 🧰 Scripts & Workflow

### Step 1: Clean GWAS Summary Statistics

**Script:** `Clean_GWAS.R`

* Reformats and filters female and male GWAS files
* Outputs cleaned `.ma` files

### Step 2: Clean pQTL Datasets

**Scripts:**

* `Clean_pQTL_CSF.R` — for CSF pQTL
* `Clean_pQTL_Brain.R` — for Brain pQTL
* Filters for ±1MB cis-SNPs and removes SNPs not in GWAS
* Outputs `.esd` and `my.flist` files

### Step 3: Combine flists

**Script:** `Make_besd_file.R`

* Combines CSF and Brain `my.flist` into one
* Runs SMR `--make-besd-dense` to generate `.besd`

### Step 4: Run SMR Analysis

**Executable:** `smr`

Example:

```bash
./smr \
  --bfile /path/to/TOPMed_merged \
  --gwas-summary AD_female.ma \
  --beqtl-summary CSF_Brain_pQTL_besd \
  --peqtl-smr 1e-3 \
  --out smr_res_female \
  --thread-num 10
```

Repeat for male:

```bash
./smr \
  --gwas-summary AD_male.ma \
  --out smr_res_male
```

---

## 📤 Output Files

* `*.ma` — Cleaned GWAS summary files
* `*.esd` — pQTL summary files per gene
* `my.flist` — Metadata for ESD files
* `.besd` — Binary eQTL summary format used for SMR
* `smr_res_female.txt`, `smr_res_male.txt` — Results of SMR

---

## 🖥️ Optional: LSF Job Submission

```bash
bsub -J SMR -n 10 -q subscription \
  -o logs/SMR.%J.out -e logs/SMR.%J.err \
  -R 'rusage[mem=40GB] span[hosts=1]' \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash SMR_submit.sh
```

---

## 📎 Notes

* All scripts are written in R and use `optparse`, `dplyr`, `vroom`, and `data.table`
* The pipeline supports both individual and combined CSF/Brain SMR analysis
* Designed for reproducibility and modular execution

---

## 📬 Contact

For questions or collaboration, please contact the repository owner.

---

**License:** MIT

**Last updated:** July 2025
