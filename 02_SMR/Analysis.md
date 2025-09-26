**AD Sex-Biased Genomics & Proteomics**

## Summary-based Mendelian Randomization (SMR) Analyses

---

## 🧹 Step 1: Clean GWAS Summary Statistics
```bash
Cd 02_SMR

Rscript analysis_codes/Clean_GWAS.R \
  --gwasdir /Path/to/GWAS/files \
  --female_gwas Stage1-2-3_EUR_Female_AD_GWAS \
  --male_gwas Stage1-2-3_EUR_Male_AD_GWAS \
  --outdir /Path/to/output/files \
  --female_ma AD_female.ma \
  --male_ma AD_male.ma
```
---

##  Step 2: Clean CSF and Brain pQTL Files
### 2A. Clean Brain pQTL files:
```bash
Rscript analysis_codes/Clean_pQTL_Brain.R \
  --pqtl_dir /Path/to/brain/pQTL/files \
  --pqtl_female Brain_female_pQTL_lifted19to38.txt \
  --pqtl_male Brain_male_pQTL_lifted19to38.txt \
  --pqtl_both Brain_joint_pQTL_lifted19to38.txt \
  --gwasdir /Path/to/GWAS/files \
  --female_gwas Stage1-2-3_EUR_Female_AD_GWAS \
  --male_gwas Stage1-2-3_EUR_Male_AD_GWAS \
  --gene analyis_codes/candidate_gene.txt \
  --output_dir /Path/to/save/output/files
```

### 2A. Clean CSF pQTL files:
```bash
Rscript analysis_codes/Clean_pQTL_CSF.R \
  --pqtl_female_dir /Female/CSF/pQTL/file/including/full/path \
  --pqtl_male_dir /Male/CSF/pQTL/file/including/full/path \
  --pqtl_both_dir /Joint/CSF/pQTL/file/including/full/path \
  --gwasdir /Path/to/GWAS/files \
  --female_gwas Stage1-2-3_EUR_Female_AD_GWAS \
  --male_gwas Stage1-2-3_EUR_Male_AD_GWAS \
  --gene /storage2/fs1/belloy2/Active/05_Projects/$USER/SMR/candidate_gene.txt \
  --output_dir /Path/to/save/output/files
```
---

## Step 3: Combine flists & Create `.besd` Format
```bash
Rscript analysis_codes/Make_besd_file.R \
  --flist_dir /Path/where/pQTL/output/saved \
  --flist_CSF CSF_pQTL.flist \
  --flist_Brain Brain_pQTL.flist \
  --flist_combined CSF_Brain_merged.flist \
  --output_dir /Path/to/save/output/files
```
---

## Step 4: Run SMR Analysis
```bash
bash analysis_codes/SMR_submit.sh
```
---

## 🧶 Step 5: Post-process SMR Results
```bash
Rscript analysis_codes/Process_SMR_results.R \
  --gene analysis_codes/candidate_gene.txt \
  --smr_res_female smr_res_female.txt \
  --smr_res_male smr_res_male.txt \
  --work_dir work/dir \
  --output SMR_results_filtered.txt
```
---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)

