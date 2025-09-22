**AD Sex-Biased Genomics & Proteomics**

## Sumstats Pre-Processing Overview

This repository documents the pre-processing of **sex-stratified GWAS summary statistics** for PWAS using **brain** and **CSF** proteogenomics datasets. Pipelines include genome build harmonization (liftover), LD-panel intersection, and LDSC-compatible cleaning.

---

## Environment Setup (Required Before Running Any Script)

Set the following environment variables before launching jobs or interactive sessions on LSF:

```bash
export LSF_DOCKER_VOLUMES="/storage1/fs1/belloy/Active:/storage1/fs1/belloy/Active \
/storage2/fs1/belloy2/Active:/storage2/fs1/belloy2/Active /scratch1/fs1/belloy:/scratch1/fs1/belloy $HOME:$HOME"

export CONDA_ENVS_DIRS="/storage1/fs1/belloy/Active/conda/envs/"
export CONDA_PKGS_DIRS="/storage1/fs1/belloy/Active/conda/pkgs/"
export PATH="/opt/conda/bin:$PATH"
export LSF_DOCKER_ENTRYPOINT=/bin/bash
```

### Optional: Launch an interactive session

Use the fusion-project Docker image to test commands interactively (optional).

```bash
bsub -Is -G compute-belloy-t1 -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -a 'docker(dmr07083/fusion-project:4.3.2)' /bin/bash
```
> **Note:** LSF memory syntax can vary by site; `mem=40GB` is accepted on our cluster. If not, convert to MB (e.g., `mem=40960`).

---

## Preprocessing for **Brain** PWAS

The brain proteome weights are in **hg19 (GRCh37)**. Sex-stratified GWAS sumstats are first **lifted from hg38 → hg19** using UCSC LiftOver, then **intersected with a 1KGP LD reference panel**, and finally cleaned with LDSC’s `munge_sumstats.py`.

- Chain file: [hg38ToHg19.over.chain.gz](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz)
- Intersection: 1KGP **EUR** panel for EU cohorts; 1KGP **AFR-admixed** panel for AFR cohorts.
- Cleaning: `munge_sumstats.py` (LDSC) to standardize/filter; remove SNPs with missing values.
  - LDSC repo: https://github.com/bulik/ldsc

### EUR LD reference panel (build GRCh37)

```bash
bsub -g /$USER/compute-belloy -J EURLD -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/EUR1KGP_GR37.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/EUR1KGP_GR37.%J.err \
  -q subscription -R 'rusage[mem=32GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/EU_all/create_LDSC_GRCh37_EUR_p1.R
```

### AFR-admixed LD reference panel (build GRCh37)

```bash
bsub -g /$USER/compute-belloy -J AFRLDpanel -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/AFR-ad_LD_ref_panel.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/AFR-ad_LD_ref_panel.%J.err \
  -q subscription -R 'rusage[mem=30GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/AFR/create_LDSC_GRCh37_AFR-admixed.R \
      --PATH_plink1.9 /usr/bin/plink1.9 \
      --PATH_plink2 /usr/bin/plink2
```

### Brain: EU_all

```bash
# Male
bsub -g /$USER/compute-belloy -J Br_PreClean_Male_EUall -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/Br_PreClean_male_EUall.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/Br_PreClean_male_EUall.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/Brain_PreCleanup_wrapper.bash \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update \
    --in_GWAS ADGC_ADSP_UKB_FinnGen_Males_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var \
    --out_GWAS EUall.cleaned.male.brain \
    --Nsize_GWAS 507675 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/PWAS/EU_all \
    --pop EUR

# Female
bsub -g /$USER/compute-belloy -J Br_PreClean_Female_EUall -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/Br_PreClean_female_EUall.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/Br_PreClean_female_EUall.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/Brain_PreCleanup_wrapper.bash \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update \
    --in_GWAS ADGC_ADSP_UKB_FinnGen_Females_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var \
    --out_GWAS EUall.cleaned.female.brain \
    --Nsize_GWAS 636154 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/PWAS/EU_all \
    --pop EUR
```

### Brain: EU_noUKB

```bash
# Male
bsub -g /$USER/compute-belloy -J Br_PreClean_Male_noUKB -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/Br_PreClean_male_noUKB.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/Br_PreClean_male_noUKB.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/Brain_PreCleanup_wrapper.bash \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update \
    --in_GWAS ADGC_ADSP_FinnGen_Males_case_control.gwama.clean.gen090.exclude_APOE_region.shared_var \
    --out_GWAS EUnoUKB.cleaned.male.brain \
    --Nsize_GWAS 230632 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/PWAS/EU_noUKB \
    --pop EUR

# Female
bsub -g /$USER/compute-belloy -J Br_PreClean_Female_noUKB -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/Br_PreClean_female_noUKB.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/Br_PreClean_female_noUKB.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/Brain_PreCleanup_wrapper.bash \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update \
    --in_GWAS ADGC_ADSP_FinnGen_Females_case_control.gwama.clean.gen090.exclude_APOE_region.shared_var \
    --out_GWAS EUnoUKB.cleaned.female.brain \
    --Nsize_GWAS 300237 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/PWAS/EU_noUKB \
    --pop EUR
```

### Brain: AFR (AFR-admixed)

```bash
# Male
bsub -g /$USER/compute-belloy -J Br_PreClean_Male_AFR -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/Br_PreClean_male_AFR.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/Br_PreClean_male_AFR.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/Brain_PreCleanup_wrapper.bash \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_sex_apoe4_GWAS_AD_AFR \
    --in_GWAS AFR_admixed_Male_case_control.gwama.clean.full.txt \
    --out_GWAS AFR.cleaned.male.brain \
    --Nsize_GWAS 2116 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/PWAS/AFR \
    --pop AFR

# Female
bsub -g /$USER/compute-belloy -J Br_PreClean_Female_AFR -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/Br_PreClean_female_AFR.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/Br_PreClean_female_AFR.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/Brain_PreCleanup_wrapper.bash \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_sex_apoe4_GWAS_AD_AFR \
    --in_GWAS AFR_admixed_Females_case_control.gwama.clean.full.txt \
    --out_GWAS AFR.cleaned.female.brain \
    --Nsize_GWAS 5149 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/PWAS/AFR \
    --pop AFR
```

---

## Preprocessing for **CSF** PWAS

The **CSF protein weights** and sex-stratified GWAS sumstats used here are both in **hg38**, so no liftover is required. We intersect sumstats with two LD panels:

1. A **non–sex-stratified** LD reference panel derived from the CSF cohort described in:
   
   Western, D., et al. (2024). *Proteogenomic analysis of human cerebrospinal fluid identifies neurologically relevant regulation and implicates causal proteins for Alzheimer’s disease*. **Nature Genetics**, 56, 2672–2684. https://doi.org/10.1038/s41588-024-01972-8

2. A **sex-stratified** LD reference panel generated analogously on sex-specific subsets of the cohort.
   
   > **TODO:** Add exact file paths and generation script references for the sex-stratified CSF LD panel.

After intersection, sumstats are cleaned with `munge_sumstats.py`, and SNPs with missing values are removed.

### CSF: EU_all

```bash
# Male
bsub -g /$USER/compute-belloy -J CSF_PreClean_Male_EUall -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/CSF_PreClean_male_EUall.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/CSF_PreClean_male_EUall.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/CSF_PreCleanup_wrapper.bash \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update \
    --in_GWAS ADGC_ADSP_UKB_FinnGen_Males_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var \
    --out_GWAS EUall.cleaned.male.CSF \
    --Nsize_GWAS 507675 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/PWAS/EU_all

# Female
bsub -g /$USER/compute-belloy -J CSF_PreClean_Female_EUall -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/CSF_PreClean_female_EUall.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/CSF_PreClean_female_EUall.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/CSF_PreCleanup_wrapper.bash \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update \
    --in_GWAS ADGC_ADSP_UKB_FinnGen_Females_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var \
    --out_GWAS EUall.cleaned.female.CSF \
    --Nsize_GWAS 636154 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/PWAS/EU_all
```

### CSF: EU_noUKB

```bash
# Male
bsub -g /$USER/compute-belloy -J CSF_PreClean_Male_noUKB -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/CSF_PreClean_male_noUKB.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/CSF_PreClean_male_noUKB.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/CSF_PreCleanup_wrapper.bash \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update \
    --in_GWAS ADGC_ADSP_FinnGen_Males_case_control.gwama.clean.gen090.exclude_APOE_region.shared_var \
    --out_GWAS EUnoUKB.cleaned.male.CSF \
    --Nsize_GWAS 230632 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/PWAS/EU_noUKB

# Female
bsub -g /$USER/compute-belloy -J CSF_PreClean_Female_noUKB -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/CSF_PreClean_female_noUKB.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/CSF_PreClean_female_noUKB.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/CSF_PreCleanup_wrapper.bash \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update \
    --in_GWAS ADGC_ADSP_FinnGen_Females_case_control.gwama.clean.gen090.exclude_APOE_region.shared_var \
    --out_GWAS EUnoUKB.cleaned.female.CSF \
    --Nsize_GWAS 300237 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/PWAS/EU_noUKB
```

### CSF: AFR (AFR-admixed)

```bash
# Male
bsub -g /$USER/compute-belloy -J CSF_PreClean_Male_AFR -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/CSF_PreClean_male_AFR.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/CSF_PreClean_male_AFR.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/CSF_PreCleanup_wrapper.bash \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_sex_apoe4_GWAS_AD_AFR \
    --in_GWAS AFR_admixed_Male_case_control.gwama.clean.full.txt \
    --out_GWAS AFR.cleaned.male.CSF \
    --Nsize_GWAS 2116 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/PWAS/AFR

# Female
bsub -g /$USER/compute-belloy -J CSF_PreClean_Female_AFR -n 1 -N \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/CSF_PreClean_female_AFR.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/CSF_PreClean_female_AFR.%J.err \
  -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  bash /storage2/fs1/belloy2/Active/04_Code/sivas/PWAS/CSF_PreCleanup_wrapper.bash \
    --in_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_sex_apoe4_GWAS_AD_AFR \
    --in_GWAS AFR_admixed_Females_case_control.gwama.clean.full.txt \
    --out_GWAS AFR.cleaned.female.CSF \
    --Nsize_GWAS 5149 \
    --out_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/PWAS/AFR
```

---

## Notes & Sanity Checks

- Ensure **`--Nsize_GWAS`** matches the effective sample size for the given sex-stratified subset.
- `--pop` controls which LD panel is used inside the wrappers (EUR vs AFR-admixed). Confirm the panel paths in the wrapper match the examples above.
- For **brain** runs, confirm liftover hg38→hg19 is applied **before** intersecting with hg19 LD panels.
- For **CSF** runs (hg38 throughout), confirm all inputs/panels are hg38 and **skip liftover**.
- Use `span[hosts=1]` in `-R` to keep a single-node allocation (useful for file IO locality).

---

## License (MIT)

Copyright (c) 2025 Sathesh K. Sivasankaran

