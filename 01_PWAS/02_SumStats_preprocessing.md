**AD Sex-Biased Genomics & Proteomics**

## Sumstats Pre-Processing Overview

This repository documents the pre-processing of sex-stratified GWAS summary statistics for PWAS using brain and CSF proteogenomic datasets.

Pipelines include:
- Genome build harmonization (hg38 → hg19 liftover when required)
- LD reference panel intersection
- Standardization/cleaning with LDSC mungestats for PWAS compatibility

---

## Preprocessing for **Brain** PWAS
The brain proteome weights are in hg19 (GRCh37). GWAS sumstats (initially in hg38) are:

- Lifted over (hg38 → hg19) using UCSC [LiftOver](https://genome-store.ucsc.edu/).
  - Chain file: [hg38ToHg19.over.chain.gz](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz).

- Intersected with ancestry-specific LD reference panels
  - LD reference: 1000Genomes EUR (for EUR cohorts) or EUR/AFR/AMR (for admixed AFR cohorts)

- Cleaned with with LDSC’s `munge_sumstats.py`.
  - Cleaning: LDSC standardizes headers and removes SNPs with missing or invalid fields

### Create EUR LD reference panel

```bash
  Rscript analysis_codes/LDSC_GRCh37_EUR.R --tDir /Path/to/the/1000/Genomes/reference
```

### Create AFR-admixed LD reference panel

```bash
Rscript analysis_codes/LDSC_GRCh37_AFR.R --tDir /Path/to/the/1000/Genomes/reference
```

### Running brain pre-cleanup
The input parameters for brain analyses are listed in Brain_PreCleanup_input.csv.

Dry-run first:
```bash
bash analysis_codes/Brain_PreCleanup_master.bash InputCSV/Brain_PreCleanup_input.csv --dry-run
```

Submit for real:
```bash
bash analysis_codes/Brain_PreCleanup_master.bash InputCSV/Brain_PreCleanup_input.csv
```
---

## Preprocessing for **CSF** PWAS

The **CSF protein weights** and sex-stratified GWAS sumstats used here are both in **hg38**, so no liftover is required. We intersect sumstats with two LD panels:

1. A **non–sex-stratified** LD reference panel derived from the CSF cohort described in:
     Western, D., et al. (2024). *Proteogenomic analysis of human cerebrospinal fluid identifies neurologically relevant regulation and implicates causal proteins for Alzheimer’s disease*. **Nature Genetics**, 56, 2672–2684. https://doi.org/10.1038/s41588-024-01972-8

2. A **sex-stratified** LD reference panel generated analogously on sex-specific subsets of the cohort.
     > **TODO:** Add exact file paths and generation script references for the sex-stratified CSF LD panel.

After intersection, sumstats are cleaned with `munge_sumstats.py` and filtered to remove missing SNPs.

### Running CSF pre-cleanup
The input parameters for CSF analyses are listed in CSF_PreCleanup_input.csv.

```bash
bash analysis_codes/CSF_PreCleanup_master.bash InputCSV/CSF_PreCleanup_input.csv
```

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)


