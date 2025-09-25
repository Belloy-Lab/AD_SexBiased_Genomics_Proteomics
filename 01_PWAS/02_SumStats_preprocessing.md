## Sumstats Pre-Processing Overview

This repository documents the pre-processing of sex-stratified GWAS summary statistics for PWAS using brain and CSF proteogenomic datasets.

Pipelines include:
- Genome build harmonization (hg38 → hg19 liftover when required)
- LD reference panel intersection
- Standardization/cleaning with LDSC mungestats for PWAS compatibility

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
bsub -Is -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -a 'docker(dmr07083/fusion-project:4.3.2)' /bin/bash
```
> **Note:** LSF memory syntax can vary by site; `mem=40GB` is accepted on our cluster. If not, convert to MB (e.g., `mem=40960`).

---

## Preprocessing for **Brain** PWAS
The brain proteome weights are in hg19 (GRCh37). GWAS sumstats (initially in hg38) are:

- Lifted over (hg38 → hg19) using UCSC [LiftOver](https://genome-store.ucsc.edu/).
  - Chain file: [hg38ToHg19.over.chain.gz](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz).

- Intersected with ancestry-specific LD reference panels
  - LD reference: 1000Genomes EUR (for EUR cohorts) or EUR/AFR/AMR (for admixed AFR cohorts)

- Cleaned with with LDSC’s `munge_sumstats.py`.
  - Cleaning: LDSC standardizes headers and removes SNPs with missing or invalid fields

### EUR LD reference panel

```bash
  Rscript PWAS/analysis_codes/LDSC_GRCh37_EUR.R --tDir /Path/to/the/1000/Genomes/reference
```

### AFR-admixed LD reference panel

```bash
Rscript PWAS/analysis_codes/LDSC_GRCh37_AFR.R --tDir /Path/to/the/1000/Genomes/reference
```

### Example input formats
#### GWAS sumstats

```bash
head -4 ADGC_ADSP_UKB_FinnGen_Males_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var
posID   SNP     CHR     BP      ALLELE1 ALLELE0 A1FREQ  BETA    SE      P       Q_V     Q_P     N_incl  META_DIR
1:713463:TATTA:T        rs531378844     1       713463  T       TATTA   0.003879        -0.0425036159867374     0.115255916834982       0.712308        4.155881        0.245116        503340  +++-
1:727233:G:A    rs151190501     1       727233  A       G       0.019182        0.0855619103535584      0.0439551962674664      0.051602        2.721039        0.256527        500236  +?++
1:732994:G:A    rs138476838     1       732994  A       G       0.014428        0.00331749103631728     0.0755537825119515      0.964962        3.203519        0.361298        503058  +--+
```

#### LD reference panel

```bash
head -4 1000g_EUR_cm.snplist
SNP     A1      A2
rs367896724     AC      A
rs555500075     TA      T
rs376342519     CCGCCGTTGCAAAGGCGCGCCG  C

head -4 1000g_EUR_cm.tab
CHR     SNP     ALT     REF     BP      posID
1       rs367896724     AC      A       10177   1:10177
1       rs555500075     TA      T       10352   1:10352
1       rs376342519     CCGCCGTTGCAAAGGCGCGCCG  C       10616   1:10616
```

### Running brain pre-cleanup
The input parameters for brain analyses are listed in Brain_PreCleanup_input.csv.

Dry-run first:
```bash
bash PWAS/analysis_codes/Brain_PreCleanup_master.bash PWAS/InputCSV/Brain_PreCleanup_input.csv --dry-run
```

Submit for real:
```bash
bash PWAS/analysis_codes/Brain_PreCleanup_master.bash PWAS/InputCSV/Brain_PreCleanup_input.csv
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
bash PWAS/analysis_codes/CSF_PreCleanup_master.bash PWAS/InputCSV/CSF_PreCleanup_input.csv
```

---

## Notes & Sanity Checks

- **Check `--Nsize_GWAS` in the input CSV**  
  Make sure the sample size column matches the effective N for the specific GWAS subset (male, female, or combined).  
  This value is passed directly to `munge_sumstats.py` and is critical for downstream LDSC/PWAS scaling.  

- **Confirm LD panel selection**  
  Each wrapper script hard-codes LD panel paths. Double-check that:  
  - Brain analyses point to the correct 1000genomes (EUR or AFR-admixed) reference.  
  - CSF analyses point to the correct sex-stratified or combined LD reference files.  

- **Genome build consistency**  
  - Brain workflow: GWAS are lifted from **hg38 → hg19** before LD intersection.  
  - CSF workflow: Everything is **hg38**, so liftover is not applied.  

- **Output naming**  
  Ensure `--out_GWAS` in the CSV uniquely identifies each run. Filenames are constructed using this prefix plus suffixes (e.g., `.CSFsexstrat.txt`, `.hg19_1KG_EURintersected.txt`).  

- **Cluster job allocation**  
  Use `span[hosts=1]` in the LSF `-R` resource string to guarantee single-node allocation. This avoids performance issues from distributed IO.  

---

**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)


