**AD Sex-Biased Genomics & Proteomics**

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
bsub -Is -G compute-belloy-t1 -q subscription -R 'rusage[mem=150GB] span[hosts=1]' -a 'docker(dmr07083/fusion-project:4.3.2)' /bin/bash
```

## PWAS Analysis Overview

We performed sex-stratified protein-wide association studies (PWAS) using FUSION (http://gusevlab.org/projects/fusion/). AD GWAS (male and female) were paired with proteogenomic prediction models (male, female, and combined) from brain and CSF.

- Primary (sex-matched) analysis
    - Male GWAS × male protein weights; female GWAS × female protein weights.

- Secondary (non-sex-matched) analysis
    - Male GWAS × combined (male+female) protein weights; female GWAS × combined (male+female) protein weights. 

- Opposite-sex (cross-sex) analysis
    - Male GWAS × female protein weights; female GWAS × male protein weights. Used as only to create sex specificity filters.


### Input file formats

#### GWAS summary statistics
Must include alleles, effect sizes (Z or beta), frequency, and sample size.

```bash
head -4 EUall.cleaned.female.brain.hg19_intersected.sumstats
A1      A2      Z       SNP     N       FRQ
T       G       -1.232  rs10875231      636154.000      0.264
A       G       -1.244  rs186077422     636154.000      0.006
T       C       -1.036  rs6678176       636154.000      0.329
```

#### Protein weights
Provided in pos files listing protein weight sets, IDs, and genomic ranges.

```bash
head -4 NGI_CSF_female_cis_weights.pos
WGT     ID      CHR     P0      P1
WGT/X10037.98__female_weight.wgt.RDat   X10037.98       19      51503097        51503097
WGT/X10361.25__female_weight.wgt.RDat   X10361.25       12      112943944       112943944
WGT/X10396.6__female_weight.wgt.RDat    X10396.6        1       150578851       150578851
```

#### LD reference panel
Standard PLINK format (.bed/.bim/.fam) separated by chromosome. Used to calculate correlations between SNPs.

### Brain Proteogenomics analysis
The input parameters for brain analyses are listed in Brain_PWAS_input.csv.

Dry-run first:
```bash
bash PWAS/analysis_codes/Brain_PWAS_master.bash PWAS/Input_CSV/Brain_PWAS_input.csv --dry-run
```

Submit for real:
```bash
bash PWAS/analysis_codes/Brain_PWAS_master.bash PWAS/Input_CSV/Brain_PWAS_input.csv
```

### CSF Proteogenomics analysis
The input parameters for CSF analyses are listed in CSF_PWAS_input.csv.

```bash
bash PWAS/analysis_codes/CSF_PWAS_master.bash PWAS/Input_CSV/CSF_PWAS_input.csv
```

### CSF HP weights analysis
The same input parameters file CSF_PWAS_input.csv is used to HP analysis.

```bash
bash PWAS/analysis_codes/CSF_PWAS_master_HP.bash PWAS/Input_CSV/CSF_PWAS_input.csv
```


### Typical Output Files
Each analysis generates the following key outputs (naming depends on --out_GWAS and input CSV):

- Log files:
  - logs/PwasMaleBrain.<jobID>.out and logs/PwasMaleBrain.<jobID>.err

- PWAS association results:
  - <out_dir>/<tissue>/<sex>/*.chr.dat

**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)

