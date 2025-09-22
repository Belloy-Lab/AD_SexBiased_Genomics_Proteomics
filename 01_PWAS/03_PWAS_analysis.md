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
bsub -Is -G compute-belloy-t1 -q subscription -R 'rusage[mem=40GB] span[hosts=1]' -a 'docker(dmr07083/fusion-project:4.3.2)' /bin/bash
```

## PWAS Analysis Overview

We performed sex-stratified protein-wide association studies (PWAS) using FUSION (http://gusevlab.org/projects/fusion/). Alzheimer’s disease (AD) GWAS (male, and female) were paired with proteogenomic prediction models (male. female and combined) from brain and CSF to test whether genetically predicted protein abundance is associated with AD in a sex-specific manner—supporting a causal role.

- Primary (sex-matched) analysis
    - Male GWAS × male protein weights, and female GWAS × female protein weights (brain or CSF). These constitute the main PWAS discoveries.

- Secondary (non-sex-matched) analysis
    - Male or female GWAS × combined (male+female) protein weights. Used to assess robustness when sex-specific proteomic models are not available.

- Opposite-sex (cross-sex) analysis
    - Male GWAS × female protein weights, and female GWAS × male protein weights. Used as a sensitivity check for sex effects.

PWAS bash wrappers create all necessary directories, manage logs, and organize outputs by tissue, sex, and ancestry panel for reproducibility.

### Brain Proteogenomics analysis
#### EU_all
```bash
# Male
bsub -g /$USER/compute-belloy -q subscription -J PwasMaleBrain -n 3 \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/PWASMAleBrain.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/PWASMAleBrain.%J.err \
  -R 'span[hosts=1] rusage[mem=100000]' -M 100000 -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/Brain_pwas.bash \
        --ss_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all \
        --ss EUall_Male.hg19_intersected.sumstats \
        --ss_sex male \
        --pw_m_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/Males \
        --pw_m train_weights.pos \
        --pw_f_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/Females \
        --pw_f train_weights.pos \
        --pw_b_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/joint_weights \
        --pw_b train_weights.pos \
        --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/1000G_Grch37/1KG_EUR/1000g_EUR_cm_ch \
        --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all \
        --out_file EUall.cleaned.brain.male
        
# female
bsub -g /$USER/compute-belloy -q subscription -J PwasFemaleBrain -n 3 \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/PwasFemaleBrain.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/PwasFemaleBrain.%J.err \
  -R 'span[hosts=1] rusage[mem=100000]' -M 100000 -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/Brain_pwas.bash \
        --ss_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all \
        --ss EUall_Female.hg19_intersected.sumstats \
        --ss_sex mafemalele \
        --pw_m_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/Males \
        --pw_m train_weights.pos \
        --pw_f_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/Females \
        --pw_f train_weights.pos \
        --pw_b_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/joint_weights \
        --pw_b train_weights.pos \
        --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/1000G_Grch37/1KG_EUR/1000g_EUR_cm_ch \
        --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all \
        --out_file EUall.cleaned.brain.female
```

#### EU_noUKB
```bash
# Male
bsub -g /$USER/compute-belloy -q subscription -J PwasMaleBrain -n 3 \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/PwasMaleBrain.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/PwasMaleBrain.%J.err \
  -R 'span[hosts=1] rusage[mem=100000]' -M 100000 -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/Brain_pwas.bash \
        --ss_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB \
        --ss EUnoUKB_Male.hg19_intersected.sumstats \
        --ss_sex male \
        --pw_m_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/Males \
        --pw_m train_weights.pos \
        --pw_f_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/Females \
        --pw_f train_weights.pos \
        --pw_b_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/joint_weights \
        --pw_b train_weights.pos \
        --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/1000G_Grch37/1KG_EUR/1000g_EUR_cm_ch \
        --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB \
        --out_file EUnoUKB.cleaned.brain.male
        
# female
bsub -g /$USER/compute-belloy -q subscription -J PwasFemaleBrain -n 3 \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/PwasFemaleBrain.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/PwasFemaleBrain.%J.err \
  -R 'span[hosts=1] rusage[mem=100000]' -M 100000 -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/Brain_pwas.bash \
        --ss_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB \
        --ss EUnoUKB_Female.hg19_intersected.sumstats \
        --ss_sex female \
        --pw_m_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/Males \
        --pw_m train_weights.pos \
        --pw_f_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/Females \
        --pw_f train_weights.pos \
        --pw_b_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/joint_weights \
        --pw_b train_weights.pos \
        --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/1000G_Grch37/1KG_EUR/1000g_EUR_cm_ch \
        --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB \
        --out_file EUnoUKB.cleaned.brain.female
```

#### AFR
```bash
# Male
bsub -g /$USER/compute-belloy -q subscription -J PwasMaleBrain -n 3 \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/PwasMaleBrain.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/PwasMaleBrain.%J.err \
  -R 'span[hosts=1] rusage[mem=100000]' -M 100000 -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/Brain_pwas.bash \
        --ss_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR \
        --ss AFR_Male.hg19_intersected.sumstats \
        --ss_sex male \
        --pw_m_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/Males \
        --pw_m train_weights.pos \
        --pw_f_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/Females \
        --pw_f train_weights.pos \
        --pw_b_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/joint_weights \
        --pw_b train_weights.pos \
        --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/1000G_Grch37/1KG_EUR/1000g_EUR_cm_ch \
        --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR \
        --out_file AFR.cleaned.brain.male
        
# female
bsub -g /$USER/compute-belloy -q subscription -J PwasFemaleBrain -n 3 \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/PwasFemaleBrain.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/PwasFemaleBrain.%J.err \
  -R 'span[hosts=1] rusage[mem=100000]' -M 100000 -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/Brain_pwas.bash \
        --ss_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR \
        --ss AFR_Female.hg19_intersected.sumstats \
        --ss_sex female \
        --pw_m_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/Males \
        --pw_m train_weights.pos \
        --pw_f_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/Females \
        --pw_f train_weights.pos \
        --pw_b_dir /storage1/fs1/belloy/Active/02_Data/03_Freeze/Wingoetal2023_protein_weights/joint_weights \
        --pw_b train_weights.pos \
        --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/1000G_Grch37/1KG_EUR/1000g_EUR_cm_ch \
        --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR \
        --out_file AFR.cleaned.brain.female
```

### CSF Proteogenomics analysis
#### EU_all
```bash
# Male
bsub -g /$USER/compute-belloy -q subscription -J PwasMaleCSF -n 3 \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/PwasMaleCSF.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/PwasMaleCSF.%J.err \
  -R 'span[hosts=1] rusage[mem=150000]' -M 150000 -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/CSF_pwas.bash \
        --ss_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all \
        --ss EUall_Male.hg19_intersected.sumstats \
        --ss_sex male \
        --pw_m_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_m NGI_CSF_male_cis_weights.pos \
        --pw_f_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_f NGI_CSF_female_cis_weights.pos \
        --pw_b_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_b NGI_CSF_cis-weights.pos \
        --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/LD/LD/CSF_proteomics_3107_samples_maf_0032_geno_1_chr_ \
        --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all \
        --out_file EUall.cleaned.CSF.male
        
# female
bsub -g /$USER/compute-belloy -q subscription -J PwasFemaleBrain -n 3 \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/PwasFemaleBrain.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/PwasFemaleBrain.%J.err \
  -R 'span[hosts=1] rusage[mem=150000]' -M 150000 -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/Brain_pwas.bash \
        --ss_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all \
        --ss EUall_Female.hg19_intersected.sumstats \
        --ss_sex female \
        --pw_m_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_m NGI_CSF_male_cis_weights.pos \
        --pw_f_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_f NGI_CSF_female_cis_weights.pos \
        --pw_b_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_b NGI_CSF_cis-weights.pos \
        --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/LD/LD/CSF_proteomics_3107_samples_maf_0032_geno_1_chr_ \
        --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all \
        --out_file EUall.cleaned.CSF.female

# HP
bsub -g /$USER$/compute-belloy -J EUCSFHPall -n 2 -N \
-o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/EUall_CSF_HP_PWAS.%J.out \
-e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/logs/EUall_CSF_HP_PWAS.%J.err \
-q subscription -R 'rusage[mem=100GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 -a 'docker(dmr07083/fusion-project:4.3.2)' \
bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/EUall_CSF_HP_PWAS.bash
```

#### EU_noUKB
```bash
# Male
bsub -g /$USER/compute-belloy -q subscription -J PwasMaleCSF -n 3 \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/PwasMaleCSF.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/PwasMaleCSF.%J.err \
  -R 'span[hosts=1] rusage[mem=150000]' -M 150000 -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/CSF_pwas.bash \
        --ss_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB \
        --ss EUnoUKB_Male.hg19_intersected.sumstats \
        --ss_sex male \
        --pw_m_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_m NGI_CSF_male_cis_weights.pos \
        --pw_f_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_f NGI_CSF_female_cis_weights.pos \
        --pw_b_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_b NGI_CSF_cis-weights.pos \
        --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/LD/LD/CSF_proteomics_3107_samples_maf_0032_geno_1_chr_ \
        --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB \
        --out_file EUnoUKB.cleaned.CSF.male
        
# female
bsub -g /$USER/compute-belloy -q subscription -J PwasFemaleCSF -n 3 \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/PwasFemaleCSF.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/PwasFemaleCSF.%J.err \
  -R 'span[hosts=1] rusage[mem=150000]' -M 150000 -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/CSF_pwas.bash \
        --ss_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB \
        --ss EUnoUKB_Female.hg19_intersected.sumstats \
        --ss_sex female \
        --pw_m_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_m NGI_CSF_male_cis_weights.pos \
        --pw_f_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_f NGI_CSF_female_cis_weights.pos \
        --pw_b_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_b NGI_CSF_cis-weights.pos \
        --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/LD/LD/CSF_proteomics_3107_samples_maf_0032_geno_1_chr_ \
        --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB \
        --out_file EUnoUKB.cleaned.CSF.female

# HP
bsub -g /$USER$/compute-belloy -J EUnoUKNCSFHP -n 1 -N \
-o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/EUnoUKNCSFHP.%J.out \
-e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/logs/EUnoUKNCSFHP.%J.err \
-q subscription -R 'rusage[mem=100GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 -a 'docker(dmr07083/fusion-project:4.3.2)' \
bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/EUnoUKB_CSF_HP_PWAS.bash
```

#### AFR
```bash
# Male
bsub -g /$USER/compute-belloy -q subscription -J PwasMaleCSF -n 3 \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/PwasMaleCSF.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/PwasMaleCSF.%J.err \
  -R 'span[hosts=1] rusage[mem=100000]' -M 100000 -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/CSF_pwas.bash \
        --ss_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR \
        --ss AFR_Male.hg19_intersected.sumstats \
        --ss_sex male \
        --pw_m_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_m NGI_CSF_male_cis_weights.pos \
        --pw_f_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_f NGI_CSF_female_cis_weights.pos \
        --pw_b_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_b NGI_CSF_cis-weights.pos \
        --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/LD/LD/CSF_proteomics_3107_samples_maf_0032_geno_1_chr_ \
        --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR \
        --out_file AFR.cleaned.CSF.male
        
# female
bsub -g /$USER/compute-belloy -q subscription -J PwasFemaleCSF -n 3 \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/PwasFemaleCSF.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/PwasFemaleCSF.%J.err \
  -R 'span[hosts=1] rusage[mem=100000]' -M 100000 -G compute-belloy-t1 -sla compute-belloy-t1 \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
    bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/CSF_pwas.bash \
        --ss_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR \
        --ss AFR_Female.hg19_intersected.sumstats \
        --ss_sex female \
        --pw_m_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_m NGI_CSF_male_cis_weights.pos \
        --pw_f_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_f NGI_CSF_female_cis_weights.pos \
        --pw_b_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
        --pw_b NGI_CSF_cis-weights.pos \
        --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/LD/LD/CSF_proteomics_3107_samples_maf_0032_geno_1_chr_ \
        --out_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR \
        --out_file AFR.cleaned.CSF.female

# HP
bsub -g /$USER/compute-belloy -J AFRCSFHP -n 1 -N \
-o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/AFRCSFHP.%J.out \
-e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR/logs/AFRCSFHP.%J.err \
-q subscription -R 'rusage[mem=60GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 -a 'docker(dmr07083/fusion-project:4.3.2)' \
bash /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/AFR_CSF_HP_PWAS.bash
```

