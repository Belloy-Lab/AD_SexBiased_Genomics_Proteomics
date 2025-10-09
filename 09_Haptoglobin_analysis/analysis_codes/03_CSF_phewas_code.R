## Haptoglobin CSF pheWAS Code

library(data.table)
library(tidyverse)


#################################################################################
### Define consistent input parameters

# Plink2 command
pl <- "plink2"
# Plink file
in_file <- "Input_data/CSF_proteomics_LONGS_PPMI_Stanford_GARFIELD_samples_geno_0.1_mac_5" # PLINK format
# HP region reference file
reg_file <- "Input_data/HP_region_variants_1MB_radius.txt"
# Covariate file
cov_file <- "Input_data/HP1_genotype_analysis_covariates.txt"
cov = fread(cov_file)

## Define Output directory for non-conditioned analysis
nc_out_dir <- "Output_files/non_conditioned/"
## Define Output directory for conditioned analysis
con_out_dir <- "Output_files/conditioned/"




#################################################################################
### LONGS only - DONE (nc: 52/52, con: 43/52)
pheno_file1 <- "Input_data/pheno_post_qc_longs_52_aptamers.txt"

check = fread(pheno_file1)
head(check)

# Non-conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file1,
                # "--pheno-name X25918.60",
                 "--covar", cov_file,
                 "--covar-name age-PC10,LONGS_MAP_660k-LONGS_ADNI_Omni2.5M,LONGS_MAP_OmniEx-LONGS_NIALOAD_Human610_Quadv1",
                 "--out", paste0(nc_out_dir, "LONGS_CSF_pQTL_PheWAS_HP_region"),
                 sep = " ")
system(command)



# Conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file1,
               #  "--pheno-name X25918.60",
                 "--covar", cov_file,
                 "--covar-name age-PC10,LONGS_MAP_660k-LONGS_ADNI_Omni2.5M,LONGS_MAP_OmniEx-LONGS_NIALOAD_Human610_Quadv1,HP1_pred",
                 "--out", paste0(con_out_dir, "LONGS_CSF_pQTL_PheWAS_HP_region_cond_HP1"),
                 sep = " ")
system(command)



#################################################################################
### LONGS + Garfield - DONE (nc: 1725/1725, con: 1449/1725)
pheno_file2 <- "Input_data/pheno_post_qc_longs_garfield_1725_aptamers.txt"


# Non-conditioned - DONE
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file2,
              #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-PC10,LONGS_MAP_660k,LONGS_MAP_660W,LONGS_DIAN_CoreEx-GARFIELD_MAP_NeuroX2",
                 "--out", paste0(nc_out_dir, "LONGS_garfield_CSF_pQTL_PheWAS_HP_region"),
                 sep = " ")
system(command)



# Conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file2,
              #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-PC10,LONGS_MAP_660k,LONGS_MAP_660W,LONGS_DIAN_CoreEx-GARFIELD_MAP_NeuroX2,HP1_pred",
                 "--out", paste0(con_out_dir, "LONGS_garfield_CSF_pQTL_PheWAS_HP_region_cond_HP1"),
                 sep = " ")
system(command)



#################################################################################
### LONGS + PPMI - DONE (nc: 25/25, con: 21/25)
pheno_file3 = "Input_data/pheno_post_qc_longs_ppmi_25_aptamers.txt"


# Non-conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file3,
                 #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-PC10,LONGS_MAP_660k-LONGS_NIALOAD_Human610_Quadv1",
                 "--out", paste0(nc_out_dir, "LONGS_PPMI_CSF_pQTL_PheWAS_HP_region"),
                 sep = " ")
system(command)



# Conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file3,
                 #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-PC10,LONGS_MAP_660k-LONGS_NIALOAD_Human610_Quadv1,HP1_pred",
                 "--out", paste0(con_out_dir, "LONGS_PPMI_CSF_pQTL_PheWAS_HP_region_cond_HP1"),
                 sep = " ")
system(command)



#################################################################################
### LONGS + Stanford - DONE (nc: 15/15, con: 12/15)
pheno_file4 = "Input_data/pheno_post_qc_longs_stanford_15_aptamers.txt"


# Non-conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file4,
                 #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-PC10,LONGS_MAP_660k-LONGS_NIALOAD_Human610_Quadv1",
                 "--out", paste0(nc_out_dir, "LONGS_Stanford_CSF_pQTL_PheWAS_HP_region"),
                 sep = " ")
system(command)



# Conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file4,
                 #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-PC10,LONGS_MAP_660k-LONGS_NIALOAD_Human610_Quadv1,HP1_pred",
                 "--out", paste0(con_out_dir, "LONGS_Stanford_CSF_pQTL_PheWAS_HP_region_cond_HP1"),
                 sep = " ")
system(command)



#################################################################################
### LONGS + Garfield + PPMI - DONE (nc: 559/559, con: 493/559)
pheno_file5 = "Input_data/pheno_post_qc_longs_garfield_ppmi_559_aptamers.txt"


# Non-conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file5,
                 #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-PC10,LONGS_MAP_660k-GARFIELD_MAP_NeuroX2",
                 "--out", paste0(nc_out_dir, "LONGS_Garfield_PPMI_CSF_pQTL_PheWAS_HP_region"),
                 sep = " ")
system(command)



# Conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file5,
                 #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-PC10,LONGS_MAP_660k-HP1_pred",
                 "--out", paste0(con_out_dir, "LONGS_Garfield_PPMI_CSF_pQTL_PheWAS_HP_region_cond_HP1"),
                 sep = " ")
system(command)




#################################################################################
### LONGS + Garfield + Stanford - DONE (nc: 963/963, con: 814/963)
pheno_file6 = "Input_data/pheno_post_qc_longs_garfield_stanford_963_aptamers.txt"


# Non-conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file6,
                 #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-PC10,LONGS_MAP_660k-GARFIELD_MAP_NeuroX2",
                 "--out", paste0(nc_out_dir, "LONGS_Garfield_Stanford_CSF_pQTL_PheWAS_HP_region"),
                 sep = " ")
system(command)



# Conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file6,
                 #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-PC10,LONGS_MAP_660k-HP1_pred",
                 "--out", paste0(con_out_dir, "LONGS_Garfield_Stanford_CSF_pQTL_PheWAS_HP_region_cond_HP1"),
                 sep = " ")
system(command)




#################################################################################
### LONGS + PPMI + Stanford - DONE (nc: 61/61, con: 52/61)
pheno_file7 = "Input_data/pheno_post_qc_longs_ppmi_stanford_61_aptamers.txt"


# Non-conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file7,
                 #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-LONGS_NIALOAD_Human610_Quadv1",
                 "--out", paste0(nc_out_dir, "LONGS_PPMI_Stanford_CSF_pQTL_PheWAS_HP_region"),
                 sep = " ")
system(command)



# Conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file7,
                 #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-LONGS_NIALOAD_Human610_Quadv1,HP1_pred",
                 "--out", paste0(con_out_dir, "LONGS_PPMI_Stanford_CSF_pQTL_PheWAS_HP_region_cond_HP1"),
                 sep = " ")
system(command)




#################################################################################
### All cohorts - DONE (nc; 3608/3608, con: /3608)
pheno_file8 = "Input_data/pheno_post_qc_all_cohorts_3608_aptamers.txt"


# Non-conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file8,
                 #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-GARFIELD_MAP_NeuroX2",
                 "--out", paste0(nc_out_dir, "All_cohorts_CSF_pQTL_PheWAS_HP_region"),
                 sep = " ")
system(command)



# Conditioned
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file8,
                 #   "--pheno-name X4792.51",
                 "--covar", cov_file,
                 "--covar-name age-HP1_pred",
                 "--out", paste0(con_out_dir, "All_cohorts_CSF_pQTL_PheWAS_HP_region_cond_HP1"),
                 sep = " ")
system(command)
