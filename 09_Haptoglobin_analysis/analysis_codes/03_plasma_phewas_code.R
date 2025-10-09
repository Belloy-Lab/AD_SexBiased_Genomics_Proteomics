## HP plasma pheWAS Code

library(data.table)
library(tidyverse)



################
# Add HP1 pred to covariate file

# Read in HP1 prediction file - rename id variable to match covariate id variable name
hp1 = fread("Input_data/HP1_predictions_f05_v2.csv") %>% 
  dplyr::rename(IID = sample_id)
nrow(hp1)
# N = 2297

# Read in plasma covariate file
cov1 = fread("Input_data/f14_covar_plasmaEURpqtl.txt") 

# Merge covariate file with HP1 prediction file
cov_hp1 = merge(cov1, hp1, by = "IID") %>% 
  dplyr::rename(fid_old = FID,
                FID = IID) %>% 
  dplyr::rename(IID = fid_old)
nrow(cov_hp1)
# N = 2297

# write out new covariate file
fwrite(cov_hp1, "Input_data/f14_covar_plasmaEURpqtl_with_HP1.txt", sep = "\t", quote = F)


#################################################################################
### Define consistent input parameters

# Plink2 command
pl <- "plink2"
# Plink file genotype file
in_file <- "Input_data/f16_genoEUR_HP_cisphewas"
# HP region reference file
reg_file <- "Input_data/HP_region_variants_1MB_radius.txt"
# Covariate file with HP1 predictions
cov_file <- "Input_data/f14_covar_plasmaEURpqtl_with_HP1.txt"
cov_file2 <- "Input_data/f14_covar_plasmaEURpqtl.txt"

## Define Output directory for non-conditioned analysis
nc_out_dir <- "Input_data/output_data/non_conditioned/"
## Define Output directory for conditioned analysis
con_out_dir <- "Input_data/output_data/conditioned/"




#################################################################################
### All Analytes
pheno_file1 <- "Input_data/f13_pheno_all6907proteins_plasmaEURpqtl.txt"

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
               #  "--pheno-name X9997.12",
                 "--covar", cov_file,
                 "--covar-name age-arrayGroup_quad660",
                 "--out", paste0(nc_out_dir, "plasma_pQTL_PheWAS_HP_region"),
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
              #   "--pheno-name X2182.54",
                 "--covar", cov_file,
                 "--covar-name age-HP1_pred",
                 "--out", paste0(con_out_dir, "plasma_pQTL_PheWAS_HP_region_cond_HP1"),
                 sep = " ")
system(command)
