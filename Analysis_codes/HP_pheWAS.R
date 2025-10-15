## PheWAS base script


#################################################################################
### Define consistent input parameters

# Plink2 command
pl <- "plink2"
# Plink file genotype file
in_file <- "/path/to/genotype/file"
# HP region reference file
reg_file <- "/path/to/HP/region/file"
# Covariate file with HP1 predictions
cov_file <- "/path/to/covariate/file"

## Define Output directory for non-conditioned analysis
nc_out_dir <- "/path/to/output/directory/for/non/conditioned/results/"
## Define Output directory for conditioned analysis
con_out_dir <- "/path/to/output/directory/for/conditioned/results/"


#################################################################################
# Non-conditioned on HP1 analyses
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file1,
                 "--covar", cov_file,
                 "--covar-name age-arrayGroup_quad660",
                 "--out", paste0(nc_out_dir, "plasma_pQTL_PheWAS_HP_region"),
                 sep = " ")
system(command)


#################################################################################
# Conditioned on HP1 analyses
command <- paste(pl, 
                 "--bfile", in_file,
                 "--glm hide-covar omit-ref log10 cols=+a1freq",
                 "--covar-variance-standardize",
                 "--keep-allele-order",
                 "--geno 0.1",
                 "--mac 10",
                 "--extract", reg_file,
                 "--pheno", pheno_file1,
                 "--covar", cov_file,
                 "--covar-name age-HP1_pred",
                 "--out", paste0(con_out_dir, "plasma_pQTL_PheWAS_HP_region_cond_HP1"),
                 sep = " ")
system(command)
