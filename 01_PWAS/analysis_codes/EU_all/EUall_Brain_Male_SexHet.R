## make tables of brain (non-sex-strat W22 and sex-strat W23) AD PWAS sex-specific genes including GWAS stats (beta, p, sex-heterogeneity)
###################### Danielle Reid - 12/11/24
#######################################################

#load the packages
library(plyr); library(dplyr)
library(readr)
library(data.table)
library(tidyr)
library(methods)
library(optparse)

#############################################################
option_list = list(
  make_option("--PD_dir", action="store", default=NA, type='character',
              help="Path to primary discovery base directory [required]"),
  make_option("--SD_dir", action="store", default=NA, type='character',
              help="Path to secondary discovery base directory [required]"),        
  make_option("--gene_list", action="store", default=NA, type='character',
              help="Filename of sex-specific gene list [required]"),
  make_option("--female_GWAS", action="store", default=NA, type='character',
              help="Path to female GWAS summary statistics [required]"),
  make_option("--male_GWAS", action="store", default=NA, type='character',
              help="Path to male GWAS summary statistics [required]"),
  make_option("--female_file", action="store", default=NA, type='character',
              help="Female file name [required]"),
  make_option("--male_file", action="store", default=NA, type='character',
              help="Male file name [required]"),
  make_option("--file_name", action="store", default=NA, type='character',
              help="Output filename [required]")            
)

opt = parse_args(OptionParser(option_list=option_list))

#############################################################

# create variables for base directory
PD_dir = opt$PD_dir; print(paste("primary discovery directory:", PD_dir))
SD_dir = opt$SD_dir; print(paste("secondary discovery directory:", SD_dir))

# create variables for sex-stratified PWAS directories and file base names
f_dir = "Female"; print(paste("constructed F-F directory: ", file.path(PD_dir, f_dir)))
m_dir = "Male"; print(paste("constructed M-M directory: ", file.path(PD_dir, m_dir)))
f_f = opt$female_file; print(paste("female file name:", f_f))
m_f = opt$male_file; print(paste("male file name:", m_f))

#load gene list, sex-strat PWASs per discovery, and GWAS summary statistics
gene_list = fread(opt$gene_list); gene_list = as.data.frame(gene_list); print("gene list loaded")
f_PWAS = fread(paste0(file.path(PD_dir, f_dir, f_f), "sex-strat-W23_weights.txt")); f_PWAS = as.data.frame(f_PWAS); print("Female Primary Discovery loaded")
m_PWAS = fread(paste0(file.path(PD_dir, m_dir, m_f), "sex-strat-W23_weights.txt")); m_PWAS = as.data.frame(m_PWAS); print("Male Primary Discovery loaded")
f_SD_PWAS = fread(paste0(file.path(SD_dir, f_dir, f_f), "non-sex-strat-W23_weights.txt")); f_SD_PWAS = as.data.frame(f_SD_PWAS); print("Female Secondary Discovery loaded")
m_SD_PWAS = fread(paste0(file.path(SD_dir, m_dir, m_f), "non-sex-strat-W23_weights.txt")); m_SD_PWAS = as.data.frame(m_SD_PWAS); print("Male Secondary Discovery loaded")

f_stats = fread(opt$female_GWAS); f_stats = as.data.frame(f_stats); print("Female GWAS summary statistics loaded")
m_stats = fread(opt$male_GWAS); m_stats = as.data.frame(m_stats); print("Male GWAS summary statistics loaded")

# add ENSG without gene name as a column
gene_list$ENSG = gsub("\\..*$", "", as.character(gene_list$ENSG_ID))

#create dataframe of ENSG_ID, ENSG, and gene name for sex-specific genes
a = subset(gene_list, topSNP != 0, select = c(ENSG_ID, ENSG, ID, NEW_hit, LOC_con, LOC_stu, CHR, P0))
a = a[order(a$CHR, a$P0), ]
names(a)[3] = "Gene"

# add sex-matched PWAS results from sex that sex-specific genes were identified
a$Male.PD.PWAS.Z = ifelse(a$Gene %in% m_PWAS$ID, m_PWAS$TWAS.Z[match(a$Gene, m_PWAS$ID)], NA)
a$Male.PD.PWAS.P = ifelse(a$Gene %in% m_PWAS$ID, m_PWAS$TWAS.P[match(a$Gene, m_PWAS$ID)], NA)
a$Male.PD.GWAS.SNP = ifelse(a$Gene %in% m_PWAS$ID, m_PWAS$BEST.GWAS.ID[match(a$Gene, m_PWAS$ID)], NA)
a$Male.PD.pQTL.SNP = ifelse(a$Gene %in% m_PWAS$ID, m_PWAS$EQTL.ID[match(a$Gene, m_PWAS$ID)], NA)

# add secondary discovery PWAS results from sex that sex-specific genes were identified
a$Male.SD.PWAS.Z = ifelse(a$Gene %in% m_SD_PWAS$ID, m_SD_PWAS$TWAS.Z[match(a$Gene, m_SD_PWAS$ID)], NA)
a$Male.SD.PWAS.P = ifelse(a$Gene %in% m_SD_PWAS$ID, m_SD_PWAS$TWAS.P[match(a$Gene, m_SD_PWAS$ID)], NA)
a$Male.SD.GWAS.SNP = ifelse(a$Gene %in% m_SD_PWAS$ID, m_SD_PWAS$BEST.GWAS.ID[match(a$Gene, m_SD_PWAS$ID)], NA)
a$Male.SD.pQTL.SNP = ifelse(a$Gene %in% m_SD_PWAS$ID, m_SD_PWAS$EQTL.ID[match(a$Gene, m_SD_PWAS$ID)], NA)

# make a new sexx-specific GWAS.SNP column so that GWAS.SNP that were only from one discovery are included altogether
a$Male.GWAS.SNP = ifelse(is.na(a$Male.PD.GWAS.SNP), a$Male.SD.GWAS.SNP, a$Male.PD.GWAS.SNP)

# add sex-strat GWAS stats for each GWAS.SNP
a$Female.GWAS.Beta = ifelse(a$Male.GWAS.SNP %in% f_stats$SNP, f_stats$BETA[match(a$Male.GWAS.SNP, f_stats$SNP)], NA)
a$Female.GWAS.SE = ifelse(a$Male.GWAS.SNP %in% f_stats$SNP, f_stats$SE[match(a$Male.GWAS.SNP, f_stats$SNP)], NA)
a$Female.GWAS.P = ifelse(a$Male.GWAS.SNP %in% f_stats$SNP, f_stats$P[match(a$Male.GWAS.SNP, f_stats$SNP)], NA)

a$Male.GWAS.Beta = ifelse(a$Male.GWAS.SNP %in% m_stats$SNP, m_stats$BETA[match(a$Male.GWAS.SNP, m_stats$SNP)], NA)
a$Male.GWAS.SE = ifelse(a$Male.GWAS.SNP %in% m_stats$SNP, m_stats$SE[match(a$Male.GWAS.SNP, m_stats$SNP)], NA)
a$Male.GWAS.P = ifelse(a$Male.GWAS.SNP %in% m_stats$SNP, m_stats$P[match(a$Male.GWAS.SNP, m_stats$SNP)], NA)

# calculate sex heterogeneity for GWAS.SNP
a$GWAS.Gender.Z <- ifelse(!is.na(a$Female.GWAS.Beta) | !is.na(a$Male.GWAS.Beta), 
                         (a$Male.GWAS.Beta - a$Female.GWAS.Beta) / sqrt((a$Male.GWAS.SE)^2 + (a$Female.GWAS.SE)^2), 
                        NA
)

a$GWAS.Gender.P <- ifelse(!is.na(a$GWAS.Gender.Z), 2*pnorm(q=abs(a$GWAS.Gender.Z), lower.tail=FALSE), NA) 

#add sex-matched PWAS results from opposite sex that sex-specific genes were identified
a$Female.PD.PWAS.Z = ifelse(a$Gene %in% f_PWAS$ID, f_PWAS$TWAS.Z[match(a$Gene, f_PWAS$ID)], NA)
a$Female.PD.PWAS.P = ifelse(a$Gene %in% f_PWAS$ID, f_PWAS$TWAS.P[match(a$Gene, f_PWAS$ID)], NA)
a$Female.PD.GWAS.SNP = ifelse(a$Gene %in% f_PWAS$ID, f_PWAS$BEST.GWAS.ID[match(a$Gene, f_PWAS$ID)], NA)
a$Female.PD.pQTL.SNP = ifelse(a$Gene %in% f_PWAS$ID, f_PWAS$EQTL.ID[match(a$Gene, f_PWAS$ID)], NA)

# add secondary discovery PWAS results from opposite sex that sex-specific genes were identified
a$Female.SD.PWAS.Z = ifelse(a$Gene %in% f_SD_PWAS$ID, f_SD_PWAS$TWAS.Z[match(a$Gene, f_SD_PWAS$ID)], NA)
a$Female.SD.PWAS.P = ifelse(a$Gene %in% f_SD_PWAS$ID, f_SD_PWAS$TWAS.P[match(a$Gene, f_SD_PWAS$ID)], NA)
a$Female.SD.GWAS.SNP = ifelse(a$Gene %in% f_SD_PWAS$ID, f_SD_PWAS$BEST.GWAS.ID[match(a$Gene, f_SD_PWAS$ID)], NA)
a$Female.SD.pQTL.SNP = ifelse(a$Gene %in% f_SD_PWAS$ID, f_SD_PWAS$EQTL.ID[match(a$Gene, f_SD_PWAS$ID)], NA)

fwrite(a, paste0(file.path(PD_dir, m_dir), "/", opt$file_name), na = NA, col.names = T, sep = '\t')



