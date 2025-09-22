## make tables of brain (non-sex-strat and sex-strat NGI CSF) AD PWAS sex-specific genes including GWAS stats (beta, p, sex-heterogeneity)
###################### Danielle Reid - 12/20/24
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
  make_option("--female_sGWAS", action="store", default=NA, type='character',
              help="Path to female NGI sex-strat GWAS summary statistics [required]"),
  make_option("--male_sGWAS", action="store", default=NA, type='character',
              help="Path to male NGI sex-strat GWAS summary statistics [required]"),
  make_option("--female_nGWAS", action="store", default=NA, type='character',
              help="Path to female NGI non-sex-strat GWAS summary statistics [required]"),
  make_option("--male_nGWAS", action="store", default=NA, type='character',
              help="Path to male NGI non-sex-strat GWAS summary statistics [required]"),            
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
f_PWAS = fread(paste0(file.path(PD_dir, f_dir, f_f), "wHP_noMHC_ext2Mb_sex-strat-CSFcis_weights.txt")); f_PWAS = as.data.frame(f_PWAS); print("Female Primary Discovery loaded")
m_PWAS = fread(paste0(file.path(PD_dir, m_dir, m_f), "wHP_noMHC_ext2Mb_sex-strat-CSFcis_weights.txt")); m_PWAS = as.data.frame(m_PWAS); print("Male Primary Discovery loaded")
f_SD_PWAS = fread(paste0(file.path(SD_dir, f_dir, f_f), "wHP_noMHC_ext2Mb_non-sex-strat-CSFcis_weights.txt")); f_SD_PWAS = as.data.frame(f_SD_PWAS); print("Female Secondary Discovery loaded")
m_SD_PWAS = fread(paste0(file.path(SD_dir, m_dir, m_f), "wHP_noMHC_ext2Mb_non-sex-strat-CSFcis_weights.txt")); m_SD_PWAS = as.data.frame(m_SD_PWAS); print("Male Secondary Discovery loaded")

fs_stats = fread(opt$female_sGWAS); fs_stats = as.data.frame(fs_stats); print("Female NGI sex-strat GWAS summary statistics loaded")
ms_stats = fread(opt$male_sGWAS); ms_stats = as.data.frame(ms_stats); print("Male NGI sex-strat GWAS summary statistics loaded")

fn_stats = fread(opt$female_nGWAS); fn_stats = as.data.frame(fn_stats); print("Female NGI non-sex-strat GWAS summary statistics loaded")
mn_stats = fread(opt$male_nGWAS); mn_stats = as.data.frame(mn_stats); print("Male NGI non-sex-strat GWAS summary statistics loaded")

# prep GWAS summary statistics
names(fs_stats)[2] = "rsID"; names(ms_stats)[2] = "rsID"; names(fn_stats)[2] = "rsID"; names(mn_stats)[2] = "rsID"

fs_stats = select(fs_stats, -posID, -tmpID, -META_DIR, -Q_V, -Q_P)
ms_stats = select(ms_stats, -posID, -tmpID, -META_DIR, -Q_V, -Q_P)
fn_stats = select(fn_stats, -posID, -tmpID, -META_DIR, -Q_V, -Q_P)
mn_stats = select(mn_stats, -posID, -tmpID, -META_DIR, -Q_V, -Q_P)

names(fs_stats)[10] = "N"; names(ms_stats)[10] = "N"; names(fn_stats)[10] = "N"; names(mn_stats)[10] = "N"

f_stats = merge(fs_stats, fn_stats, by = c("rsID", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "P", "N", "SNP"), all = TRUE)
m_stats = merge(ms_stats, mn_stats, by = c("rsID", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "P", "N", "SNP"), all = TRUE)

#create dataframe of analyte, gene, CHR, TSS, etc for sex-specific genes
a = subset(gene_list, topSNP != 0, select = c(ID, EntrezGeneSymbol, Gene, NSNP, NWGT, NEW_hit, LOC_con, LOC_stu, CHR, TSS))
a = a[order(a$CHR, a$TSS), ]

# Add sex-matched PWAS results from sex that sex-specific genes were identified
a <- a %>%
  left_join(m_PWAS %>%
              select(EntrezGeneSymbol, ID, TWAS.Z, TWAS.P, BEST.GWAS.ID, EQTL.ID),
            by = c("EntrezGeneSymbol", "ID")) %>%
  rename(
    Male.PD.PWAS.Z = TWAS.Z,
    Male.PD.PWAS.P = TWAS.P,
    Male.PD.GWAS.SNP = BEST.GWAS.ID,
    Male.PD.pQTL.SNP = EQTL.ID
)

# Add secondary discovery PWAS results from sex that sex-specific genes were identified
a <- a %>%
  left_join(m_SD_PWAS %>%
              select(EntrezGeneSymbol, ID, TWAS.Z, TWAS.P, BEST.GWAS.ID, EQTL.ID),
            by = c("EntrezGeneSymbol", "ID")) %>%
  rename(
    Male.SD.PWAS.Z = TWAS.Z,
    Male.SD.PWAS.P = TWAS.P,
    Male.SD.GWAS.SNP = BEST.GWAS.ID,
    Male.SD.pQTL.SNP = EQTL.ID
)

# add sex-strat GWAS stats for each PD.GWAS.SNP
a$Female.PD.GWAS.Beta = ifelse(a$Male.PD.GWAS.SNP %in% f_stats$SNP, f_stats$BETA[match(a$Male.PD.GWAS.SNP, f_stats$SNP)], NA)
a$Female.PD.GWAS.SE = ifelse(a$Male.PD.GWAS.SNP %in% f_stats$SNP, f_stats$SE[match(a$Male.PD.GWAS.SNP, f_stats$SNP)], NA)
a$Female.PD.GWAS.P = ifelse(a$Male.PD.GWAS.SNP %in% f_stats$SNP, f_stats$P[match(a$Male.PD.GWAS.SNP, f_stats$SNP)], NA)

a$Male.PD.GWAS.Beta = ifelse(a$Male.PD.GWAS.SNP %in% m_stats$SNP, m_stats$BETA[match(a$Male.PD.GWAS.SNP, m_stats$SNP)], NA)
a$Male.PD.GWAS.SE = ifelse(a$Male.PD.GWAS.SNP %in% m_stats$SNP, m_stats$SE[match(a$Male.PD.GWAS.SNP, m_stats$SNP)], NA)
a$Male.PD.GWAS.P = ifelse(a$Male.PD.GWAS.SNP %in% m_stats$SNP, m_stats$P[match(a$Male.PD.GWAS.SNP, m_stats$SNP)], NA)

# calculate sex heterogeneity for GWAS.SNP
a$PD.GWAS.Gender.Z <- ifelse(!is.na(a$Female.PD.GWAS.Beta) | !is.na(a$Male.PD.GWAS.Beta), 
                         (a$Male.PD.GWAS.Beta - a$Female.PD.GWAS.Beta) / sqrt((a$Male.PD.GWAS.SE)^2 + (a$Female.PD.GWAS.SE)^2), 
                        NA
)

a$PD.GWAS.Gender.P <- ifelse(!is.na(a$PD.GWAS.Gender.Z), 2*pnorm(q=abs(a$PD.GWAS.Gender.Z), lower.tail=FALSE), NA) 

# add sex-strat GWAS stats for each SD.GWAS.SNP
a$Female.SD.GWAS.Beta = ifelse(a$Male.SD.GWAS.SNP %in% f_stats$SNP, f_stats$BETA[match(a$Male.SD.GWAS.SNP, f_stats$SNP)], NA)
a$Female.SD.GWAS.SE = ifelse(a$Male.SD.GWAS.SNP %in% f_stats$SNP, f_stats$SE[match(a$Male.SD.GWAS.SNP, f_stats$SNP)], NA)
a$Female.SD.GWAS.P = ifelse(a$Male.SD.GWAS.SNP %in% f_stats$SNP, f_stats$P[match(a$Male.SD.GWAS.SNP, f_stats$SNP)], NA)

a$Male.SD.GWAS.Beta = ifelse(a$Male.SD.GWAS.SNP %in% m_stats$SNP, m_stats$BETA[match(a$Male.SD.GWAS.SNP, m_stats$SNP)], NA)
a$Male.SD.GWAS.SE = ifelse(a$Male.SD.GWAS.SNP %in% m_stats$SNP, m_stats$SE[match(a$Male.SD.GWAS.SNP, m_stats$SNP)], NA)
a$Male.SD.GWAS.P = ifelse(a$Male.SD.GWAS.SNP %in% m_stats$SNP, m_stats$P[match(a$Male.SD.GWAS.SNP, m_stats$SNP)], NA)

# calculate sex heterogeneity for GWAS.SNP
a$SD.GWAS.Gender.Z <- ifelse(!is.na(a$Female.SD.GWAS.Beta) | !is.na(a$Male.SD.GWAS.Beta), 
                         (a$Male.SD.GWAS.Beta - a$Female.SD.GWAS.Beta) / sqrt((a$Male.SD.GWAS.SE)^2 + (a$Female.SD.GWAS.SE)^2), 
                        NA
)

a$SD.GWAS.Gender.P <- ifelse(!is.na(a$SD.GWAS.Gender.Z), 2*pnorm(q=abs(a$SD.GWAS.Gender.Z), lower.tail=FALSE), NA) 

# Add sex-matched PWAS results from opposite sex that sex-specific genes were identified
a <- a %>%
  left_join(f_PWAS %>%
              select(EntrezGeneSymbol, ID, TWAS.Z, TWAS.P, BEST.GWAS.ID, EQTL.ID),
            by = c("EntrezGeneSymbol", "ID")) %>%
  rename(
    Female.PD.PWAS.Z = TWAS.Z,
    Female.PD.PWAS.P = TWAS.P,
    Female.PD.GWAS.SNP = BEST.GWAS.ID,
    Female.PD.pQTL.SNP = EQTL.ID
)

# Add secondary discovery PWAS results from opposite sex that sex-specific genes were identified
a <- a %>%
  left_join(f_SD_PWAS %>%
              select(EntrezGeneSymbol, ID, TWAS.Z, TWAS.P, BEST.GWAS.ID, EQTL.ID),
            by = c("EntrezGeneSymbol", "ID")) %>%
  rename(
    Female.SD.PWAS.Z = TWAS.Z,
    Female.SD.PWAS.P = TWAS.P,
    Female.SD.GWAS.SNP = BEST.GWAS.ID,
    Female.SD.pQTL.SNP = EQTL.ID
)

fwrite(a, paste0(file.path(PD_dir, m_dir), "/", opt$file_name), na = NA, col.names = T, sep = '\t')



