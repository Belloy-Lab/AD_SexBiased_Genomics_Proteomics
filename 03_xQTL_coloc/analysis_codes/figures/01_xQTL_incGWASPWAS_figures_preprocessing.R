### Preprocessing Output data from xQTL analysis to create figure

library(data.table)
library(tidyverse)
library(patchwork)
library(grid)
library(reshape2)
library(gridExtra)
library(scales)
library(optparse)
option_list <- list(
  make_option("--work_dir", type = "character", help = "Working directory"),
  make_option("--abf_dir", type = "character", help = "Directory of ABF CSV files"),
  make_option("--adj_abf_dir", type = "character", help = "Directory of adjusted ABF CSV files"),
  make_option("--susie_dir", type = "character", help = "Directory of SuSiE CSV files"),
  make_option("--out_GP_qtl", type = "character", help = "Output CSV for annotated GWAS and PWAS combined results"),
  make_option("--out_pwas", type = "character", help = "Output CSV for PWAS subset"),
  make_option("--out_gwas", type = "character", help = "Output CSV for GWAS subset"),
  make_option("--out_mqtl", type = "character", help = "Output CSV for mQTL subset"),
  make_option("--out_haqtl", type = "character", help = "Output CSV for haQTL subset"),
  make_option("--out_caqtl", type = "character", help = "Output CSV for caQTL subset")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$work_dir)) stop("--work_dir is required")
dir.create(opt$work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$work_dir)


###############################################################################################################################################
## Append Default ABF files
abf_dir = file.path(opt$work_dir, opt$abf_dir)
abf_files = list.files(abf_dir, pattern = ".csv$", full.names = T)

abf = data.frame()
for (file in abf_files){
  
  df = fread(file)
  
  abf = rbind(abf, df)
}

n_distinct(abf$locus) # 57
table(abf$discovery, abf$locus) # GBA, INPP5D, and TSPAN14 observed in both GWAS and PWAS


# Create new variables for further preprocessing steps
abf = abf %>%
  arrange(desc(PP4)) %>% 
  dplyr::rename(abf_PP4 = PP4) %>% 
  mutate(id = paste(dset, locus, gene_name, molecular_trait_id, qtl_type, discovery, sep = "-:-"),
         ad_sample = str_extract(dset, "[^_]+$"))


###############################################################################################################################################
## Append Adjusted ABF files
adj_dir = file.path(opt$work_dir, opt$adj_abf_dir) # "/storage2/fs1/belloy2/Active/05_Projects/sivas/xQTL/Adj_ABF/"
adj_files = list.files(adj_dir, pattern = ".csv$", full.names = T)

adj = data.frame()
for (file in adj_files){
  
  df = fread(file)
  
  adj = rbind(adj, df)
}

n_distinct(adj$locus) # 40 - exclused GWAS only loci
table(adj$discovery, adj$locus) 

# Create new variables for further preprocessing steps - filter to only include gene matching locus for adjusted
adj = adj %>%
  filter(locus == gene_name) %>% 
  arrange(desc(PP4)) %>% 
  dplyr::rename(adj_PP4 = PP4) %>% 
  mutate(dset = gsub("_boosted", "", dset),
         id = paste(dset, locus, gene_name, molecular_trait_id, qtl_type, discovery, sep = "-:-"),
         ad_sample = str_extract(dset, "[^_]+$"))


###############################################################################################################################################
## Append SuSiE files
susie_dir = file.path(opt$work_dir, opt$susie_dir)
susie_files = list.files(susie_dir, pattern = ".csv$", full.names = T)

susie = data.frame()
for (file in susie_files){
  
  df = fread(file)
  
  susie = rbind(susie, df)
}

n_distinct(susie$locus) # 57
table(susie$discovery, susie$locus) 

# Create new variables for further preprocessing steps
susie = susie %>%
  arrange(desc(PP4)) %>% 
  dplyr::rename(susie_PP4 = PP4) %>% 
  mutate(id = paste(dset, locus, gene_name, molecular_trait_id, qtl_type, discovery, sep = "-:-"),
         ad_sample = str_extract(dset, "[^_]+$"))


###############################################################################################################################################
### Merge appended files together based on id variable - extract meaningful variables and rename
## Combine SuSiE and ABF Files
abf_susie = full_join(abf, susie, by = "id") %>% 
  mutate(susie_PP4 = ifelse(is.na(susie_PP4), 0, susie_PP4),
         abf_PP4 = ifelse(is.na(abf_PP4), 0, abf_PP4),
         best_PP4 = pmax(susie_PP4, abf_PP4),
         best_method = ifelse(susie_PP4 >= abf_PP4, "susie", "abf"),
         n_snps = ifelse(is.na(n_snps.x), n_snps.y, n_snps.x)) %>% 
  separate(id, into = c("dset", "others"), sep = "-:-", extra = "merge") %>%
  separate(others, into = c("locus", "gene_name", "molecular_trait_id2", "qtl_type", "discovery"), sep = "-:-") %>% 
  dplyr::rename(hit1 = hit1.y, hit2 = hit2.y) %>% 
  mutate(molecular_trait_id.y = as.character(molecular_trait_id.y),
         ad_sample = case_when(!is.na(ad_sample.x) ~ ad_sample.x,
                               T ~ ad_sample.y),
         molecular_trait_id = case_when(!is.na(molecular_trait_id.x) ~ molecular_trait_id.x,
                                        T ~ molecular_trait_id.y),
         CHR = case_when(!is.na(CHROM.x) ~ CHROM.x,
                         T ~ CHROM.y),
         BP = case_when(!is.na(BP.x) ~ BP.x,
                        T ~ BP.y),
         id = paste(dset, locus, gene_name, molecular_trait_id, qtl_type, discovery, sep = "-:-")) %>% 
  dplyr::select(CHR, BP, dset, locus, gene_name, molecular_trait_id, abf_PP4, susie_PP4, best_PP4, best_method, qtl_type, hit1, hit2, n_snps, ad_sample, discovery, id)

# Run checks
n_distinct(abf_susie$locus) # 57
table(abf_susie$discovery, abf_susie$locus) 


#######################################################
## Combine adjusted abf values to abf_susie merged file
abf_susie2 = full_join(abf_susie, adj, by = "id") %>% 
  mutate(adj_PP4 = ifelse(is.na(adj_PP4), 0, adj_PP4),
         best_PP4 = pmax(best_PP4, adj_PP4),
         best_method = ifelse(adj_PP4 >= best_PP4, "adj", best_method)) %>% 
  separate(id, into = c("dset", "others"), sep = "-:-", extra = "merge") %>%
  separate(others, into = c("locus", "gene_name", "molecular_trait_id2", "qtl_type", "discovery"), sep = "-:-") %>% 
  dplyr::rename(hit1 = hit1.x, hit2 = hit2.x, n_snps = n_snps.x) %>% 
  mutate(ad_sample = case_when(!is.na(ad_sample.x) ~ ad_sample.x,
                               T ~ ad_sample.y),
         molecular_trait_id = case_when(!is.na(molecular_trait_id.x) ~ molecular_trait_id.x,
                                        T ~ molecular_trait_id.y),
         BP = case_when(!is.na(BP.x) ~ BP.x,
                        T ~ BP.y)) %>% 
  dplyr::select(CHR, BP, dset, locus, gene_name, molecular_trait_id, abf_PP4, susie_PP4, adj_PP4, best_PP4, best_method, qtl_type, hit1, hit2, n_snps, ad_sample, discovery)

# Run checks
n_distinct(abf_susie2$locus) # 57
table(abf_susie2$discovery, abf_susie2$locus) 





###############################################################################################################################################
### Seperate by QTL type - create tissue_type variable for table/figure

#######################################################
## eQTL
eqtl = abf_susie2 %>% 
  filter(qtl_type == "eQTL") %>% 
  mutate(study = str_extract(dset, ".*(?=_chr_)"),
         study = case_when(
           grepl("^eQTL_catalogue_quach_Pam3CSK", study) ~ "eQTL_catalogue_quach_Pam3CSK4",
           TRUE ~ study)) %>% 
  arrange(study)
unique(eqtl$study) # 87 unique eQTL datasets

# Create tissue_type variable based on index position of study variable - tissue_type variable used to create groups for final figure
eqtl = eqtl %>% 
  mutate(tissue_type = case_when(
    study == unique(eqtl$study)[1] ~ "BrainMeta (Multiple regions)",
    study %in% unique(eqtl$study)[2:6] ~ "MetaBrain Brain (5 regions)",
    study %in% unique(eqtl$study)[7:9] ~ "GTEx_v8 Other (35 tissues)",
    study %in% unique(eqtl$study)[10:11] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(eqtl$study)[12:13] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(eqtl$study)[14] ~ "GTEx_v8 Blood",
    study == unique(eqtl$study)[15] ~ "Blueprint Monocytes",
    study %in% unique(eqtl$study)[16:17] ~ "Braineac (2 regions)",
    study == unique(eqtl$study)[18] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(eqtl$study)[19] ~ "GTEx_v8 Brain (13 regions)",
    study == unique(eqtl$study)[20] ~ "CEDAR Monocytes", 
    study %in% unique(eqtl$study)[21:22] ~ "GTEx_v8 Brain (13 regions)",
    study == unique(eqtl$study)[23] ~ "CommonMind DLPFC", 
    study == unique(eqtl$study)[24] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(eqtl$study)[25] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(eqtl$study)[26:28] ~ "GTEx_v8 Other (35 tissues)",
    study %in% unique(eqtl$study)[29:32] ~ "Fairfax 2014 Monocytes (4 conditions)",
    study == unique(eqtl$study)[33] ~ "GTEx_v8 Other (35 tissues)",
    study %in% unique(eqtl$study)[34:36] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(eqtl$study)[37:38] ~ "Kasela 2017 T-cells (2 types)",
    study %in% unique(eqtl$study)[39:44] ~ "GTEx_v8 Other (35 tissues)",
    study %in% unique(eqtl$study)[45:47] ~ "Nedelec Macrophages (3 conditions)",
    study == unique(eqtl$study)[48] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(eqtl$study)[49:52] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(eqtl$study)[53] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(eqtl$study)[54:58] ~ "Quach Monocytes (5 conditions)",
    study %in% unique(eqtl$study)[59:63] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(eqtl$study)[64] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(eqtl$study)[65:66] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(eqtl$study)[67] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(eqtl$study)[68:74] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(eqtl$study)[75] ~ "Young 2019 Microglia",
    study == unique(eqtl$study)[76] ~ "eQTLgen Blood",
    study == unique(eqtl$study)[77] ~ "Fujita 2024 Astrocytes",
    study == unique(eqtl$study)[78] ~ "Fujita 2024 Endothelial Cells",
    study == unique(eqtl$study)[79] ~ "Fujita 2024 Excitatory Neurons",
    study == unique(eqtl$study)[80] ~ "Fujita 2024 Inhibitory Neurons",
    study == unique(eqtl$study)[81] ~ "Fujita 2024 Microglia",
    study == unique(eqtl$study)[82] ~ "Fujita 2024 Oligodendrocytes",
    study == unique(eqtl$study)[83] ~ "Fujita 2024 Oligodendrocyte Progenitor Cells",
    study == unique(eqtl$study)[84] ~ "Kosoy 2022 Microglia", 
    study == unique(eqtl$study)[85] ~ "Wingo (All) DLPFC", 
    study == unique(eqtl$study)[86] ~ "Wingo (Female) DLPFC",
    study == unique(eqtl$study)[87] ~ "Wingo (Male) DLPFC"))


#######################################################
## pQTL
pqtl = abf_susie2 %>% 
  filter(qtl_type == "pQTL") %>% 
  mutate(study = str_extract(dset, ".*(?=_chr_)")) %>% 
  arrange(study)
unique(pqtl$study) # 8 unique pQTL datasets

# Create tissue_type variable based on index position of study variable - tissue_type variable used to create groups for final figure
pqtl = pqtl %>% 
  mutate(tissue_type = case_when(
    study == unique(pqtl$study)[1] ~ "ARIC Plasma", 
    study == unique(pqtl$study)[2] ~ "NGI (All) CSF", 
    study == unique(pqtl$study)[3] ~ "NGI (Female) CSF", 
    study == unique(pqtl$study)[4] ~ "NGI (Male) CSF", 
    study == unique(pqtl$study)[5] ~ "UKB Plasma", 
    study == unique(pqtl$study)[6] ~ "Wingo (All) DLPFC", 
    study == unique(pqtl$study)[7] ~ "Wingo (Female) DLPFC", 
    study == unique(pqtl$study)[8] ~ "Wingo (Male) DLPFC"))


#######################################################
## sQTL
sqtl = abf_susie2 %>% 
  filter(qtl_type == "sQTL") %>% 
  mutate(study = str_extract(dset, ".*(?=_chr_)")) %>% 
  arrange(study)
unique(sqtl$study) # 62 unique sQTL datasets

# Create tissue_type variable based on index position of study variable - tissue_type variable used to create groups for final figure
sqtl = sqtl %>% 
  mutate(tissue_type = case_when(
    study %in% unique(sqtl$study)[1:3] ~ "GTEx_v8 Other (35 tissues)",
    study %in% unique(sqtl$study)[4:5] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(sqtl$study)[6:7] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(sqtl$study)[8] ~ "GTEx_v8 Blood",
    study == unique(sqtl$study)[9] ~ "Blueprint Monocytes",
    study %in% unique(sqtl$study)[10:11] ~ "Braineac (2 regions)",
    study == unique(sqtl$study)[12] ~ "GTEx_v8 Other (35 tissues)",
    study %in% unique(sqtl$study)[13:15] ~ "GTEx_v8 Brain (13 regions)",
    study == unique(sqtl$study)[16] ~ "CommonMind DLPFC",
    study == unique(sqtl$study)[17] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(sqtl$study)[18] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(sqtl$study)[19:22] ~ "GTEx_v8 Other (35 tissues)",
    study %in% unique(sqtl$study)[23:25] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(sqtl$study)[26:31] ~ "GTEx_v8 Other (35 tissues)",
    study %in% unique(sqtl$study)[32:34] ~ "Nedelec Macrophages (3 conditions)",
    study == unique(sqtl$study)[35] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(sqtl$study)[36:39] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(sqtl$study)[40] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(sqtl$study)[41:45] ~ "Quach Monocytes (5 conditions)",
    study %in% unique(sqtl$study)[46:50] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(sqtl$study)[51] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(sqtl$study)[52:53] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(sqtl$study)[54] ~ "GTEx_v8 Brain (13 regions)",
    study %in% unique(sqtl$study)[55:61] ~ "GTEx_v8 Other (35 tissues)",
    study == unique(sqtl$study)[62] ~ "Young 2019 Microglia"))


#######################################################
## mQTL
mqtl = abf_susie2 %>% 
  filter(qtl_type == "mQTL") %>% 
  mutate(study = str_extract(dset, ".*(?=_chr_)")) %>% 
  arrange(study)
unique(mqtl$study) # 2 unique mQTL datasets

# Create tissue_type variable based on index position of study variable - tissue_type variable used to create groups for final figure
mqtl = mqtl %>% 
  mutate(tissue_type = case_when(
    study == unique(mqtl$study)[1] ~ "BrainMeta (Multiple regions)",
    study == unique(mqtl$study)[2] ~ "xQTLServe DLPFC"))


#######################################################
## haQTL
haqtl = abf_susie2 %>% 
  filter(qtl_type == "haQTL") %>% 
  mutate(study = str_extract(dset, ".*(?=_chr_)")) %>% 
  arrange(study)
unique(haqtl$study) # 1 unique haQTL datasets

# Create tissue_type variable based on index position of study variable - tissue_type variable used to create groups for final figure
haqtl = haqtl %>% 
  mutate(tissue_type = "xQTLServe DLPFC")


#######################################################
## caQTL
caqtl = abf_susie2 %>% 
  filter(qtl_type == "caQTL") %>% 
  mutate(study = str_extract(dset, ".*(?=_chr_)")) %>% 
  arrange(study)
unique(caqtl$study) # 2 unique caQTL datasets

# Create tissue_type variable based on index position of study variable - tissue_type variable used to create groups for final figure
caqtl = caqtl %>% 
  mutate(tissue_type = case_when(
    study == unique(caqtl$study)[1] ~ "BrainMeta (Multiple regions)",
    study == unique(caqtl$study)[2] ~ "Kosoy 2022 Microglia"))


###############################################################################################################################################
## Combine all Results
annotated = rbind(eqtl, pqtl, sqtl, mqtl, haqtl, caqtl) %>% 
  arrange(desc(best_PP4))

## Run checks - Make sure no missing values in relevant columns
sum(is.na(annotated$best_PP4))
sum(is.na(annotated$tissue_type))
sum(is.na(annotated$qtl_type))
sum(is.na(annotated$study))
summary(annotated$best_PP4)

n_distinct(annotated$locus) # 57
table(annotated$discovery, annotated$locus) 


#####################################################################
## Combine PMFBP1 and HP as one locus
# Extract HP locus from GWAS findings - currently annotated as PMFBP1
hp = annotated %>% 
  filter(locus == "PMFBP1") %>% 
  mutate(locus = "HP")

# Extract all non-HP loci
non_hp = annotated %>% 
  filter(locus != "PMFBP1")

# Make sure some of hp and non_hp equate to full annotated dataframe
nrow(hp) + nrow(non_hp) == nrow(annotated) # TRUE

# rbind hp and non_hp to create new annotated_hp dataframe
annotated_hp = rbind(non_hp, hp)

n_distinct(annotated_hp$locus) # 56 - Drops PMFBP1
table(annotated_hp$discovery, annotated_hp$locus) # GBA, INPP5D, TSPAN14, and HP have GWAS and PWAS results


#####################################################################
## Update WDFY1 locus to exclude sex-het findings
annotated_hp_wdfy1 = annotated_hp %>% 
  filter(ad_sample != "SexHet")


# Write out comprehensive results
fwrite(annotated_hp_wdfy1, file.path(opt$work_dir, "results", opt$out_GP_qtl))

#####################################################################
## Write out separate GWAS and PWAS results seperately - optional
# PWAS
pwas = annotated_hp_wdfy1 %>% 
  filter(discovery != "GWAS") %>% 
  arrange(desc(best_PP4))
fwrite(pwas, file.path(opt$work_dir, "results", opt$out_pwas))

# GWAS
gwas = annotated_hp_wdfy1 %>% 
  filter(discovery == "GWAS") %>% 
  arrange(desc(best_PP4))
fwrite(gwas, file.path(opt$work_dir, "results", opt$out_gwas))

###############################################################################################################################################
### Make filters to m/ca/ha QTL results to exclude low n_snps results and get top hits for each locus/dset/discovery combination

#############################################################################
## mQTL
summary(mqtl$n_snps) # Median = 158

mqtl2 = annotated_hp_wdfy1 %>% 
  filter(qtl_type == "mQTL",
         n_snps >= median(mqtl$n_snps)) %>% 
  group_by(dset, locus) %>% 
  filter(best_PP4 == max(best_PP4)) %>% 
  arrange(CHR, BP)

n_distinct(mqtl2$locus) # 56

# Write out top mQTL results
fwrite(mqtl2, file.path(opt$work_dir, "results", opt$out_mqtl))

#############################################################################
## haQTL
summary(haqtl$n_snps) # Median = 2452

haqtl2 = annotated_hp_wdfy1 %>% 
  filter(qtl_type == "haQTL",
         n_snps >= 2452) %>% 
  group_by(dset, locus) %>% 
  filter(best_PP4 == max(best_PP4)) %>% 
  arrange(CHR, BP)

n_distinct(haqtl2$locus) # 56

# Write out top haQTL results
fwrite(haqtl2, file.path(opt$work_dir, "results", opt$out_haqtl))

#############################################################################
## caQTL
summary(caqtl$n_snps) # Median = 183

caqtl2 = annotated_hp_wdfy1 %>% 
  filter(qtl_type == "caQTL",
         n_snps >= 183) %>% 
  group_by(dset, locus) %>% 
  filter(best_PP4 == max(best_PP4)) %>% 
  arrange(CHR, BP)

n_distinct(caqtl2$locus) # 52 - no caQTL data fro USP4, USP19, MANF, and HLA

# Write out top caQTL results
fwrite(caqtl2, file.path(opt$work_dir, "results", opt$out_caqtl))

###############################################################################################################################################

# End of figure preprocessing script
# Reference figure code for further instructions on creating xQTL overview figure

# sessionInfo()
# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 22631)
# 
# Matrix products: default
# 
# 
# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# time zone: America/Chicago
# tzcode source: internal
# 
# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] scales_1.4.0      gridExtra_2.3     reshape2_1.4.4    patchwork_1.3.1   lubridate_1.9.4  
#  [6] forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.0.4       readr_2.1.5      
# [11] tidyr_1.3.1       tibble_3.2.1      ggplot2_3.5.2     tidyverse_2.0.0   data.table_1.17.0
# 
# loaded via a namespace (and not attached):
#  [1] gtable_0.3.6       compiler_4.4.1     tidyselect_1.2.1   Rcpp_1.0.14        R6_2.6.1          
#  [6] plyr_1.8.9         generics_0.1.4     pillar_1.10.2      RColorBrewer_1.1-3 tzdb_0.5.0        
# [11] rlang_1.1.6        stringi_1.8.7      timechange_0.3.0   cli_3.6.5          withr_3.0.2       
# [16] magrittr_2.0.3     rstudioapi_0.17.1  hms_1.1.3          lifecycle_1.0.4    vctrs_0.6.5       
# [21] glue_1.8.0         farver_2.1.2       tools_4.4.1        pkgconfig_2.0.3   



