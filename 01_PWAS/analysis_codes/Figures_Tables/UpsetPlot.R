#!/usr/bin/env Rscript
# Process PWAS results to make upset plot.
# written by Danielle M Reid

#load library
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(haven)
  library(gridExtra)
  library(ggpattern)
  library(UpSetR)
  library(stringr)
  library(biomaRt)
  library(CMplot)
  library(readxl)
  library(ggpubr)
  library(RColorBrewer)
  library(scales)
  library(patchwork)
  library(forcats)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ReactomePA)
  library(DOSE)
  library(openxlsx)
})

# Define CLI options
option_list <- list(
  # Main directories
  make_option("--out_dir", type = "character", help = "Output directory for final plots and merged files"),

  # File name suffixes (partials)
  make_option("--fb_pf", type = "character", help = "Filename of Female Brain PWAS result (combined weights) with full path"),
  make_option("--fb_if", type = "character", help = "Filename of Female Brain gene novelty info with full pat"),
  make_option("--fc_pf", type = "character", help = "Filename of Female CSF PWAS result (combined weights) with full path"),
  make_option("--fc_if", type = "character", help = "Filename of Female CSF gene novelty info with full path"),
  make_option("--mb_pf", type = "character", help = "Filename of Male Brain PWAS result (combined weights) with full path"),
  make_option("--mb_if", type = "character", help = "Filename of Male Brain gene novelty info with full path"),
  make_option("--mc_pf", type = "character", help = "Filename of Male CSF PWAS result (combined weights) with full path"),
  make_option("--mc_if", type = "character", help = "Filename of Male CSF gene novelty info with full path"),

  # files (for filtering, z-score matching, and top gene selection)
  make_option("--fb", type = "character", help = "Female Brain (sex-stratified) weights file with full path"),
  make_option("--fnb", type = "character", help = "Female Brain (non-stratified) weights file with full path"),
  make_option("--mb", type = "character", help = "Male Brain (sex-stratified) weights file with full path"),
  make_option("--mnb", type = "character", help = "Male Brain (non-stratified) weights file with full path"),
  make_option("--fc", type = "character", help = "Female CSF (sex-stratified) weights file with full path"),
  make_option("--fnc", type = "character", help = "Female CSF (non-stratified) weights file with full path"),
  make_option("--mc", type = "character", help = "Male CSF (sex-stratified) weights file with full path"),
  make_option("--mnc", type = "character", help = "Male CSF (non-stratified) weights file with full path"),

  # Cross-sex weight analysis files
  make_option("--mfb", type = "character", help = "Brain male GWAS with female protein weights PWAS file with full path"),
  make_option("--mfc", type = "character", help = "CFS male GWAS with female protein weights PWAS file with full path"),
  make_option("--fmb", type = "character", help = "Brain female GWAS with male protein weights PWAS file with full path"),
  make_option("--fmc", type = "character", help = "CSF female GWAS with male protein weights PWAS file with full path")
)
# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Assign arguments to variables
dir <- opt$out_dir
#dir = "/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/temp2/"

############################################
# FEMALE
brain_f_PWAS = fread(opt$fb_pf, header = T, sep = '\t')
csf_f_PWAS = fread(opt$fc_pf, header = T, sep = '\t')
#brain_f_PWAS = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/Sex/Females_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-strat_sex-strat_W23_weights.txt", header = T, sep = '\t')
#csf_f_PWAS = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/Sex/Females_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_weights.txt", header = T, sep = '\t')
brain_f_PWAS = as.data.frame(brain_f_PWAS); csf_f_PWAS = as.data.frame(csf_f_PWAS)

brain_f_info = fread(opt$fb_if, header = T, sep = '\t')
csf_f_info = fread(opt$fc_if, header = T, sep = '\t')
#brain_f_info = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/Sex/Females_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-strat_sex-strat_W23_top-female-specific-genes_novelty.txt", header = T, sep = '\t')
#csf_f_info = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/Sex/Females_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_top-female-specific-genes_novelty.txt", header = T, sep = '\t')
brain_f_info = as.data.frame(brain_f_info); csf_f_info = as.data.frame(csf_f_info)

# Initialize NEW_hit column in csf_f_PWAS
csf_f_PWAS$NEW_hit <- NA

# Iterate through csf_f_info and update csf_f_PWAS
for (i in seq_len(nrow(csf_f_info))) {
  match_idx <- which(
    csf_f_PWAS$ID == csf_f_info$ID[i] &
    csf_f_PWAS$Gene == csf_f_info$Gene[i] &
    csf_f_PWAS$NWGT == csf_f_info$NWGT[i] &
    csf_f_PWAS$TSS == csf_f_info$TSS[i] &
    csf_f_PWAS$Discovery == csf_f_info$Discovery[i] &
    csf_f_PWAS$specific == csf_f_info$specific[i])
  
  if (length(match_idx) > 0) {
    csf_f_PWAS$NEW_hit[match_idx] <- csf_f_info$NEW_hit[i]
  }
}

rm(match_idx, i)

# initialize NEW_hit column in brain_f_PWAS
brain_f_PWAS$NEW_hit <- NA

# Iterate through brain_f_info and update brain_f_PWAS
for (i in seq_len(nrow(brain_f_info))) {
  match_idx <- which(
    brain_f_PWAS$Gene == brain_f_info$Gene[i] &
    brain_f_PWAS$NWGT == brain_f_info$NWGT[i] &
    brain_f_PWAS$P0 == brain_f_info$P0[i] &
    brain_f_PWAS$Discovery == brain_f_info$Discovery[i] &
    brain_f_PWAS$specific == brain_f_info$specific[i])
  
  if (length(match_idx) > 0) {
    brain_f_PWAS$NEW_hit[match_idx] <- brain_f_info$NEW_hit[i]
  }
}

# harmonize dataframes so that they can be concatenated
b = subset(brain_f_PWAS, select = c(Gene, CHR, TWAS.Z, TWAS.P, Start_hg38, fdr_p:NEW_hit))
c = subset(csf_f_PWAS, select = c(ID, Gene, CHR, TSS, TWAS.Z, TWAS.P:NEW_hit))

colnames(b) = c("Gene", "CHR", "TWAS.Z", "TWAS.P", "P0", "fdr_p", "specific", "topSNP", "Discovery", "NEW_hit")
b$Analyte = NA
colnames(c) = c("Analyte", "Gene", "CHR", "P0", "TWAS.Z", "TWAS.P", "fdr_p", "specific", "topSNP", "Discovery", "NEW_hit")

b$Tissue = "B"
c$Tissue = "C"

#create new data frame by concatenating brain_f_PWAS and csf_f_PWAS dfs
new_pwas_dat = rbind(b, c)
new_pwas_dat$P0 = as.numeric(new_pwas_dat$P0)

# remove CEBPZOS as sex-specific gene in specific & topSNP
new_pwas_dat$specific[new_pwas_dat$Gene == "CEBPZOS"] <- 0
new_pwas_dat$topSNP[new_pwas_dat$Gene == "CEBPZOS"] <- 0

# save
mov_col = "Analyte"
new_pwas_dat = new_pwas_dat %>% relocate(all_of(mov_col), .before =  "Gene")
fwrite(new_pwas_dat, paste0(dir, "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-female-specific-genes.txt"), na = NA, col.names = T, quote = F, row.names = F, sep = '\t')

rm(list = setdiff(ls(), "opt"))

############################################
# MALES
dir <- opt$out_dir
#dir = "/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/temp2/"

brain_m_PWAS = fread(opt$mb_pf, header = T, sep = '\t')
csf_m_PWAS = fread(opt$mc_pf, header = T, sep = '\t')
#brain_m_PWAS = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/Sex/Males_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-strat_sex-strat_W23_weights.txt", header = T, sep = '\t')
#csf_m_PWAS = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/Sex/Males_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_weights.txt", header = T, sep = '\t')
brain_m_PWAS = as.data.frame(brain_m_PWAS); csf_m_PWAS = as.data.frame(csf_m_PWAS)

brain_m_info = fread(opt$mb_if, header = T, sep = '\t')
csf_m_info = fread(opt$mc_if, header = T, sep = '\t')
#brain_m_info = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/Sex/Males_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-strat_sex-strat_W23_top-male-specific-genes_novelty.txt", header = T, sep = '\t')
#csf_m_info = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/Sex/Males_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_top-male-specific-genes_novelty.txt", header = T, sep = '\t')
brain_m_info = as.data.frame(brain_m_info); csf_m_info = as.data.frame(csf_m_info)

# Initialize NEW_hit column in csf_m_PWAS
csf_m_PWAS$NEW_hit <- NA

# Iterate through csf_m_info and update csf_m_PWAS
for (i in seq_len(nrow(csf_m_info))) {
  match_idx <- which(
    csf_m_PWAS$ID == csf_m_info$ID[i] &
    csf_m_PWAS$Gene == csf_m_info$Gene[i] &
    csf_m_PWAS$NWGT == csf_m_info$NWGT[i] &
    csf_m_PWAS$TSS == csf_m_info$TSS[i] &
    csf_m_PWAS$Discovery == csf_m_info$Discovery[i] &
    csf_m_PWAS$specific == csf_m_info$specific[i])
  
  if (length(match_idx) > 0) {
    csf_m_PWAS$NEW_hit[match_idx] <- csf_m_info$NEW_hit[i]
  }
}

rm(match_idx, i)

# initialize NEW_hit column in brain_m_PWAS
brain_m_PWAS$NEW_hit <- NA

# Iterate through brain_m_info and update brain_m_PWAS
for (i in seq_len(nrow(brain_m_info))) {
  match_idx <- which(
    brain_m_PWAS$Gene == brain_m_info$Gene[i] &
    brain_m_PWAS$NWGT == brain_m_info$NWGT[i] &
    brain_m_PWAS$P0 == brain_m_info$P0[i] &
    brain_m_PWAS$Discovery == brain_m_info$Discovery[i] &
    brain_m_PWAS$specific == brain_m_info$specific[i])
  
  if (length(match_idx) > 0) {
    brain_m_PWAS$NEW_hit[match_idx] <- brain_m_info$NEW_hit[i]
  }
}

# harmonize dataframes so that they can be concatenated
b = subset(brain_m_PWAS, select = c(Gene, CHR, TWAS.Z, TWAS.P, Start_hg38, fdr_p:NEW_hit))
c = subset(csf_m_PWAS, select = c(ID, Gene, CHR, TSS, TWAS.Z, TWAS.P:NEW_hit))

colnames(b) = c("Gene", "CHR", "TWAS.Z", "TWAS.P", "P0", "fdr_p", "specific", "topSNP", "Discovery", "NEW_hit")
b$Analyte = NA
colnames(c) = c("Analyte", "Gene", "CHR", "P0", "TWAS.Z", "TWAS.P", "fdr_p", "specific", "topSNP", "Discovery", "NEW_hit")

b$Tissue = "B"
c$Tissue = "C"

#create new data frame by concatenating brain_m_PWAS and csf_m_PWAS dfs
new_pwas_dat = rbind(b, c)
new_pwas_dat$P0 = as.numeric(new_pwas_dat$P0)

mov_col = "Analyte"
new_pwas_dat = new_pwas_dat %>% relocate(all_of(mov_col), .before =  "Gene")
fwrite(new_pwas_dat, paste0(dir, "ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-male-specific-genes.txt"), na = NA, col.names = T, quote = F, row.names = F, sep = '\t')

rm(list = setdiff(ls(), "opt"))

####################################################################################
#### make upset plot version of sex-specific venn diagrams adding non significant genes (p > 0.05) for each tissue considering collapsing like topSNP
## genes for females in brain p<0.05 collapsing p-values between same tissue discoveries, similar to how i previously identified topSNP genes
## load PWAS result files - genome coordinates are in hg38 and gene names have been updated
dir <- opt$out_dir
#dir = "/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/temp2/"
fb = fread(opt$fb)
fnb = fread(opt$fnb)
mb = fread(opt$mb)
mnb = fread(opt$mnb)
mfb <- fread(opt$mfb)
#fb = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/Sex/Females_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_sex-strat-W23_weights.txt")
#fnb = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/NonSex/Females_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-sex-strat-W23_weights.txt")
#mb = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/Sex/Males_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_sex-strat-W23_weights.txt")
#mnb = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/NonSex/Males_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-sex-strat-W23_weights.txt")
#mfb = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/Sex/MaleGWAS_FemaleProtWGT_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_sex-strat-W23_female-weights.txt")
fb = as.data.frame(fb); fnb = as.data.frame(fnb); mb = as.data.frame(mb); mnb = as.data.frame(mnb); mfb = as.data.frame(mfb)

## sig genes F-F p < 0.05; remove rows where F-F p > 0.05
fb1 = fb[which(fb$TWAS.P < 0.05), ]
# filter 1- not sig in M-M W23 or if sig in M-M W23 PWAS z-score opposite sign
fb1.1 <- fb1[!fb1$ENSG_ID %in% mb[mb$TWAS.P < 0.05, "ENSG_ID"] |
            (fb1$ENSG_ID %in% mb[mb$TWAS.P < 0.05, "ENSG_ID"] & sign(fb1$TWAS.Z) != sign(mb$TWAS.Z[match(fb1$ENSG_ID, mb$ENSG_ID)])),
]
# filter 2- not sig in M-N W23 or if sig in M-N W23 PWAS z-score opposite sign
fb2.1 <- fb1.1[!fb1.1$ENSG_ID %in% mnb[mnb$TWAS.P < 0.05, "ENSG_ID"] |
                (fb1.1$ENSG_ID %in% mnb[mnb$TWAS.P < 0.05, "ENSG_ID"] & sign(fb1.1$TWAS.Z) != sign(mnb$TWAS.Z[match(fb1.1$ENSG_ID, mnb$ENSG_ID)])),
]
#filter 3- not sig in M-F W23 PWAS or if sig in M-F W23 PWAS z-score opposite sign
fb3.1 <- fb2.1[!fb2.1$ENSG_ID %in% mfb[mfb$TWAS.P < 0.05, "ENSG_ID"] |
                (fb2.1$ENSG_ID %in% mfb[mfb$TWAS.P < 0.05, "ENSG_ID"] & sign(fb2.1$TWAS.Z) != sign(mfb$TWAS.Z[match(fb2.1$ENSG_ID, mfb$ENSG_ID)])),
]
## sig genes in F-N p < 0.05; remove rows where f_n p > 0.05
fnb1 = fnb[which(fnb$TWAS.P < 0.05), ]
# filter 1- not sig in M-M W23 or if sig in M-M W23 PWAS z-score opposite sign
fnb1.1 <- fnb1[!fnb1$ENSG_ID %in% mb[mb$TWAS.P < 0.05, "ENSG_ID"] |
            (fnb1$ENSG_ID %in% mb[mb$TWAS.P < 0.05, "ENSG_ID"] & sign(fnb1$TWAS.Z) != sign(mb$TWAS.Z[match(fnb1$ENSG_ID, mb$ENSG_ID)])),
]
# filter 2- not sig in M-N W23 or if sig in M-N W23 PWAS z-score opposite sign
fnb2.1 <- fnb1.1[!fnb1.1$ENSG_ID %in% mnb[mnb$TWAS.P < 0.05, "ENSG_ID"] |
                (fnb1.1$ENSG_ID %in% mnb[mnb$TWAS.P < 0.05, "ENSG_ID"] & sign(fnb1.1$TWAS.Z) != sign(mnb$TWAS.Z[match(fnb1.1$ENSG_ID, mnb$ENSG_ID)])),
]
#filter 3- not sig in M-F W23 PWAS or if sig in M-F W23 PWAS z-score opposite sign
fnb3.1 <- fnb2.1[!fnb2.1$ENSG_ID %in% mfb[mfb$TWAS.P < 0.05, "ENSG_ID"] |
                (fnb2.1$ENSG_ID %in% mfb[mfb$TWAS.P < 0.05, "ENSG_ID"] & sign(fnb2.1$TWAS.Z) != sign(mfb$TWAS.Z[match(fnb2.1$ENSG_ID, mfb$ENSG_ID)])),
]

#make new subset of meta F-F W23 and meta F non-sex-strat W23 data with only necessary columns; add specific column to each with default value of 0 then update to matching ENSG_IDs in corresponding specific gene df
f_dat = as.data.frame(fb)
f_non_dat = as.data.frame(fnb)
f_dat$specific = 0
f_non_dat$specific = 0
f_dat$specific[f_dat$ENSG_ID %in% fb3.1$ENSG_ID & f_dat$TWAS.P %in% fb3.1$TWAS.P] = 1; print(paste("number of f_dat$specific == 1 genes prior to determining common genes = ", sum(f_dat$specific==1)))
f_non_dat$specific[f_non_dat$ENSG_ID %in% fnb3.1$ENSG_ID & f_non_dat$TWAS.P %in% fnb3.1$TWAS.P] = 2; print(paste("number of f_non_dat$specific == 2 genes prior to determining common genes = ", sum(f_non_dat$specific==2)))
#update specific column with a value of 3 for ENSG_IDs that pass the filters in both dataset that match
common_ids = intersect(fb3.1$ENSG_ID, fnb3.1$ENSG_ID); print(paste("number of genes in common between f_dat and f_non_dat = ", length(common_ids)))
f_dat$specific[f_dat$ENSG_ID %in% common_ids] = 3; print(paste("number of genes in f_dat that are shared with f_non_dat = ", sum(f_dat$specific==3)))
f_non_dat$specific[f_non_dat$ENSG_ID %in% common_ids] = 3; print(paste("number of genes in f_non_dat that are shared with f_dat = ", sum(f_non_dat$specific==3)))
#get top SNPs - identify sex-specific genes in both datasets that are unique to one dataset
f_non_only <- setdiff(fnb3.1$ENSG_ID, fb3.1$ENSG_ID); print(paste("number of f_non_only genes = ", length(f_non_only)))
f_sex_only <- setdiff(fb3.1$ENSG_ID, fnb3.1$ENSG_ID); print(paste("number of f_sex_genes only = ", length(f_sex_only)))
#add topSNP column to both datasets meta F-F W23 and meta F non-sex-strat W23, then update topSNP column based on genes unique to a dataset
f_dat$topSNP <- 0
f_non_dat$topSNP <- 0
f_dat$topSNP[f_dat$ENSG_ID %in% f_sex_only] <- 1
f_non_dat$topSNP[f_non_dat$ENSG_ID %in% f_non_only] <- 2
#get top SNPs - identify sex-specific genes that are in both datasets and select the gene with lowest p value for labeling
fs_merge <- merge(fb3.1, fnb3.1, by =  c("ENSG_ID", "ID", "Gene"), suffixes = c("_sex", "_non"))
fs_merge$topSNP <- ifelse(fs_merge$TWAS.P_sex < fs_merge$TWAS.P_non, 1,
                            ifelse(fs_merge$TWAS.P_sex > fs_merge$TWAS.P_non, 2, 3)
)
print(paste("number of genes fs_merge$topSNP == 1 prior to updating f_dat = ", sum(fs_merge$topSNP==1)))
print(paste("number of genes fs_merge$topSNP == 2 prior to updating f_non_dat = ", sum(fs_merge$topSNP==2)))
print(paste("number of genes fs_merge$topSNP == 3 prior to updating f_dat = ", sum(fs_merge$topSNP==3)))
#change value of topSNP in f_dat and f_non_dat for common genes with lowest p value for labeling
#genes that have the same p value in f_sex and f_non that are common genes, for labelling purposes add to topSNP in f_dat
f_dat$topSNP[f_dat$ENSG_ID %in% fs_merge$ENSG_ID & f_dat$ENSG_ID %in% fs_merge$ENSG_ID[fs_merge$topSNP ==1]] <- 3
print(paste("number of genes f_dat$topSNP == 1 after initial update = ", sum(f_dat$topSNP==1))); print(paste("number of genes f_dat$topSNP == 3 after initial update = ", sum(f_dat$topSNP==3)))
f_dat$topSNP[f_dat$ENSG_ID %in% fs_merge$ENSG_ID & f_dat$ENSG_ID %in% fs_merge$ENSG_ID[fs_merge$topSNP ==3]] <- 3
print(paste("final number of genes f_dat$topSNP == 1 after update = ", sum(f_dat$topSNP==1))); print(paste("final number of genes f_dat$topSNP == 3 after update = ", sum(f_dat$topSNP==3)))
f_non_dat$topSNP[f_non_dat$ENSG_ID %in% fs_merge$ENSG_ID & f_non_dat$ENSG_ID %in% fs_merge$ENSG_ID[fs_merge$topSNP ==2]] <- 3
print(paste("final number of genes f_non_dat$topSNP == 2 after update = ", sum(f_non_dat$topSNP==2))); print(paste("final number of genes f_non_dat$topSNP == 3 after update = ", sum(f_non_dat$topSNP==3)))
print(paste("total number of topSNPs = ", sum(sum(f_dat$topSNP==1)+sum(f_dat$topSNP==3)+sum(f_non_dat$topSNP==2)+sum(f_non_dat$topSNP==3))))
#create new data frame by concatenating f_sex and f_non dfs
f_dat$Discovery = "Primary"
f_non_dat$Discovery = "Secondary"
new_pwas_dat = rbind(f_dat, f_non_dat)
new_pwas_dat$P1 = as.numeric(new_pwas_dat$P1)
new_pwas_dat$TWAS.P = as.numeric(new_pwas_dat$TWAS.P)
fwrite(new_pwas_dat, paste0(dir, "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_p0.05_noMHC_ext2Mb_non-strat_sex-strat-W23_weights.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')

rm(list = setdiff(ls(), "opt"))

## genes for females in CSF cis p<0.05 collapsing p-values between same tissue discoveries, similar to how i previously identified topSNP genes
## load PWAS result files
dir <- opt$out_dir
#dir = "/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/temp2/"
fc = fread(opt$fc)
fnc = fread(opt$fnc)
mc = fread(opt$mc)
mnc = fread(opt$mnc)
mfc <- fread(opt$mfc)
#fc = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/Sex/Females_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_sex-strat-CSFcis_weights.txt")
#fnc = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/NonSex/Females_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-sex-strat-CSFcis_weights.txt")
#mc = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/Sex/Males_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_sex-strat-CSFcis_weights.txt")
#mnc = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/NonSex/Males_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_non-sex-strat-CSF_weights-raw.txt")
#mfc = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/Sex/MaleGWAS_FemaleProtWGT_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_sex-strat-CSFcis_female-weights.txt")
fc = as.data.frame(fc); fnc = as.data.frame(fnc); mc = as.data.frame(mc); mnc = as.data.frame(mnc); mfc = as.data.frame(mfc)

## sig genes in F-F CSF cis; remove rows where F-F CSF cis p > 0.05
fc1 = fc[which(fc$TWAS.P < 0.05), ]
# filter 1- not sig in M-M CSF cis or if sig in M-M CSF cis PWAS z-score opposite sign
fc1.1 <- fc1[!fc1$ID %in% mc[mc$TWAS.P < 0.05, "ID"] |
            (fc1$ID %in% mc[mc$TWAS.P < 0.05, "ID"] & sign(fc1$TWAS.Z) != sign(mc$TWAS.Z[match(fc1$ID, mc$ID)])),
]
# filter 2- not sig in M-N CSF cis or if sig in M-N CSF cis PWAS z-score opposite sign
fc2.1 <- fc1.1[!fc1.1$ID %in% mnc[mnc$TWAS.P < 0.05, "ID"] |
                (fc1.1$ID %in% mnc[mnc$TWAS.P < 0.05, "ID"] & sign(fc1.1$TWAS.Z) != sign(mnc$TWAS.Z[match(fc1.1$ID, mnc$ID)])),
]
#filter 3- not sig in M-F CSF cis PWAS or if sig in M-F CSF cis PWAS z-score opposite sign
fc3.1 <- fc2.1[!fc2.1$ID %in% mfc[mfc$TWAS.P < 0.05, "ID"] |
                (fc2.1$ID %in% mfc[mfc$TWAS.P < 0.05, "ID"] & sign(fc2.1$TWAS.Z) != sign(mfc$TWAS.Z[match(fc2.1$ID, mfc$ID)])),
]
## sig genes in F-N p < 0.05; remove rows where f_n p > 0.05
fnc1 = fnc[which(fnc$TWAS.P < 0.05), ]
# filter 1- not sig in M-M CSF cis or if sig in M-M CSF cis PWAS z-score opposite sign
fnc1.1 <- fnc1[!fnc1$ID %in% mc[mc$TWAS.P < 0.05, "ID"] |
            (fnc1$ID %in% mc[mc$TWAS.P < 0.05, "ID"] & sign(fnc1$TWAS.Z) != sign(mc$TWAS.Z[match(fnc1$ID, mc$ID)])),
]
# filter 2- not sig in M-N CSF cis or if sig in M-N CSF cis PWAS z-score opposite sign
fnc2.1 <- fnc1.1[!fnc1.1$ID %in% mnc[mnc$TWAS.P < 0.05, "ID"] |
                (fnc1.1$ID %in% mnc[mnc$TWAS.P < 0.05, "ID"] & sign(fnc1.1$TWAS.Z) != sign(mnc$TWAS.Z[match(fnc1.1$ID, mnc$ID)])),
]
#filter 3- not sig in M-F CSF cis PWAS or if sig in M-F CSF cis PWAS z-score opposite sign
fnc3.1 <- fnc2.1[!fnc2.1$ID %in% mfc[mfc$TWAS.P < 0.05, "ID"] |
                (fnc2.1$ID %in% mfc[mfc$TWAS.P < 0.05, "ID"] & sign(fnc2.1$TWAS.Z) != sign(mfc$TWAS.Z[match(fnc2.1$ID, mfc$ID)])),
]

#make new subset of meta F-F CSF cis and meta F non-sex-strat CSF cis data with only necessary columns; add specific column to each with default value of 0 then update to matching IDs in corresponding specific gene df
f_dat = as.data.frame(fc)
f_non_dat = as.data.frame(fnc)
f_dat$specific = 0
f_non_dat$specific = 0
f_dat$specific[f_dat$ID %in% fc3.1$ID & f_dat$TWAS.P %in% fc3.1$TWAS.P] = 1; print(paste("number of f_dat$specific == 1 genes prior to determining common genes = ", sum(f_dat$specific==1)))
f_non_dat$specific[f_non_dat$ID %in% fnc3.1$ID & f_non_dat$TWAS.P %in% fnc3.1$TWAS.P] = 2; print(paste("number of f_non_dat$specific == 2 genes prior to determining common genes = ", sum(f_non_dat$specific==2)))
#update specific column with a value of 3 for IDs that pass the filters in both dataset that match
common_ids = merge(fc3.1, fnc3.1, by = c("ID", "Gene", "CHR", "TSS")); print(paste("number of genes in common between f_dat and f_non_dat = ", length(common_ids)))
f_dat$specific[f_dat$ID %in% common_ids$ID & f_dat$Gene %in% common_ids$Gene] = 3; print(paste("number of genes in f_dat that are shared with f_non_dat = ", sum(f_dat$specific==3)))
f_non_dat$specific[f_non_dat$ID %in% common_ids$ID & f_non_dat$Gene %in% common_ids$Gene] = 3; print(paste("number of genes in f_non_dat that are shared with f_dat = ", sum(f_non_dat$specific==3)))
#get top SNPs - identify sex-specific genes in both datasets that are unique to one dataset
f_non_only <- f_non_dat[f_non_dat$specific==2,]; print(paste("number of f_non_only genes = ", length(f_non_only)))
f_sex_only <- f_dat[f_dat$specific==1,]; print(paste("number of f_sex_genes only = ", length(f_sex_only)))
#add topSNP column to both datasets meta F-F CSF cis and meta F non-sex-strat CSF cis, then update topSNP column based on genes unique to a dataset
f_dat$topSNP <- 0
f_non_dat$topSNP <- 0
f_dat$topSNP[f_dat$ID %in% f_sex_only$ID & f_dat$Gene %in% f_sex_only$Gene] <- 1
f_non_dat$topSNP[f_non_dat$ID %in% f_non_only$ID & f_non_dat$Gene %in% f_non_only$Gene & f_non_dat$TWAS.P %in% f_non_only$TWAS.P] <- 2
#get top SNPs - identify sex-specific genes that are in both datasets and select the gene with lowest p value for labeling
fs_non_dat <- dplyr::select(fnc3.1, ID, Gene, CHR, TSS, TWAS.P, fdr_p)
fs_dat <- dplyr::select(fc3.1, ID, Gene, CHR, TSS, TWAS.P, fdr_p)
fs_merge <- merge(fc3.1, fnc3.1, by = c("ID", "Gene", "CHR", "TSS"), suffixes = c("_sex", "_non"))
fs_merge$topSNP <- ifelse(fs_merge$TWAS.P_sex < fs_merge$TWAS.P_non, 1,
                            ifelse(fs_merge$TWAS.P_sex > fs_merge$TWAS.P_non, 2, 3)
)
print(paste("number of genes fs_merge$topSNP == 1 prior to updating f_dat = ", sum(fs_merge$topSNP==1)))
print(paste("number of genes fs_merge$topSNP == 2 prior to updating f_non_dat = ", sum(fs_merge$topSNP==2)))
print(paste("number of genes fs_merge$topSNP == 3 prior to updating f_dat = ", sum(fs_merge$topSNP==3)))
#change value of topSNP in f_dat and f_non_dat for common genes with lowest p value for labeling
#genes that have the same p value in f_sex and f_non that are common genes, for labelling purposes add to topSNP in f_dat
f_dat$topSNP[f_dat$TWAS.P %in% fs_merge$TWAS.P_sex & f_dat$ID %in% fs_merge$ID[fs_merge$topSNP ==1]] <- 3
print(paste("number of genes f_dat$topSNP == 1 after initial update = ", sum(f_dat$topSNP==1))); print(paste("number of genes f_dat$topSNP == 3 after initial update = ", sum(f_dat$topSNP==3)))
f_non_dat$topSNP[f_non_dat$TWAS.P %in% fs_merge$TWAS.P_non & f_non_dat$ID %in% fs_merge$ID[fs_merge$topSNP ==2]] <- 3
print(paste("final number of genes f_dat$topSNP == 1 after update = ", sum(f_dat$topSNP==1))); print(paste("final number of genes f_dat$topSNP == 3 after update = ", sum(f_dat$topSNP==3)))
f_dat$topSNP[f_dat$TWAS.P %in% fs_merge$TWAS.P_sex & f_dat$ID %in% fs_merge$ID[fs_merge$topSNP ==3]] <- 3
print(paste("final number of genes f_non_dat$topSNP == 2 after update = ", sum(f_non_dat$topSNP==2))); print(paste("final number of genes f_non_dat$topSNP == 3 after update = ", sum(f_non_dat$topSNP==3)))
print(paste("total number of topSNPs = ", sum(sum(f_dat$topSNP==1)+sum(f_dat$topSNP==3)+sum(f_non_dat$topSNP==2)+sum(f_non_dat$topSNP==3))))
#create new data frame by concatenating f_sex and f_non dfs
f_dat$Discovery = "Primary"
f_non_dat$Discovery = "Secondary"
new_pwas_dat = rbind(f_dat, f_non_dat)
new_pwas_dat$TSS = as.numeric(new_pwas_dat$TSS)
new_pwas_dat$TWAS.P = as.numeric(new_pwas_dat$TWAS.P)
fwrite(new_pwas_dat, paste0(dir, "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_p0.05_noMHC_ext2Mb_non-sex-strat_sex-strat-CSFcis_weights.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')

rm(list = setdiff(ls(), "opt"))

## load combined PWAS results of top sex-specific genes- for set3 and set6
dir <- opt$out_dir
#dir = "/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/temp2/"

f = fread(paste0(dir, "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-female-specific-genes.txt"), sep = '\t')
f = as.data.frame(f)

# load brain PWAS results of genes for females in brain p<0.05 following the 3 sex-specific filters - for set2
fb <- fread(paste0(dir, "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_p0.05_noMHC_ext2Mb_non-strat_sex-strat-W23_weights.txt"))
fb = as.data.frame(fb)

# load CSF cis PWAS results of genes for females in CSF p<0.05 following the 3 sex-specific filters - for set5
fc <- fread(paste0(dir, "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_p0.05_noMHC_ext2Mb_non-sex-strat_sex-strat-CSFcis_weights.txt"))
fc = as.data.frame(fc)

## genes for females in brain fdr_p < 0.05 following the 3 sex-specific filters
set3 <- subset(f, topSNP > 0 & Tissue == "B", select = "Gene")
set3 <- unlist(set3, use.names = FALSE)

# create set2 (genes for females in brain p<0.05 following the 3 sex-specific filters) 
set2 <- subset(fb, topSNP > 0, select = "Gene")
set2 <- unlist(set2, use.names = FALSE)

## genes for females in CSF cis fdr_p < 0.05 following the 3 sex-specific filters
set6 <- subset(f, topSNP > 0 & Tissue == "C", select = "Gene")
set6 <- unlist(set6, use.names = FALSE)

# create set5 (genes for females in CSF cis p<0.05 following the 3 sex-specific filters) 
set5 <- subset(fc, topSNP > 0, select = "Gene")
set5 <- unlist(set5, use.names = FALSE)

# make sure list of genes are unique in each set thus far
set2 <- unique(set2)
set3 <- unique(set3)
set5 <- unique(set5)
set6 <- unique(set6)

## genes for females in brain p > 0.05; remove genes found in the brain sets (set2 and set3) 
x = subset(f, TWAS.P > 0.05 & Tissue == "B", select = "Gene")
x <- x[!x$Gene %in% union(set2, set3), ]
set1 = unique(x)

## genes for females in CSF cis p > 0.05; remove genes found in the CSF cis sets (set5 and set6)
y = subset(f, TWAS.P > 0.05 & Tissue == "C", select = "Gene")
y <- y[!y$Gene %in% union(set5, set6), ]
set4 = unique(y)

set1 <- unique(set1)
set2 <- unique(set2)
set3 <- unique(set3)
set4 <- unique(set4)
set5 <- unique(set5)
set6 <- unique(set6)

setsf <- union(set1, c(set2, set3, set4, set5, set6))
setsf = unique(setsf)

All = data.frame(
    Gene = setsf
)

All$`Brain_P>0.05` = 0
All$`Brain_P<0.05` = 0
All$`Brain_P.fdr<0.05` = 0
All$`CSF_P.fdr<0.05` = 0
All$`CSF_P<0.05` = 0
All$`CSF_P>0.05` = 0

All$`Brain_P>0.05`[All$Gene %in% set1] = 1
All$`Brain_P<0.05`[All$Gene %in% set2] = 1
All$`Brain_P.fdr<0.05`[All$Gene %in% set3] = 1
All$`CSF_P.fdr<0.05`[All$Gene %in% set6] = 1
All$`CSF_P<0.05`[All$Gene %in% set5] = 1
All$`CSF_P>0.05`[All$Gene %in% set4] = 1

All$Brain_CSF_P <- 0

# get z scores for fb and fc where p < 0.05 to then check if replicated genes across tissue have concordant z-score sign across tissue
fb1 <- subset(fb, topSNP > 0, select = c("Gene", "TWAS.Z"))
fb1$Zsign = sign(fb1$TWAS.Z)
fb1$Count = 1
fb1 = dplyr::select(fb1, -TWAS.Z)
fb1 = unique(fb1)

fc1 <- subset(fc, topSNP > 0, select = c("Gene", "TWAS.Z"))
fc1$Zsign = sign(fc1$TWAS.Z)
fc1$Count = 1
fc1 = dplyr::select(fc1, -TWAS.Z)
fc1 = unique(fc1)

combined_zsign <- merge(fb1, fc1, by = "Gene", suffixes = c(".fb", ".fc"))
valid_genes <- combined_zsign$Gene[combined_zsign$Zsign.fb == combined_zsign$Zsign.fc]

# Update the Brain_CSF_P column based on combined conditions
All$Brain_CSF_P <- ifelse(
  All$`Brain_P<0.05` == 1 & All$`CSF_P<0.05` == 1 & All$Gene %in% valid_genes, 1, 0
)

fwrite(All, paste0(dir, "Overlap_W23_NGI-CSF_female_fdr-p3.txt"), col.names = T, row.names=F,quote=F, sep='\t')

colnames(All) = c("Gene", "Brain_P>0.05", "Brain_P<0.05", "Brain_P.fdr<0.05_&_Sex-specific", "CSF_P.fdr<0.05_&_Sex-specific", "CSF_P<0.05", "CSF_P>0.05", "Brain_CSF_P")

pdf(paste0(dir, "Overlap_W23_NGI-CSF_female_fdr-p3.pdf"), height = 5, onefile = FALSE)
upset(All, sets = c("CSF_P.fdr<0.05_&_Sex-specific", "CSF_P<0.05", "CSF_P>0.05", "Brain_P.fdr<0.05_&_Sex-specific", "Brain_P<0.05", "Brain_P>0.05"),
      keep.order = T, order.by = "freq", mainbar.y.label = "CSF and Brain Gene/Protein Replication",
      sets.x.label = "Proteins/Genes", text.scale = c(1.4, 1.4, 1.4, 1.4, 1.5, 1.4)
)
dev.off()
# set size axis tick labels are too packed and you cannot customize tick labeling other than size
# and you cannot add space to the left margin to allow more appropriate spacing between ticks

#vague set names so that the histogram is wider
# Brain1 = "Brain_P>0.05"
# Brain2 = "Brain_P<0.05"
# Brain3 = "Brain_P.fdr<0.05_&_Sex-specific"
# CSF3 = "CSF_P.fdr<0.05_&_Sex-specific"
# CSF2 = "CSF_P<0.05"
# CSF1 = "CSF_P>0.05"
colnames(All) = c("Gene", "Brain1", "Brain2", "Brain3", "CSF3", "CSF2", "CSF1", "Brain_CSF_P")

pdf(paste0(dir, "Overlap_W23_NGI-CSF_female_fdr-p3.v2.pdf"), height = 5, onefile = FALSE)
upset(All, sets = c("CSF3", "CSF2", "CSF1", "Brain3", "Brain2", "Brain1"),
      keep.order = T, order.by = "freq", mainbar.y.label = "CSF and Brain Gene/Protein Replication",
      sets.x.label = "Proteins/Genes", text.scale = c(1.4, 1.4, 1.4, 1.4, 1.5, 1.4)
)
dev.off()

rm(list = setdiff(ls(), "opt"))

#####################
#### MALES
## genes for males in brain p<0.05 collapsing p-values between same tissue discoveries, similar to how i previously identified topSNP genes
## load PWAS result files - genome coordinates are in hg38 and gene names have been updated
dir <- opt$out_dir
#dir = "/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/temp2/"
mb = fread(opt$mb)
mnb = fread(opt$mnb)
fb = fread(opt$fb)
fnb = fread(opt$fnb)
fmb <- fread(opt$fmb)
#mb = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/Sex/Males_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_sex-strat-W23_weights.txt")
#mnb = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/NonSex/Males_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-sex-strat-W23_weights.txt")
#fb = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/Sex/Females_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_sex-strat-W23_weights.txt")
#fnb = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/NonSex/Females_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-sex-strat-W23_weights.txt")
#fmb = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Brain/Sex/FemaleGWAS_MaleProtWGT_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_sex-strat-W23_male-weights.txt")
fb = as.data.frame(fb); fnb = as.data.frame(fnb); mb = as.data.frame(mb); mnb = as.data.frame(mnb); fmb = as.data.frame(fmb)

##sig genes M-M; remove rows where M-M fdr_p > 0.05
mb1 = mb[which(mb$TWAS.P < 0.05), ]
#filter 1- not sig in F-F W23 PWAS or if sig in F-F W23 PWAS z-score opposite sign
mb1.1 <- mb1[!mb1$ENSG_ID %in% fb[fb$TWAS.P < 0.05, "ENSG_ID"] |
             (mb1$ENSG_ID %in% fb[fb$TWAS.P < 0.05, "ENSG_ID"]& sign(mb1$TWAS.Z) != sign(fb$TWAS.Z[match(mb1$ENSG_ID, fb$ENSG_ID)])), 
]
#filter 2- not sig in F-N W23 PWAS or if sig in F-N W23 PWAS z-score opposite
mb2.1 <- mb1.1[!mb1.1$ENSG_ID %in% fnb[fnb$TWAS.P < 0.05, "ENSG_ID"] | 
             (mb1.1$ENSG_ID %in% fnb[fnb$TWAS.P < 0.05, "ENSG_ID"] & sign(mb1.1$TWAS.Z) != sign(fnb$TWAS.Z[match(mb1.1$ENSG_ID, fnb$ENSG_ID)])), 
]
#filter 3- not sig in F-M W23 PWAS or if sig in F-M W23 PWAS z-score opposite sign
mb3.1 <- mb2.1[!mb2.1$ENSG_ID %in% fmb[fmb$TWAS.P < 0.05, "ENSG_ID"] | 
             (mb2.1$ENSG_ID %in% fmb[fmb$TWAS.P < 0.05, "ENSG_ID"] & sign(mb2.1$TWAS.Z) != sign(fmb$TWAS.Z[match(mb2.1$ENSG_ID, fmb$ENSG_ID)])), 
]
## sig genes in M-N p < 0.05; remove rows where m_n p > 0.05
mnb1 = mnb[which(mnb$fdr_p<0.05),]
#filter 1- not sig in F-F W23 PWAS or if sig has opposite z-score
mnb1.1 <- mnb1[!mnb1$ENSG_ID %in% fb[fb$TWAS.P < 0.05, "ENSG_ID"] | 
             (mnb1$ENSG_ID %in% fb[fb$TWAS.P < 0.05, "ENSG_ID"] & sign(mnb1$TWAS.Z) != sign(fb$TWAS.Z[match(mnb1$ENSG_ID, fb$ENSG_ID)])), 
]
#filter 2- not sig in F non-sex-strat W23 PWAS or if sig z score opposite sign
mnb2.1 <- mnb1.1[!mnb1.1$ENSG_ID %in% fnb[fnb$TWAS.P<0.05,"ENSG_ID"] |
              (mnb1.1$ENSG_ID %in% fnb[fnb$TWAS.P < 0.05,"ENSG_ID"] & sign(mnb1.1$TWAS.Z) != sign(fnb$TWAS.Z[match(mnb1.1$ENSG_ID, fnb$ENSG_ID)])),
]
#filter 3- not sig in F-M W23 PWAS or if sig opposite z-score
mnb3.1 <- mnb2.1[!mnb2.1$ENSG_ID %in% fmb[fmb$TWAS.P<0.05,"ENSG_ID"] |
              (mnb2.1$ENSG_ID %in% fmb[fmb$TWAS.P < 0.05,"ENSG_ID"] & sign(mnb2.1$TWAS.Z) != sign(fmb$TWAS.Z[match(mnb2.1$ENSG_ID, fmb$ENSG_ID)])),
]

#make new subset of meta M-M W23 and meta M non-sex-strat W23 data with only necessary columns; add specific column to each with default value of 0 then update to matching ENSG_IDs in corresponding specific gene df
m_dat = as.data.frame(mb)
m_non_dat = as.data.frame(mnb)
m_dat$specific = 0
m_non_dat$specific = 0
m_dat$specific[m_dat$ENSG_ID %in% mb3.1$ENSG_ID & m_dat$TWAS.P %in% mb3.1$TWAS.P] = 1; print(paste("number of m_dat$specific == 1 genes prior to determining common genes = ", sum(m_dat$specific==1)))
m_non_dat$specific[m_non_dat$ENSG_ID %in% mnb3.1$ENSG_ID & m_non_dat$TWAS.P %in% mnb3.1$TWAS.P] = 2; print(paste("number of m_non_dat$specific == 2 genes prior to determining common genes = ", sum(m_non_dat$specific==2)))
#update specific column with a value of 3 for ENSG_IDs that pass the filters in both dataset that match
common_ids = intersect(mb3.1$ENSG_ID, mnb3.1$ENSG_ID); print(paste("number of genes in common between m_dat and m_non_dat = ", length(common_ids)))
m_dat$specific[m_dat$ENSG_ID %in% common_ids] = 3; print(paste("number of genes in m_dat that are shared with m_non_dat = ", sum(m_dat$specific==3)))
m_non_dat$specific[m_non_dat$ENSG_ID %in% common_ids] = 3; print(paste("number of genes in m_non_dat that are shared with m_dat = ", sum(m_non_dat$specific==3)))
#get top SNPs - identify sex-specific genes in both datasets that are unique to one dataset
m_non_only <- setdiff(mnb3.1$ENSG_ID, mb3.1$ENSG_ID); print(paste("number of m_non_only genes = ", length(m_non_only)))
m_sex_only <- setdiff(mb3.1$ENSG_ID, mnb3.1$ENSG_ID); print(paste("number of m_sex_genes only = ", length(m_sex_only)))
#add topSNP column to both datasets meta M-M W23 and meta M non-sex-strat W23, then update topSNP column based on genes unique to a dataset
m_dat$topSNP <- 0
m_non_dat$topSNP <- 0
m_dat$topSNP[m_dat$ENSG_ID %in% m_sex_only] <- 1
m_non_dat$topSNP[m_non_dat$ENSG_ID %in% m_non_only] <- 2
#get top SNPs - identify sex-specific genes that are in both datasets and select the gene with lowest p value for labeling
ms_merge <- merge(mb3.1, mnb3.1, by =  c("ENSG_ID", "ID", "Gene"), suffixes = c("_sex", "_non"))
ms_merge$topSNP <- ifelse(ms_merge$TWAS.P_sex < ms_merge$TWAS.P_non, 1,
                            ifelse(ms_merge$TWAS.P_sex > ms_merge$TWAS.P_non, 2, 3)
)
print(paste("number of genes ms_merge$topSNP == 1 prior to updating m_dat = ", sum(ms_merge$topSNP==1)))
print(paste("number of genes ms_merge$topSNP == 2 prior to updating m_non_dat = ", sum(ms_merge$topSNP==2)))
print(paste("number of genes ms_merge$topSNP == 3 prior to updating m_dat = ", sum(ms_merge$topSNP==3)))
#change value of topSNP in m_dat and m_non_dat for common genes with lowest p value for labeling
#genes that have the same p value in m_sex and m_non that are common genes, for labelling purposes add to topSNP in m_dat
m_dat$topSNP[m_dat$ENSG_ID %in% ms_merge$ENSG_ID & m_dat$ENSG_ID %in% ms_merge$ENSG_ID[ms_merge$topSNP ==1]] <- 3
print(paste("number of genes m_dat$topSNP == 1 after initial update = ", sum(m_dat$topSNP==1))); print(paste("number of genes m_dat$topSNP == 3 after initial update = ", sum(m_dat$topSNP==3)))
m_dat$topSNP[m_dat$ENSG_ID %in% ms_merge$ENSG_ID & m_dat$ENSG_ID %in% ms_merge$ENSG_ID[ms_merge$topSNP ==3]] <- 3
print(paste("final number of genes m_dat$topSNP == 1 after update = ", sum(m_dat$topSNP==1))); print(paste("final number of genes m_dat$topSNP == 3 after update = ", sum(m_dat$topSNP==3)))
m_non_dat$topSNP[m_non_dat$ENSG_ID %in% ms_merge$ENSG_ID & m_non_dat$ENSG_ID %in% ms_merge$ENSG_ID[ms_merge$topSNP ==2]] <- 3
print(paste("final number of genes m_non_dat$topSNP == 2 after update = ", sum(m_non_dat$topSNP==2))); print(paste("final number of genes m_non_dat$topSNP == 3 after update = ", sum(m_non_dat$topSNP==3)))
print(paste("total number of topSNPs = ", sum(sum(m_dat$topSNP==1)+sum(m_dat$topSNP==3)+sum(m_non_dat$topSNP==2)+sum(m_non_dat$topSNP==3))))
#create new data frame by concatenating m_sex and m_non dfs
m_dat$Discovery = "Primary"
m_non_dat$Discovery = "Secondary"
new_pwas_dat = rbind(m_dat, m_non_dat)
new_pwas_dat$P1 = as.numeric(new_pwas_dat$P1)
new_pwas_dat$TWAS.P = as.numeric(new_pwas_dat$TWAS.P)
fwrite(new_pwas_dat, paste0(dir, "ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_p0.05_noMHC_ext2Mb_non-strat_sex-strat-W23_weights.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')

rm(list = setdiff(ls(), "opt"))

## genes for males in CSF cis p<0.05 collapsing p-values between same tissue discoveries, similar to how i previously identified topSNP genes
## load PWAS result files
dir <- opt$out_dir
#dir = "/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/temp2/"
mc = fread(opt$mc)
mnc = fread(opt$mnc)
fc = fread(opt$fc)
fnc = fread(opt$fnc)
fmc <- fread(opt$fmc)
#mc = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/Sex/Males_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_sex-strat-CSFcis_weights.txt")
#mnc = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/NonSex/Males_results/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-sex-strat-CSFcis_weights.txt")
#fc = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/Sex/Females_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_sex-strat-CSFcis_weights.txt")
#fnc = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/NonSex/Females_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-sex-strat-CSFcis_weights.txt")
#fmc = fread("/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/CSF/Sex/FemaleGWAS_MaleProtWGT_results/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_sex-strat-CSFcis_male-weights.txt")
fc = as.data.frame(fc); fnc = as.data.frame(fnc); mc = as.data.frame(mc); mnc = as.data.frame(mnc); fmc = as.data.frame(fmc)

#sig genes M-M; remove rows where M-M fdr_p > 0.05
mc1 = mc[which(mc$TWAS.P<0.05),]
#filter 1- not sig in F-F CSF cis PWAS or if sig in F-F CSF cis PWAS z-score opposite sign
mc1.1 <- mc1[!mc1$ID %in% fc[fc$TWAS.P < 0.05, "ID"] |
            (mc1$ID %in% fc[fc$TWAS.P < 0.05, "ID"] & sign(mc1$TWAS.Z) != sign(fc$TWAS.Z[match(mc1$ID, fc$ID)])),
]
#filter 2- not sig in F-N CSF PWAS or if sig in F-N CSF PWAS z-score opposite
mc2.1 <- mc1.1[!mc1.1$ID %in% fnc[fnc$TWAS.P < 0.05, "ID"] | 
             (mc1.1$ID %in% fnc[fnc$TWAS.P < 0.05, "ID"] & sign(mc1.1$TWAS.Z) != sign(fnc$TWAS.Z[match(mc1.1$ID, fnc$ID)])), 
]
#filter 3- not sig in F-M CSF cis PWAS or if sig in F-M CSF cis PWAS z-score opposite sign
mc3.1 <- mc2.1[!mc2.1$ID %in% fmc[fmc$TWAS.P < 0.05, "ID"] | 
             (mc2.1$ID %in% fmc[fmc$TWAS.P < 0.05, "ID"] & sign(mc2.1$TWAS.Z) != sign(fmc$TWAS.Z[match(mc2.1$ID, fmc$ID)])), 
]
#sig for males non-sex-strat CSF cis; removes rows where fdr_p >0.05
mnc1 = mnc[which(mnc$TWAS.P<0.05),]
#filter 1- not sig in F-F CSF cis PWAS or if sig has opposite z-score
mnc1.1 <- mnc1[!mnc1$ID %in% fc[fc$TWAS.P < 0.05, "ID"] |
            (mnc1$ID %in% fc[fc$TWAS.P < 0.05, "ID"] & sign(mnc1$TWAS.Z) != sign(fc$TWAS.Z[match(mnc1$ID, fc$ID)])),
]
#filter 2- not sig in F non-sex-strat CSF cis PWAS or if sig z score opposite sign
mnc2.1 <- mnc1.1[!mnc1.1$ID %in% fnc[fnc$TWAS.P < 0.05, "ID"] |
                (mnc1.1$ID %in% fnc[fnc$TWAS.P < 0.05, "ID"] & sign(mnc1.1$TWAS.Z) != sign(fnc$TWAS.Z[match(mnc1.1$ID, fnc$ID)])),
]
#filter 3- not sig in F-M CSF cis PWAS or if sig opposite z-score
mnc3.1 <- mnc2.1[!mnc2.1$ID %in% fmc[fmc$TWAS.P < 0.05, "ID"] |
                (mnc2.1$ID %in% fmc[fmc$TWAS.P < 0.05, "ID"] & sign(mnc2.1$TWAS.Z) != sign(fmc$TWAS.Z[match(mnc2.1$ID, fmc$ID)])),
]

#make new subset of meta M-M CSF cis and meta M non-sex-strat CSF cis data with only necessary columns; add specific column to each with default value of 0 then update to matching IDs in corresponding specific gene df
m_dat = as.data.frame(mc)
m_non_dat = as.data.frame(mnc)
m_dat$specific = 0
m_non_dat$specific = 0
m_dat$specific[m_dat$ID %in% mc3.1$ID & m_dat$TWAS.P %in% mc3.1$TWAS.P] = 1; print(paste("number of m_dat$specific == 1 genes prior to determining common genes = ", sum(m_dat$specific==1)))
m_non_dat$specific[m_non_dat$ID %in% mnc3.1$ID & m_non_dat$TWAS.P %in% mnc3.1$TWAS.P] = 2; print(paste("number of m_non_dat$specific == 2 genes prior to determining common genes = ", sum(m_non_dat$specific==2)))
#update specific column with a value of 3 for IDs that pass the filters in both dataset that match
common_ids = merge(mc3.1, mnc3.1, by = c("ID", "Gene", "CHR", "TSS")); print(paste("number of genes in common between m_dat and m_non_dat = ", length(common_ids)))
m_dat$specific[m_dat$ID %in% common_ids$ID & m_dat$Gene %in% common_ids$Gene] = 3; print(paste("number of genes in m_dat that are shared with m_non_dat = ", sum(m_dat$specific==3)))
m_non_dat$specific[m_non_dat$ID %in% common_ids$ID & m_non_dat$Gene %in% common_ids$Gene] = 3; print(paste("number of genes in m_non_dat that are shared with m_dat = ", sum(m_non_dat$specific==3)))
#get top SNPs - identify sex-specific genes in both datasets that are unique to one dataset
m_non_only <- m_non_dat[m_non_dat$specific==2,]; print(paste("number of m_non_only genes = ", length(m_non_only)))
m_sex_only <- m_dat[m_dat$specific==1,]; print(paste("number of m_sex_genes only = ", length(m_sex_only)))
#add topSNP column to both datasets meta M-M CSF cis and meta M non-sex-strat CSF cis, then update topSNP column based on genes unique to a dataset
m_dat$topSNP <- 0
m_non_dat$topSNP <- 0
m_dat$topSNP[m_dat$ID %in% m_sex_only$ID & m_dat$Gene %in% m_sex_only$Gene] <- 1
m_non_dat$topSNP[m_non_dat$ID %in% m_non_only$ID & m_non_dat$Gene %in% m_non_only$Gene & m_non_dat$TWAS.P %in% m_non_only$TWAS.P] <- 2
#get top SNPs - identify sex-specific genes that are in both datasets and select the gene with lowest p value for labeling
ms_merge <- merge(mc3.1, mnc3.1, by = c("ID", "Gene", "CHR", "TSS"), suffixes = c("_sex", "_non"))
ms_merge$topSNP <- ifelse(ms_merge$TWAS.P_sex < ms_merge$TWAS.P_non, 1,
                            ifelse(ms_merge$TWAS.P_sex > ms_merge$TWAS.P_non, 2, 3)
)
print(paste("number of genes ms_merge$topSNP == 1 prior to updating m_dat = ", sum(ms_merge$topSNP==1)))
print(paste("number of genes ms_merge$topSNP == 2 prior to updating m_non_dat = ", sum(ms_merge$topSNP==2)))
print(paste("number of genes ms_merge$topSNP == 3 prior to updating m_dat = ", sum(ms_merge$topSNP==3)))
#change value of topSNP in m_dat and m_non_dat for common genes with lowest p value for labeling
#genes that have the same p value in m_sex and m_non that are common genes, for labelling purposes add to topSNP in m_dat
m_dat$topSNP[m_dat$TWAS.P %in% ms_merge$TWAS.P_sex & m_dat$ID %in% ms_merge$ID[ms_merge$topSNP ==1]] <- 3
print(paste("number of genes m_dat$topSNP == 1 after initial update = ", sum(m_dat$topSNP==1))); print(paste("number of genes m_dat$topSNP == 3 after initial update = ", sum(m_dat$topSNP==3)))
m_non_dat$topSNP[m_non_dat$TWAS.P %in% ms_merge$TWAS.P_non & m_non_dat$ID %in% ms_merge$ID[ms_merge$topSNP ==2]] <- 3
print(paste("final number of genes m_dat$topSNP == 1 after update = ", sum(m_dat$topSNP==1))); print(paste("final number of genes m_dat$topSNP == 3 after update = ", sum(m_dat$topSNP==3)))
m_dat$topSNP[m_dat$TWAS.P %in% ms_merge$TWAS.P_sex & m_dat$ID %in% ms_merge$ID[ms_merge$topSNP ==3]] <- 3
print(paste("final number of genes m_non_dat$topSNP == 2 after update = ", sum(m_non_dat$topSNP==2))); print(paste("final number of genes m_non_dat$topSNP == 3 after update = ", sum(m_non_dat$topSNP==3)))
print(paste("total number of topSNPs = ", sum(sum(m_dat$topSNP==1)+sum(m_dat$topSNP==3)+sum(m_non_dat$topSNP==2)+sum(m_non_dat$topSNP==3))))
#create new data frame by concatenating m_sex and m_non dfs
m_dat$Discovery = "Primary"
m_non_dat$Discovery = "Secondary"
new_pwas_dat = rbind(m_dat, m_non_dat)
new_pwas_dat$TSS = as.numeric(new_pwas_dat$TSS)
new_pwas_dat$TWAS.P = as.numeric(new_pwas_dat$TWAS.P)
fwrite(new_pwas_dat, paste0(dir, "ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_p0.05_noMHC_ext2Mb_non-sex-strat_sex-strat-CSFcis_weights.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')

rm(list = setdiff(ls(), "opt"))

## load combined PWAS results of top sex-specific genes- for set3 and set6
dir <- opt$out_dir
#dir = "/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/temp2/"
m = fread(paste0(dir, "ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-male-specific-genes.txt"), sep = '\t')
m = as.data.frame(m)

# load brain PWAS results of genes for males in brain p<0.05 following the 3 sex-specific filters - for set2
mb <- fread(paste0(dir, "ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_p0.05_noMHC_ext2Mb_non-strat_sex-strat-W23_weights.txt"))
mb = as.data.frame(mb)

# load CSF cis PWAS results of genes for males in brain p<0.05 following the 3 sex-specific filters - for set5
mc <- fread(paste0(dir, "ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_p0.05_noMHC_ext2Mb_non-sex-strat_sex-strat-CSFcis_weights.txt"))
mc = as.data.frame(mc)

## genes for males in brain fdr_p < 0.05 following the 3 sex-specific filters
set3 <- subset(m, topSNP > 0 & Tissue == "B", select = "Gene")
set3 <- unlist(set3, use.names = FALSE)

# create set2 (genes for males in brain p<0.05 following the 3 sex-specific filters) 
set2 <- subset(mb, topSNP > 0, select = "Gene")
set2 <- unlist(set2, use.names = FALSE)

## genes for males in CSF cis fdr_p < 0.05 following the 3 sex-specific filters
set6 <- subset(m, topSNP > 0 & Tissue == "C", select = "Gene")
set6 <- unlist(set6, use.names = FALSE)

# create set5 (genes for males in CSF cis p<0.05 following the 3 sex-specific filters) 
set5 <- subset(mc, topSNP > 0, select = "Gene")
set5 <- unlist(set5, use.names = FALSE)

# make sure list of genes are unique in each set thus far
set2 <- unique(set2)
set3 <- unique(set3)
set5 <- unique(set5)
set6 <- unique(set6)

## genes for males in brain p > 0.05; remove genes found in the brain sets (set2 and set3) 
x = subset(m, TWAS.P > 0.05 & Tissue == "B", select = "Gene")
x <- x[!x$Gene %in% union(set2, set3), ]
set1 = unique(x)

## genes for males in CSF cis p > 0.05; remove genes found in the CSF cis sets (set5 and set6)
y = subset(m, TWAS.P > 0.05 & Tissue == "C", select = "Gene")
y <- y[!y$Gene %in% union(set5, set6), ]
set4 = unique(y)

set1 <- unique(set1)
set2 <- unique(set2)
set3 <- unique(set3)
set4 <- unique(set4)
set5 <- unique(set5)
set6 <- unique(set6)

setsm <- union(set1, c(set2, set3, set4, set5, set6))
setsm = unique(setsm)

All = data.frame(
    Gene = setsm
)

All$`Brain_P>0.05` = 0
All$`Brain_P<0.05` = 0
All$`Brain_P.fdr<0.05` = 0
All$`CSF_P.fdr<0.05` = 0
All$`CSF_P<0.05` = 0
All$`CSF_P>0.05` = 0

All$`Brain_P>0.05`[All$Gene %in% set1] = 1
All$`Brain_P<0.05`[All$Gene %in% set2] = 1
All$`Brain_P.fdr<0.05`[All$Gene %in% set3] = 1
All$`CSF_P.fdr<0.05`[All$Gene %in% set6] = 1
All$`CSF_P<0.05`[All$Gene %in% set5] = 1
All$`CSF_P>0.05`[All$Gene %in% set4] = 1

All$Brain_CSF_P <- 0

# get z scores for mb and mc where p < 0.05 to then check if replicated genes across tissue have concordant z-score sign across tissue
mb1 <- subset(mb, topSNP > 0, select = c("Gene", "TWAS.Z"))
mb1$Zsign = sign(mb1$TWAS.Z)
mb1$Count = 1
mb1 = dplyr::select(mb1, -TWAS.Z)
mb1 = unique(mb1)

mc1 <- subset(mc, topSNP > 0, select = c("Gene", "TWAS.Z"))
mc1$Zsign = sign(mc1$TWAS.Z)
mc1$Count = 1
mc1 = dplyr::select(mc1, -TWAS.Z)
mc1 = unique(mc1)

combined_zsign <- merge(mb1, mc1, by = "Gene", suffixes = c(".mb", ".mc"))
valid_genes <- combined_zsign$Gene[combined_zsign$Zsign.mb == combined_zsign$Zsign.mc]

# Update the Brain_CSF_P column based on combined conditions
All$Brain_CSF_P <- ifelse(
  All$`Brain_P<0.05` == 1 & All$`CSF_P<0.05` == 1 & All$Gene %in% valid_genes, 1, 0
)

fwrite(All, paste0(dir, "Overlap_W23_NGI-CSF_male_fdr-p3.txt"), col.names = T, row.names=F,quote=F, sep='\t')

colnames(All) = c("Gene", "Brain_P>0.05", "Brain_P<0.05", "Brain_P.fdr<0.05_&_Sex-specific", "CSF_P.fdr<0.05_&_Sex-specific", "CSF_P<0.05", "CSF_P>0.05", "Brain_CSF_P")

pdf(paste0(dir, "Overlap_W23_NGI-CSF_male_fdr-p3.pdf"), height = 5, onefile = FALSE)
upset(All, sets = c("CSF_P.fdr<0.05_&_Sex-specific", "CSF_P<0.05", "CSF_P>0.05", "Brain_P.fdr<0.05_&_Sex-specific", "Brain_P<0.05", "Brain_P>0.05"),
      keep.order = T, order.by = "freq", mainbar.y.label = "CSF and Brain Gene/Protein Replication",
      sets.x.label = "Proteins/Genes", text.scale = c(1.4, 1.4, 1.4, 1.4, 1.5, 1.4)
)
dev.off()

colnames(All) = c("Gene", "Brain1", "Brain2", "Brain3", "CSF3", "CSF2", "CSF1", "Brain_CSF_P")

pdf(paste0(dir, "Overlap_W23_NGI-CSF_male_fdr-p3.v2.pdf"), height = 5, onefile = FALSE)
upset(All, sets = c("CSF3", "CSF2", "CSF1", "Brain3", "Brain2", "Brain1"),
      keep.order = T, order.by = "freq", mainbar.y.label = "CSF and Brain Gene/Protein Replication",
      sets.x.label = "Proteins/Genes", text.scale = c(1.4, 1.4, 1.4, 1.4, 1.5, 1.4)
)
dev.off()
