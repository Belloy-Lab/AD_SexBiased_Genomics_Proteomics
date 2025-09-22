#### make sex-specific gene lists
###################### Danielle Reid
#################################################################################################################

#load the packages
library(plyr)
library(dplyr)
library(readr)
library(data.table)
library(tidyr)
library(methods)
library(optparse)
library(stringr)
library(biomaRt)

#############################################################
option_list = list(
  make_option("--brain_dir", action="store", default=NA, type='character',
              help="Path to Brain primary discovery base directory [required]"),
  make_option("--csf_dir", action="store", default=NA, type='character',
              help="Path to CSF primary discovery base directory [required]"),
  make_option("--out_dir", action="store", default=NA, type='character',
              help="Path to stora sex specific files [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

#############################################################

# create variables for base directory
## create variables for directories and file base names
# the brain PWAS coordinates have already been lifted to hg38 and gene names have been updated across tissue
brain_dir = opt$brain_dir; print(paste("Brain primary discovery directory:", brain_dir)) # /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/Brain/Sex/"
csf_dir = opt$csf_dir; print(paste("CSF primary discovery directory:", csf_dir)) # /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/CSF/Sex/"
dir = opt$csf_dir; print(paste("CSF primary discovery directory:", csf_dir)) # "/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/Final_plots/"

f_dir = "Female/"
m_dir = "Male/"
f_Bf = "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS"
m_Bf = "ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS"
f_Cf = "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP"
m_Cf = "ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_wHP"

brain_PWAS = fread(paste0(brain_dir, f_dir, f_Bf, "_noMHC_ext2Mb_non-strat_sex-strat_W23_weights.txt"), header = T, sep = '\t')
csf_PWAS = fread(paste0(csf_dir, f_dir, f_Cf, "_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_weights.txt"), header = T, sep = '\t')
brain_PWAS = as.data.frame(brain_PWAS); csf_PWAS = as.data.frame(csf_PWAS)

## prepping to plot combined brain and CSF sex-specific results

# add novelty info
csf_info = fread(paste0(csf_dir, f_dir, f_Cf, "_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_top-female-specific-genes_novelty.txt"), header = T, sep = '\t')
brain_info = fread(paste0(brain_dir, f_dir, f_Bf, "_noMHC_ext2Mb_non-strat_sex-strat_W23_top-female-specific-genes_novelty.txt"), header = T, sep = '\t')
csf_info = as.data.frame(csf_info); brain_info = as.data.frame(brain_info)

# Initialize NEW_hit column in csf_PWAS
csf_PWAS$NEW_hit <- NA

# Iterate through csf_info and update csf_PWAS
for (i in seq_len(nrow(csf_info))) {
  match_idx <- which(
    csf_PWAS$ID == csf_info$ID[i] &
    csf_PWAS$Gene == csf_info$Gene[i] &
    csf_PWAS$NWGT == csf_info$NWGT[i] &
    csf_PWAS$TSS == csf_info$TSS[i] &
    csf_PWAS$Discovery == csf_info$Discovery[i] &
    csf_PWAS$specific == csf_info$specific[i]
  )
  if (length(match_idx) > 0) {
    csf_PWAS$NEW_hit[match_idx] <- csf_info$NEW_hit[i]
  }
}

rm(match_idx, i)

# initialize NEW_hit column in brain_PWAS
brain_PWAS$NEW_hit <- NA

# Iterate through brain_info and update brain_PWAS
for (i in seq_len(nrow(brain_info))) {
  match_idx <- which(
    brain_PWAS$Gene == brain_info$Gene[i] &
    brain_PWAS$NWGT == brain_info$NWGT[i] &
    brain_PWAS$P0 == brain_info$P0[i] &
    brain_PWAS$Discovery == brain_info$Discovery[i] &
    brain_PWAS$specific == brain_info$specific[i]
  )
  if (length(match_idx) > 0) {
    brain_PWAS$NEW_hit[match_idx] <- brain_info$NEW_hit[i]
  }
}

# harmonize dataframes so that they can be concatenated
b = subset(brain_PWAS, select = c(Gene, CHR, TWAS.Z, TWAS.P, Start_hg38, fdr_p:NEW_hit))
c = subset(csf_PWAS, select = c(ID, Gene, CHR, TSS, TWAS.Z, TWAS.P:NEW_hit))

colnames(b) = c("Gene", "CHR", "TWAS.Z", "TWAS.P", "P0", "fdr_p", "specific", "topSNP", "Discovery", "NEW_hit")
b$Analyte = NA
colnames(c) = c("Analyte", "Gene", "CHR", "P0", "TWAS.Z", "TWAS.P", "fdr_p", "specific", "topSNP", "Discovery", "NEW_hit")

b$Tissue = "B"
c$Tissue = "C"

#create new data frame by concatenating brain_PWAS and csf_PWAS dfs
new_pwas_dat = rbind(b, c)
new_pwas_dat$P0 = as.numeric(new_pwas_dat$P0)

#define gap size between chromosomes
gap_size <- 30000000

# Combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
data_cum = new_pwas_dat %>%
group_by(CHR) %>%
dplyr::summarize(max_P0 = max(P0)) %>%
dplyr::mutate(P0_add = lag(cumsum(max_P0) + (1:n() * gap_size), default = 0)) %>% 
dplyr::select(CHR, P0_add)

# Merge above data with original data frame and calculate cumulative bp position for each SNP by adding relative position and adding factor together creating P0_cum column
new_pwas_dat = new_pwas_dat %>%
inner_join(data_cum, by = "CHR") %>%
dplyr::mutate(P0_cum = P0 + P0_add)

# Create a new column for indicating color directly
new_pwas_dat$label_color <- ifelse(new_pwas_dat$NEW_hit == "Y", "#0A09FA",
                                   ifelse(new_pwas_dat$NEW_hit == "N", "#000000", NA))

# update color to add green for genes that have previous evidence in the locus region but not for the specific gene
new_pwas_dat$label_color <- ifelse(!is.na(new_pwas_dat$NEW_hit) & 
                                    (new_pwas_dat$Gene == "TOM1L2" | 
                                    new_pwas_dat$Gene == "INPP5D" |
                                    new_pwas_dat$Gene == "RABEP1" | 
                                    new_pwas_dat$Gene == "CD84" | 
                                    new_pwas_dat$Gene == "ENO3"), 
                                    "#007404",
                                        new_pwas_dat$label_color)

#create new column for indicating nudge_y directly
new_pwas_dat$nudge_y <- ifelse(new_pwas_dat$topSNP == 1, 0.4,
                               ifelse(new_pwas_dat$topSNP == 2, 1, 
                                      ifelse(new_pwas_dat$topSNP == 3, 2, NA)))

# remove CEBPZOS as sex-specific gene in specific & topSNP
new_pwas_dat$specific[new_pwas_dat$Gene == "CEBPZOS"] <- 0
new_pwas_dat$topSNP[new_pwas_dat$Gene == "CEBPZOS"] <- 0

# save

mov_col = "Analyte"
new_pwas_dat = new_pwas_dat %>% relocate(all_of(mov_col), .before =  "Gene")
fwrite(new_pwas_dat, paste0(dir, "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-female-specific-genes.txt"), na = NA, col.names = T, quote = F, row.names = F, sep = '\t')

rm(b, brain_info, brain_PWAS, c, csf_info, csf_PWAS, data_cum, i, match_idx, new_pwas_dat)

############################################
# MALES

brain_PWAS = fread(paste0(brain_dir, m_dir, m_Bf, "_noMHC_ext2Mb_non-strat_sex-strat_W23_weights.txt"), header = T, sep = '\t')
csf_PWAS = fread(paste0(csf_dir, m_dir, m_Cf, "_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_weights.txt"), header = T, sep = '\t')
brain_PWAS = as.data.frame(brain_PWAS); csf_PWAS = as.data.frame(csf_PWAS)

## prepping to plot combined brain and CSF sex-specific results

# add novelty info
csf_info = fread(paste0(csf_dir, m_dir, m_Cf, "_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_top-male-specific-genes_novelty.txt"), header = T, sep = '\t')
brain_info = fread(paste0(brain_dir, m_dir, m_Bf, "_noMHC_ext2Mb_non-strat_sex-strat_W23_top-male-specific-genes_novelty.txt"), header = T, sep = '\t')
csf_info = as.data.frame(csf_info); brain_info = as.data.frame(brain_info)

# Initialize NEW_hit column in csf_PWAS
csf_PWAS$NEW_hit <- NA

# Iterate through csf_info and update csf_PWAS
for (i in seq_len(nrow(csf_info))) {
  match_idx <- which(
    csf_PWAS$ID == csf_info$ID[i] &
    csf_PWAS$Gene == csf_info$Gene[i] &
    csf_PWAS$NWGT == csf_info$NWGT[i] &
    csf_PWAS$TSS == csf_info$TSS[i] &
    csf_PWAS$Discovery == csf_info$Discovery[i] &
    csf_PWAS$specific == csf_info$specific[i])
  
  if (length(match_idx) > 0) {
    csf_PWAS$NEW_hit[match_idx] <- csf_info$NEW_hit[i]
  }
}

rm(match_idx, i)

# initialize NEW_hit column in brain_PWAS
brain_PWAS$NEW_hit <- NA

# Iterate through brain_info and update brain_PWAS
for (i in seq_len(nrow(brain_info))) {
  match_idx <- which(
    brain_PWAS$Gene == brain_info$Gene[i] &
    brain_PWAS$NWGT == brain_info$NWGT[i] &
    brain_PWAS$P0 == brain_info$P0[i] &
    brain_PWAS$Discovery == brain_info$Discovery[i] &
    brain_PWAS$specific == brain_info$specific[i])
  
  if (length(match_idx) > 0) {
    brain_PWAS$NEW_hit[match_idx] <- brain_info$NEW_hit[i]
  }
}

# harmonize dataframes so that they can be concatenated
b = subset(brain_PWAS, select = c(Gene, CHR, TWAS.Z, TWAS.P, Start_hg38, fdr_p:NEW_hit))
c = subset(csf_PWAS, select = c(ID, Gene, CHR, TSS, TWAS.Z, TWAS.P:NEW_hit))

colnames(b) = c("Gene", "CHR", "TWAS.Z", "TWAS.P", "P0", "fdr_p", "specific", "topSNP", "Discovery", "NEW_hit")
b$Analyte = NA
colnames(c) = c("Analyte", "Gene", "CHR", "P0", "TWAS.Z", "TWAS.P", "fdr_p", "specific", "topSNP", "Discovery", "NEW_hit")

b$Tissue = "B"
c$Tissue = "C"

#create new data frame by concatenating brain_PWAS and csf_PWAS dfs
new_pwas_dat = rbind(b, c)
new_pwas_dat$P0 = as.numeric(new_pwas_dat$P0)

# Combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
data_cum = new_pwas_dat %>%
group_by(CHR) %>%
dplyr::summarize(max_P0 = max(P0)) %>%
dplyr::mutate(P0_add = lag(cumsum(max_P0) + (1:n() * gap_size), default = 0)) %>% 
dplyr::select(CHR, P0_add)

# Merge above data with original data frame and calculate cumulative bp position for each SNP by adding relative position and adding factor together creating P0_cum column
new_pwas_dat = new_pwas_dat %>%
inner_join(data_cum, by = "CHR") %>%
dplyr::mutate(P0_cum = P0 + P0_add)

# Create a new column for indicating color directly
new_pwas_dat$label_color <- ifelse(new_pwas_dat$NEW_hit == "Y", "#0A09FA",
    ifelse(new_pwas_dat$NEW_hit == "N", "#000000", NA))

# update color to add green for genes that have previous evidence in the locus region, etc but not for the specific gene
new_pwas_dat$label_color <- ifelse(!is.na(new_pwas_dat$NEW_hit) & 
    (new_pwas_dat$Gene == "PHB1"), "#007404", new_pwas_dat$label_color)

#create new column for indicating nudge_y directly
new_pwas_dat$nudge_y <- ifelse(new_pwas_dat$topSNP == 1, 0.3,
	ifelse(new_pwas_dat$topSNP == 2, 0.3, ifelse(new_pwas_dat$topSNP == 3, 0.2, NA)))

mov_col = "Analyte"
new_pwas_dat = new_pwas_dat %>% relocate(all_of(mov_col), .before =  "Gene")
fwrite(new_pwas_dat, paste0(dir, "ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-male-specific-genes.txt"), na = NA, col.names = T, quote = F, row.names = F, sep = '\t')
