#!/usr/bin/env Rscript
### Master Script for xQTL pipeline ABF - on consistent script

## libraries
library(data.table)
library(plyr)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(locuscomparer)
library(coloc)
library(foreach)
library(doParallel)
library(stringr)
library(R.utils)
library(Rfast)

## Extract reference information: study, tissue, filepath, qtl_type from reference file that is used in bash script
args <- commandArgs(trailingOnly = TRUE)
cat("Arguments:", args, "\n")
study <- args[1]
tissue <- args[2]
filepath <- args[3]
qtl_type <- args[4]

# Clarify all variables are read in properly
cat("Study:", study, "\n")
cat("Tissue:", tissue, "\n")
cat("Filepath:", filepath, "\n")
filepath <- trimws(filepath)
cat("Normalized Filepath (after trimws):", filepath, "\n")
cat(file.exists(filepath), "\n")

###################################
## read in AD samples
Male_rb = fread("GWAS_inputs/Stage1-2-3_EUR_Male_AD_GWAS")
Male_rb = Male_rb %>%
  mutate(Z = BETA/SE,
         s = 37435/507902)

Female_rb = fread("GWAS_inputs/Stage1-2-3_EUR_Female_AD_GWAS")
Female_rb = Female_rb %>%
  mutate(Z = BETA/SE,
         s = 59170/636266)
 
SexHet_rb = fread("GWAS_inputs/Stage1-2-3_EUR_SexHet_AD_GWAS")
SexHet_rb = SexHet_rb %>%
  mutate(Z = BETA/SE,
         s = 96605/1144168) %>%
  dplyr::rename(other_P = P,
                P = SEX_HET_P)

###################################
## read in hg38 gene list with position boundaries
gtf <- rtracklayer::import("reference_files/gencode.v42.basic.annotation.gtf.gz")
ref_seq=as.data.table(gtf)

# Define variables not read in from bash script for abf
gene_list = "Job_scripts/pwas_gwas_xQTL_reference_genes.csv"
LC_dir = "ABF_LCplots/"
LC_threshold = 0.70
out_dir = "ABF/"

## Note: For BrainMeta Datasets - these are curated for the PWAS/GWAS paper - if needed to run QTL on BrainMeta you will need to create a new version of the summary stats
if (study == "BrainMeta" & qtl_type == "caQTL") {
  source("functions/other_code/BrainMeta_caqtl_bulk_coloc.R")
} else if (study == "BrainMeta" & qtl_type == "mQTL") {
  source("functions/other_code/BrainMeta_mqtl_bulk_coloc.R")
} else if (study == "kosoy" & qtl_type == "caQTL") {
  source("functions/other_code/kosoy_caqtl_bulk_coloc.R")
} else if (study == "xQTLServe" & qtl_type == "haQTL") {
  source("functions/other_code/xQTLServe_haqtl_bulk_coloc.R")
} else if (study == "xQTLServe" & qtl_type == "mQTL") {
  source("functions/other_code/xQTLServe_mqtl_bulk_coloc.R")
}  else if (study == "ARIC" & qtl_type == "pQTL") {
  source("functions/other_code/ARIC_pqtl_bulk_coloc.R")
} else if (study == "UKB" & qtl_type == "pQTL") {
  source("functions/other_code/UKB_pqtl_bulk_coloc.R")
} else if (study == "NGI" & tissue == "all_CSF" & qtl_type == "pQTL") {
  source("functions/other_code/NGI_all_pqtl_bulk_coloc.R")
} else if (study == "NGI" & tissue == "female_CSF" & qtl_type == "pQTL") {
  source("functions/other_code/NGI_female_pqtl_bulk_coloc.R")
} else if (study == "NGI" & tissue == "male_CSF" & qtl_type == "pQTL") {
  source("functions/other_code/NGI_male_pqtl_bulk_coloc.R")
} 


