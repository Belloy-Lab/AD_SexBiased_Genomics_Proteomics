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

################################################################################
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

################################################################################
#Read in GWAS summary statistics files
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

# Source bulk colocalization abf
source("functions/bulk_colocalization_abf.R")
