#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(vroom)
  library(dplyr)
})

# Define command-line options
option_list <- list(
  make_option("--gwasdir", type = "character", default = NULL,
              help = "Directory path to GWAS files", metavar = "path"),
  make_option("--female_gwas", type = "character", default = NULL,
              help = "Female GWAS summary statistics filename.", metavar = "ADGC_ADSP_UKB_FinnGen_Females_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var"),
  make_option("--male_gwas", type = "character", default = NULL,
              help = "Male GWAS summary statistics filename.", metavar = "ADGC_ADSP_UKB_FinnGen_Males_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var"),

  make_option("--outdir", type = "character", default = ".",
              help = "Directory to save output .ma files [default: current working directory]", metavar = "outdir"),
  make_option("--female_ma", type = "character", default = NULL,
              help = "Female output base filename (.ma)", metavar = "AD_female.ma"),
  make_option("--male_ma", type = "character", default = NULL,
              help = "Male output base filename (.ma)", metavar = "AD_male.ma")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Validate inputs
if (is.null(opt$gwasdir) || is.null(opt$female_gwas) || is.null(opt$male_gwas) || 
    is.null(opt$outdir) || is.null(opt$female_ma) || is.null(opt$male_ma)) {
  stop("All arguments --gwasdir, --female_gwas, --male_gwas, --outdir, --female_ma, and --male_ma must be provided.\n", call. = FALSE)
}

# trait names
trait_name <- c("Female", "Male")

# trait dir
trait_dir <- c(file.path(opt$gwasdir, opt$female_gwas),
               file.path(opt$gwasdir, opt$male_gwas))

# trait save_name
trait_save_name <- c(file.path(opt$outdir, opt$female_ma),
                     file.path(opt$outdir, opt$male_ma))

for (i in 1:length(trait_name)){
  AD_GWAS <- vroom(trait_dir[i])
  AD_GWAS <- as.data.frame(AD_GWAS)
  AD_GWAS <- select(AD_GWAS, SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P, N_incl)
  AD_GWAS <- distinct(AD_GWAS,AD_GWAS$SNP,.keep_all = T)
  colnames(AD_GWAS) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
  AD_GWAS <- AD_GWAS[,c(1:8)]

  write.table(AD_GWAS, trait_save_name[i],
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, 
              fileEncoding = "UTF-8",  # Save with UTF-8 encoding
              eol = "\n")          # Use UNIX line endings
}
