#!/usr/bin/env Rscript

suppressPackageStartupMessages({
 library(dplyr)
})

# Define command-line options
option_list <- list(
  make_option("--flist_dir", type = "character", help = "Directory path to both CSF and Brain flists"),
  make_option("--flist_CSF", type = "character", help = "CSF flist file name"),
  make_option("--flist_Brain", type = "character", help = "Brain flist file name"),
  make_option("--flist_combined", type = "character", help = "User defined combined flist file name"),
  make_option("--output_dir", type = "character", default = ".", help = "Output directory [default: current dir]")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

my_flist_NGI <- read.table(file.path(opt$flist_dir, opt$flist_CSF), header = T)
my_flist_wingo <- read.table(file.path(opt$flist_dir, opt$flist_Brain), header = T)
my_flist <- rbind(my_flist_NGI,my_flist_wingo)

# Save as Linux-readable format
write.table(my_flist, file = file.path(opt$output_dir, opt$flist_combined),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8", eol = "\n")

# data has been prepared, now running the smr to create besd-format files
smr_exe <- "SMR/amalysis_codes/smr"
com_flist <- file.path(opt$output_dir, opt$flist_combined)
smr_output <- file.path(opt$output_dir, "CSF_Brain_pQTL_besd")

# Make binary executable (if needed)
system(sprintf("chmod +x %s", smr_exe))

# Run SMR to generate BESD files
system(sprintf("%s --eqtl-flist %s --make-besd-dense --out %s", smr_exe, com_flist, smr_output))