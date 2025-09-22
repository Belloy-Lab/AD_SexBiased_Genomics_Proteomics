#!/usr/bin/env Rscript

library(data.table)
library(optparse)

# ---- optparse ----
option_list <- list(
  make_option("--in_dir",    type = "character", help = "Directory containing the GWAS summary stats"),
  make_option("--gwas_file", type = "character", help = "GWAS summary stats filename (within --in_dir)"),
  make_option("--out_dir",   type = "character", help = "Output directory"),
  make_option("--out_file",  type = "character", help = "Output filename to write"),
  make_option("--LD_dir",    type = "character", help = "Directory containing merged LD panel table"),
  make_option("--LD_file",   type = "character", help = "Merged LD panel filename (e.g., chr1-22.txt)")
)
opt <- parse_args(OptionParser(option_list = option_list))

# assign variables to match your original code style
dir     <- opt$in_dir
ref_dir <- opt$LD_dir
out_dir <- opt$out_dir

# resolved paths
gwas_path <- file.path(opt$in_dir, opt$gwas_file)
ld_path   <- file.path(opt$LD_dir,  opt$LD_file)
out_path  <- file.path(opt$out_dir, opt$out_file)

# ensure output dir exists
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- GWAS ----
gwas <- fread(gwas_path, header = TRUE, sep = "\t")
gwas <- as.data.frame(gwas)
print("loaded GWAS hg38 meta sumstats"); print(nrow(gwas))

gwas <- na.omit(gwas)
print("nrow of meta sumstats after removing NAs"); print(nrow(gwas))
names(gwas)[2] <- "ID"

# match your original key: CHR:BP:ALLELE0:ALLELE1
gwas$tmpID <- paste(gwas$CHR, gwas$BP, gwas$ALLELE0, gwas$ALLELE1, sep=":")

# ---- LD panel ----
ref_dat <- fread(ld_path, header = FALSE, sep = "\t")
header_names <- c("CHR", "ID", "POS", "BP", "A2", "A1")
colnames(ref_dat) <- header_names

# build both allele-order keys
ref_dat$tmpID  <- paste(ref_dat$CHR, ref_dat$BP, ref_dat$A1, ref_dat$A2, sep=":")
ref_dat$tmpID2 <- paste(ref_dat$CHR, ref_dat$BP, ref_dat$A2, ref_dat$A1, sep=":")

# ---- intersections ----
gwas1 <- gwas[which(gwas$tmpID %in% ref_dat$tmpID), ]
print("nrow of meta sumstats after intersection 1"); print(nrow(gwas1))

gwas2 <- gwas[which(gwas$tmpID %in% ref_dat$tmpID2), ]
print("nrow of meta sumstats after intersection 2"); print(nrow(gwas2))

intersect <- rbind(gwas1, gwas2)
print("total nrow of meta sumstats after intersection"); print(nrow(intersect))

# attach SNP by matching primary key; fallback to swapped key
intersect$SNP <- ifelse(intersect$tmpID %in% ref_dat$tmpID,
                        ref_dat$ID[match(intersect$tmpID, ref_dat$tmpID)],
                        ref_dat$ID[match(intersect$tmpID, ref_dat$tmpID2)])

# ---- output ----
fwrite(intersect, out_path, col.names = TRUE, sep = "\t")
