## VEP Pre-processing 

# Libraries
library(tidyverse)
library(data.table)
library(dplyr)
library(optparse)

#############################################################
option_list = list(
  make_option("--VEP_dir", type = "character", default = NA,
              help = "Path to the working directory for VEP analysis [required]"),

  make_option("--GWAS_dir", type = "character", default = NA,
              help = "Directory containing GWAS significant SNP lists [required]"),
  make_option("--GWAS_file", type = "character", default = NA,
              help = "Filename of the GWAS significant SNP list in CSV format [required]"),
  make_option("--GWAS_LD_dir", type = "character", default = NA,
              help = "Directory containing LD matrices for GWAS SNPs [required]"),

  make_option("--PWAS_dir", type = "character", default = NA,
              help = "Directory containing significant PWAS results [required]"),

  make_option("--PWAS_F_Brain", type = "character", default = NA,
              help = "PWAS LD file for significant genes in female brain [required]"),
  make_option("--PWAS_M_Brain", type = "character", default = NA,
              help = "PWAS LD file for significant genes in male brain [required]"),
  make_option("--PWAS_F_CSF", type = "character", default = NA,
              help = "PWAS LD file for significant genes in female CSF [required]"),
  make_option("--PWAS_M_CSF", type = "character", default = NA,
              help = "PWAS LD file for significant genes in male CSF [required]"),

  make_option("--PWAS_X3504_FC", type = "character", default = NA,
              help = "LD file for X3504.1 gene from female CSF PWAS [required]"),
  make_option("--PWAS_TOM1L2_FB", type = "character", default = NA,
              help = "LD file for TOM1L2 gene from female brain PWAS [required]")
)

opt = parse_args(OptionParser(option_list=option_list))
###################################################################################
### AD GWAS Variant list
gwas = fread(file.path(opt$VEP_Dir, opt$GWAS_dir, opt$GWAS_file), header = T, sep = '\t')

# Create Top hit file
top = gwas %>% 
  filter(!is.na(TOP)) %>% # Filter to top hits
  filter(GENE != "CR1" & GENE != "MINDY2" & GENE != "CTSH" & GENE != "MS4A4A")  %>% # Filter out variants that did not pass sex-het
  dplyr::select(SNP, posID, CHR, BP, ALLELE1, ALLELE0, A1FREQ, IND, GENE, NEW_hit, LOC_con) # Select relevant columns

# Write out top variants
fwrite(top, file.path(opt$VEP_Dir, opt$GWAS_dir, "GWAS_top_20_hits.csv"))

# Define location of LD matrices
ld_dir = file.path(opt$VEP_Dir, opt$GWAS_LD_dir) # "GWAS/LD/"

# Set empty dataframe for GWAS variants
gwas_vars = data.frame()

# Loop through each row in top df - extract variant of interest -
for (i in 1:nrow(top)){
  # Define variables needed
  ind = top$IND[i]
  locus = top$LOC_con[i]
  rsID = top$SNP[i]

  # read in LD matrix based off of ind variable
  ld_matrix = fread(paste0(ld_dir, "/LD_matrixIND=", ind, ".ld.csv")) %>% 
    select(V1, rsID) %>% 
    rename(rsIDs = V1, r2 = rsID) %>% 
    mutate(locus_index = locus)

  # Append variants to final list
  gwas_vars = rbind(gwas_vars, ld_matrix)
}

# Filter gwas_vars to only include variatsn with r2 > 0.8
gwas_vars_filt = gwas_vars %>% 
  filter(r2 > 0.8) %>% 
  rename(SNP = rsIDs)

# Merge comprehensive results with variants that passed r2 filter
gwas_merged = inner_join(gwas, gwas_vars_filt, by = "SNP") %>% 
  arrange(CHR, BP)

# Check to make sure unique rsIDs match nrow of gwas_merged
nrow(gwas_merged) == n_distinct(gwas_merged$SNP)

# Write out merged/filtered results
fwrite(gwas_merged, "Variants_list/gwas_final_variants.csv")

####################################################################################################################################################################################
# Convert gwas_merged file to vcf format for VEP analysis

# Create the VCF header
vcf_header <- c(
  "##fileformat=VCFv4.0",
  "##source=YourDataSource",
  "##reference=hg38",
  "##INFO=<ID=ID,Number=1,Type=String,Description=\"Identifier\">",
  "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|Cytoband|ENSP|SWISSPROT|TREMBL|UNIPARC|NEAREST|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
)

# Select relavent columns from GWAS merged
csv_data = gwas_merged %>% 
  dplyr::select(posID, SNP, CHR, BP, ALLELE1, ALLELE0)

# Format the CSV data for VCF
vcf_data <- csv_data %>%
  mutate(CHROM = CHR,
         POS = BP,
         ID = SNP,
         REF = ALLELE0,
         ALT = ALLELE1,
         QUAL = ".",
         FILTER = ".",
         INFO = ".") %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)


# Write VCF header and data to a file
vcf_output_file <- "Variants_list/gwas_final_variants.vcf"
writeLines(vcf_header, vcf_output_file)
write.table(vcf_data, file = vcf_output_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, append = TRUE)

rm(list = ls())
########################################

### PWAS
pwas_sources <- list(
  list(path = file.path(op$VEP_dir, opt$PWAS_dir, opt$PWAS_F_Brain), discovery = "Female Brain"),
  list(path = file.path(op$VEP_dir, opt$PWAS_dir, opt$PWAS_M_Brain), discovery = "Male Brain"),
  list(path = file.path(op$VEP_dir, opt$PWAS_dir, opt$PWAS_F_CSF), discovery = "Female CSF"),
  list(path = file.path(op$VEP_dir, opt$PWAS_dir, opt$PWAS_M_CSF), discovery = "Male CSF"),
  list(path = file.path(op$VEP_dir, opt$PWAS_dir, opt$PWAS_X3504_FC), discovery = "Female CSF"),
  list(path = file.path(op$VEP_dir, opt$PWAS_dir, opt$PWAS_TOM1L2_FB), discovery = "Female Brain")
)

# Function to process one source
process_source <- function(src) {
  df <- fread(src$path)
  genes <- unique(df$Gene_name2)

# build a data.frame for each gene and bind them together
result <- lapply(genes, function(gene) {
  df %>%
    filter(Gene_name2 == gene, top_pQTL_r2 != "") %>%
    rename(
      SNP       = top_pQTL_r2_snp,
      r2        = top_pQTL_r2,
      locus_index = Gene_name2,
      CHR       = top_pQTL_r2_chr,
      BP        = top_pQTL_r2_bp,
      ALLELE1   = top_pQTL_r2_A1,
      ALLELE0   = top_pQTL_r2_A2
    ) %>%
    mutate(
      discovery = src$discovery,
      posID     = paste(CHR, BP, ALLELE0, ALLELE1, sep = ":")
    ) %>%
    select(CHR, BP, ALLELE1, ALLELE0, posID, SNP, r2, locus_index, discovery)
}) %>%
  bind_rows() %>%
  arrange(CHR, BP)

# sanity check
stopifnot(n_distinct(df$Gene_name2) == n_distinct(result$locus_index))
  result
}

# Process all sources and merge
pwas_list   <- lapply(pwas_sources, process_source)
pwas_merged <- bind_rows(pwas_list) %>%
  filter(!locus_index %in% c("CLIC1", "PACSIN1")) %>%  # drop problematic genes
  arrange(CHR, BP)

# Final integrity checks
stopifnot(n_distinct(pwas_merged$SNP) == nrow(pwas_merged))
stopifnot(n_distinct(pwas_merged$locus_index) == 37)
stopifnot(
  nrow(filter(pwas_merged, r2 == 1)) ==
  n_distinct(filter(pwas_merged, r2 == 1)$locus_index)
)

# Write out merged/filtered results
fwrite(pwas_merged, "Variants_list/pwas_final_variants.csv")

####################################################################################################################################################################################
# Convert pwas_merged file to vcf format for VEP analysis
vcf_header <- c(
  "##fileformat=VCFv4.0",
  "##source=YourDataSource",
  "##reference=hg38",
  "##INFO=<ID=ID,Number=1,Type=String,Description=\"Identifier\">",
  "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|Cytoband|ENSP|SWISSPROT|TREMBL|UNIPARC|NEAREST|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
)

# Select relavent columns from GWAS merged
csv_data = pwas_merged %>% 
  dplyr::select(posID, SNP, CHR, BP, ALLELE1, ALLELE0)

# Format the CSV data for VCF
vcf_data <- csv_data %>%
  mutate(CHROM = CHR,
         POS = BP,
         ID = SNP,
         REF = ALLELE0,
         ALT = ALLELE1,
         QUAL = ".",
         FILTER = ".",
         INFO = ".") %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)

# Write VCF header and data to a file
vcf_output_file <- "Variants_list/pwas_final_variants.vcf"
writeLines(vcf_header, vcf_output_file)
write.table(vcf_data, file = vcf_output_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, append = TRUE)
## End of Code for PWAS Preprocessing
