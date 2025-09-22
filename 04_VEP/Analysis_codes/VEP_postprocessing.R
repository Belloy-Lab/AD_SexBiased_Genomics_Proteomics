## VEP Post-processing

# Libraries
library(tidyverse)
library(data.table)
library(patchwork)
library(stringr)
library(optparse)

#############################################################
option_list = list(
  make_option("--VEP_dir", type = "character", default = NA,
              help = "Path to the working directory for VEP annotation and related inputs [required]"),
  make_option("--Results_dir", type = "character", default = NA,
              help = "Directory to store all post-processing outputs [required]"),
  make_option("--GWAS_VEP", type = "character", default = NA,
              help = "GWAS VEP annotation output (.txt) downloaded from the Ensembl VEP web portal [required]"),
  make_option("--GWAS_var", type = "character", default = NA,
              help = "CSV file containing top GWAS variants [required]"),
  make_option("--PWAS_VEP", type = "character", default = NA,
              help = "PWAS VEP annotation output (.txt) downloaded from the Ensembl VEP web portal [required]"),
  make_option("--PWAS_var", type = "character", default = NA,
              help = "CSV file containing top PWAS variants [required]")
)

opt = parse_args(OptionParser(option_list=option_list))
#####################################################################################
## GWAS
## Read in GWAS VEP results
gvep = fread(file.path(opt$VEP_dir, opt$GWAS_VEP)) 

gvep = gvep %>% 
  rename(SNP = `#Uploaded_variation`)
n_distinct(gvep$SNP)
# N = 977 unique SNPs


## Read in GWAS input list
gwas = fread(file.path(opt$VEP_dir, opt$GWAS_var)) %>% 
  group_by(locus_index) %>%
  mutate(TOP = if_else(TOP == "", unique(TOP[TOP != ""]), TOP)) %>% # Fill in empty TOP observations to include top variant for each locus
  dplyr::select(locus_index, TOP, SNP, r2, CHR, BP, ALLELE1, ALLELE0, A1FREQ) # select relavent columns for supplemntary table


####################################
## Merge GWAS results with VEP output for GWAS
g_merge = merge(gwas, gvep, by = "SNP") %>% 
  arrange(CHR, BP)
n_distinct(g_merge$SNP)
# N = 977 unique SNPs


# Check to see number of rows/SNPs with corresponding effect alleles
g_a1_match = g_merge %>% 
  filter(ALLELE1 == Allele) %>% 
  arrange(CHR, BP)
n_distinct(g_a1_match$SNP)
# N = 924 unique SNPs


# Check to see number of rows/SNPs with non-corresponding effect alleles
g_a1_nonmatch = g_merge %>% 
  filter(ALLELE1 != Allele) %>% 
  distinct(SNP, .keep_all = TRUE) %>% # distinct to remove duplicate SNPs
  arrange(CHR, BP)
n_distinct(g_a1_nonmatch$SNP)
# N = 53 unique SNPs


# Write out 3 tables for supplementary - talk to Michael about non-match (53 varaints for MAPT, 3 variants for NCK2)
fwrite(g_merge, file.path(opt$Results_dir, "gwas_vep_output_all_snps.csv"))
fwrite(g_a1_match, file.path(opt$Results_dir, "gwas_vep_output_matching_EA.csv"))
fwrite(g_a1_nonmatch, file.path(opt$Results_dir, "gwas_vep_output_non-matching_EA.csv"))

rm(list = ls())
## End of Code for GWAS Postprocessing

##########################################################################################################################################################################
## PWAS
## Read in PWAS VEP results
pvep = fread(file.path(opt$VEP_dir, opt$PWAS_VEP)) %>% 
  rename(SNP = `#Uploaded_variation`)
n_distinct(pvep$SNP)
# N = 1592 unique SNPs


## Read in PWAS input list
pwas = fread(file.path(opt$VEP_dir, opt$PWAS_var)) %>% 
  group_by(locus_index) %>%
  mutate(TOP = if_else(r2 == 1, SNP, NA_character_)) %>%
  fill(TOP, .direction = "downup") %>%
  ungroup()%>% # Fill in empty TOP observations to include top variant for each locus
  dplyr::select(discovery, locus_index, TOP, SNP, r2, CHR, BP, ALLELE1, ALLELE0) # select relavent columns for supplemntary table



####################################
## Merge PWAS results with VEP output for PWAS
p_merge = merge(pwas, pvep, by = "SNP") %>% 
  arrange(CHR, BP)
n_distinct(p_merge$SNP)
# N = 1592 unique SNPs


# Check to see number of rows/SNPs with corresponding effect alleles
p_a1_match = p_merge %>% 
  filter(ALLELE1 == Allele) %>% 
  arrange(CHR, BP)
n_distinct(p_a1_match$SNP)
# N = 1591 unique SNPs


# Check to see number of rows/SNPs with non-corresponding effect alleles
p_a1_nonmatch = p_merge %>% 
  filter(ALLELE1 != Allele) %>% 
  distinct(SNP, .keep_all = TRUE) %>% # distinct to remove duplicate SNPs
  arrange(CHR, BP)
n_distinct(p_a1_nonmatch$SNP)
# N = 1 unique SNPs


# Write out 3 tables for supplementary - talk to Michael about non-match (53 varaints for MAPT, 3 variants for NCK2)
fwrite(p_merge, file.path(opt$Results_dir, "pwas_vep_output_all_snps.csv"))
fwrite(p_a1_match, file.path(opt$Results_dir, "pwas_vep_output_matching_EA.csv"))
fwrite(p_a1_nonmatch, file.path(opt$Results_dir, "pwas_vep_output_non-matching_EA.csv"))

rm(list = ls())

##########################################################################################################################################################################
# Create Figure for PWAS VEP results (Version 1 - non-collapsed loci)
p_merge = fread(file.path(opt$Results_dir, "pwas_vep_output_all_snps.csv"))

p_fig = p_merge %>% 
  filter(IMPACT != "LOW") %>% # Remove LOW impact variants
  filter(SYMBOL == "-" | SYMBOL == locus_index) # Filter to only include variants with gene symbols that match locus index or are unannotated

n_distinct(p_fig$SNP)
# N = 856 unique SNPs

n_distinct(p_fig$locus_index)
# N = 34 unique SNPs
# Lose MINDY1, RASA4B, and RCN2


## Create count df for plotting matrix
# Get unique locus_index and consequence values
locus_index_unique = unique(p_fig$locus_index)
consequence_unique = unique(p_fig$Consequence)

# Create a cross join dataframe with all combinations of locus_index and consequence
p_counts = expand.grid(locus_index = locus_index_unique, Consequence = consequence_unique)

# Add the count column with observations 0 or 1
p_counts = p_counts %>%
  left_join(p_fig %>% 
              distinct(locus_index, Consequence) %>% 
              mutate(count = 1),  
            by = c("locus_index", "Consequence")) %>%
  mutate(count = ifelse(is.na(count), 0, count))

## Get CHR and BP for each locus_index for plotting positions
c_bp = p_fig %>% 
  distinct(locus_index, .keep_all = TRUE) %>% 
  dplyr::select(locus_index, CHR, BP)

# Merge CHR and BP into count df
p_counts = merge(p_counts, c_bp, by = "locus_index") %>% 
  arrange(CHR, BP)


# Write out count df - this file has 476 rows (34 x 14)
fwrite(p_counts, file.path(opt$Results_dir, "pwas_vep_output_counts.csv"))

####################################################################################################################
# In excel - manuualy add rows for MINDY1, RASA4B, and RCN2 - with 0 observations for all consequence in these 3 loci
# New file has 518 rows (37 x 14)
p_fig2 = fread(file.path(opt$Results_dir, "pwas_vep_output_counts.csv"))


## Need to collapse consequence observations to single consequences
# First, split Consequence into individual consequences where commas exist
split_consequences = p_fig2 %>%
  separate_rows(Consequence, sep = ",") %>%
  group_by(locus_index, Consequence) %>%
  summarize(count = max(count), CHR = max(CHR), BP = max(BP), .groups = 'drop')

# Get the unique consequences to ensure we're working with the correct values
unique_consequences = unique(split_consequences$Consequence)

# Prepare combinations of unique locus_index and unique consequences
collapse_df = expand.grid(locus_index = unique(p_fig2$locus_index), Consequence = unique_consequences)

# Aggregate counts per those combinations
collapse_df = collapse_df %>%
  left_join(split_consequences, by = c("locus_index", "Consequence")) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  arrange(CHR, BP)

# Order locus_index according to CHR and BP
ordered_locus_index <- collapse_df %>%
  arrange(CHR, BP) %>%
  pull(locus_index) %>%
  unique()

# Update collapse_df with ordered locus_index
collapse_df <- collapse_df %>%
  mutate(locus_index = factor(locus_index, levels = ordered_locus_index))

# Calculate total counts for each Consequence
total_counts <- collapse_df %>%
  group_by(Consequence) %>%
  summarize(total_count = sum(count))

# Merge total_counts with collapse_df for plotting purposes
collapse_df <- collapse_df %>%
  left_join(total_counts, by = "Consequence") %>%
  mutate(locus_index = factor(locus_index, levels = unique(locus_index)))


# Create the matrix plot using ggplot2
p_matrix = ggplot(collapse_df, aes(x = locus_index, y = Consequence, fill = as.factor(count))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none",
        plot.margin = unit(c(1, 4, 1, 1), "cm")) + # Adjust margins to fit the text
  labs(
    x = "Locus Index",
    y = "Consequence",
    title = "Matrix of Consequences by Locus Index") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

print(p_matrix)

# Save version 1a of matrix
file_path <- file.path(opt$Results_dir, "PWAS_VEP_matrix_1a.jpg")
ggsave(filename = file_path, plot = p_matrix, device = "jpg", width = 11, height = 9, dpi = 320)

rm(list = ls())

####################################################################################################################
# Create Figure (Version 2 - collapsed loci)
# Collapse PWAS Genes in loci
p_merge = fread(file.path(opt$Results_dir, "pwas_vep_output_all_snps.csv")) %>% 
  mutate(locus_name = case_when(
    locus_index == "TDRKH" ~ "TDRKH/MINDY1",
    locus_index == "MINDY1" ~ "TDRKH/MINDY1",
    locus_index == "USP4" ~ "USP4/USP19",
    locus_index == "USP19" ~ "USP4/USP19",
    locus_index == "ANXA11" ~ "TSPAN14/ANXA11",
    locus_index == "TSPAN14" ~ "TSPAN14/ANXA11",
    locus_index == "ETFA" ~ "ETFA/SCAPER/RCN2",
    locus_index == "SCAPER" ~ "ETFA/SCAPER/RCN2",
    locus_index == "RCN2" ~ "ETFA/SCAPER/RCN2",
    locus_index == "ENO3" ~ "ENO3/RABEP1",
    locus_index == "RABEP1" ~ "ENO3/RABEP1",
    locus_index == "SLC44A2" ~ "SLC44A2/CARM1",
    locus_index == "CARM1" ~ "SLC44A2/CARM1",
    TRUE ~ locus_index
  ))

# Apply filter to keep all genes within locus_name matching Symbol
p_fig = p_merge %>% 
  filter(IMPACT != "LOW") %>% # Remove LOW impact variants
  filter(SYMBOL == "-" | sapply(SYMBOL, function(sym) any(str_detect(locus_name, sym))))

n_distinct(p_fig$SNP)
# N = 902 unique SNPs

n_distinct(p_fig$locus_name)
# N = 29 unique SNPs
# Lose RASA4B 

# Get unique locus_index and consequence values
locus_index_unique = unique(p_fig$locus_name)
consequence_unique = unique(p_fig$Consequence)

# Create a cross join dataframe with all combinations of locus_index and consequence
p_counts = expand.grid(locus_name = locus_index_unique, Consequence = consequence_unique)

# Add the count column with observations 0 or 1
p_counts = p_counts %>%
  left_join(p_fig %>% 
              distinct(locus_name, Consequence) %>% 
              mutate(count = 1),  
            by = c("locus_name", "Consequence")) %>%
  mutate(count = ifelse(is.na(count), 0, count))

## Get CHR and BP for each locus_index for plotting positions
c_bp = p_fig %>% 
  distinct(locus_name, .keep_all = TRUE) %>% 
  dplyr::select(locus_name, CHR, BP)

# Merge CHR and BP into count df
p_counts = merge(p_counts, c_bp, by = "locus_name") %>% 
  arrange(CHR, BP)

# Write out count df - this file has 406 rows (29 x 14)
fwrite(p_counts, file.path(opt$Results_dir, "pwas_vep_output_counts_collapsed_loci.csv"))

rm(list = ls())

####################################################################################################################
# In excel - manually add rows for RASA4B - with 0 observations for all consequence in thic locus
# New file has 420 rows (30 x 14)
p_fig2 = fread(file.path(opt$Results_dir, "pwas_vep_output_counts_collapsed_loci.csv"))

## Need to collapse consequence observations to single consequences
# First, split Consequence into individual consequences where commas exist
split_consequences = p_fig2 %>%
  separate_rows(Consequence, sep = ",") %>%
  group_by(locus_name, Consequence) %>%
  summarize(count = max(count), CHR = max(CHR), BP = max(BP), .groups = 'drop')

# Get the unique consequences to ensure we're working with the correct values
unique_consequences = unique(split_consequences$Consequence)

# Prepare combinations of unique locus_name and unique consequences
collapse_df = expand.grid(locus_name = unique(p_fig2$locus_name), Consequence = unique_consequences)

# Aggregate counts per those combinations
collapse_df = collapse_df %>%
  left_join(split_consequences, by = c("locus_name", "Consequence")) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  arrange(CHR, BP)

# Order locus_name according to CHR and BP
ordered_locus_name <- collapse_df %>%
  arrange(CHR, BP) %>%
  pull(locus_name) %>%
  unique()

# Update collapse_df with ordered locus_name
collapse_df <- collapse_df %>%
  mutate(locus_name = factor(locus_name, levels = ordered_locus_name))

# Calculate total counts for each Consequence
total_counts <- collapse_df %>%
  group_by(Consequence) %>%
  summarize(total_count = sum(count))

# Merge total_counts with collapse_df for plotting purposes
collapse_df <- collapse_df %>%
  left_join(total_counts, by = "Consequence") %>%
  mutate(locus_name = factor(locus_name, levels = unique(locus_name)))

# Create the matrix plot using ggplot2
p_matrix = ggplot(collapse_df, aes(x = locus_name, y = Consequence, fill = as.factor(count))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none",
        plot.margin = unit(c(1, 4, 1, 1), "cm")) + # Adjust margins to fit the text
  labs(
    x = "PWAS Locus",
    y = "Consequence",
    title = "Matrix of Consequences by Locus Index") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

print(p_matrix)

# Save version 2a of matrix
file_path <- file.path(opt$Results_dir, "PWAS_VEP_matrix_2a.jpg")
ggsave(filename = file_path, plot = p_matrix, device = "jpg", width = 11, height = 9, dpi = 320)

## End of Code for PWAS Postprocessing