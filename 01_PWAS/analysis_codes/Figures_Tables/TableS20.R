## table S20: AD PWAS EU all sex specific hits - comparison of those hits to related findings when excluding UKB
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

#############################################################
option_list = list(
  make_option("--dir", action="store", default=NA, type='character',
              help="Path for sex-specific top gene [required]"),
  make_option("--female_top_gene", action="store", default=NA, type='character',
              help="top female  specific  genes file [required]"),
  make_option("--male_top_gene", action="store", default=NA, type='character',
              help="top male  specific  genes file [required]"),
  make_option("--noUKB_brain_dir", action="store", default=NA, type='character',
              help="Path for EU_noUKB brain sex-specific results location [required]"),
  make_option("--noUKB_CSF_dir", action="store", default=NA, type='character',
              help="Path for EU_noUKB CSF sex-specific results location [required]"),
  make_option("--noUKB_female_primary_brain_PWASfile", action="store", default=NA, type='character',
              help="EU_noUKB primary female brain PWAS 'weights' file [required]"),
  make_option("--noUKB_male_primary_brain_PWASfile", action="store", default=NA, type='character',
              help="EU_noUKB primary male brain PWAS 'weights' file [required]"),
  make_option("--noUKB_female_primary_CSF_PWASfile", action="store", default=NA, type='character',
              help="EU_noUKB primary female CSF PWAS 'weights' file [required]"),
  make_option("--noUKB_male_primary_CSF_PWASfile", action="store", default=NA, type='character',
              help="EU_noUKB primary male CSF PWAS 'weights' file [required]"),
  make_option("--Known_locus_file", action="store", default=NA, type='character',
              help="AD_Risk_Loci_consensus_2024.txt document from the supporting folder with full path location [required]")
#  make_option("--Annotation_file", action="store", default=NA, type='character',
#              help="Hg38_Biomart07082025.txt file from the supporting folder with full path location [required]")
)

opt = parse_args(OptionParser(option_list=option_list))
#############################################################

dir = opt$dir # "/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/Final_plots/"
x = fread(paste0(dir, opt$female_top_genel), header = T, sep = '\t') # /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/Final_plots/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-female-specific-genes.txt
y = fread(paste0(dir, opt$male_top_gene), header = T, sep = '\t') # /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/Final_plots/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-male-specific-genes.txt
x = as.data.frame(x); y = as.data.frame(y)

f_dir = "Female/"; m_dir = "Male/"


# EU_noUKB - load sensitivity check meta-analysis PWAS results for both discoveries across sexes
xb <- fread(paste0(opt$noUKB_brain_dir, f_dir, opt$noUKB_female_primary_brain_PWASfile)) # "/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_noUKB/Brain/Sex/Female/ADGC_ADSP_FinnGen_Females_cc.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-strat_sex-strat_W23_weights.txt
yb <- fread(paste0(opt$noUKB_brain_dir, m_dir, opt$noUKB_male_primary_brain_PWASfile)) #"/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_noUKB/Brain/Sex/Male/ADGC_ADSP_FinnGen_Males_cc.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-strat_sex-strat_W23_weights.txt
xc <- fread(paste0(opt$noUKB_CSF_dir, f_dir, opt$noUKB_female_primary_CSF_PWASfile)) # "/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_noUKB/CSF/Sex/Female/ADGC_ADSP_FinnGen_Females_cc.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_weights.txt
yc <- fread(paste0(opt$noUKB_CSF_dir, m_dir, opt$noUKB_male_primary_CSF_PWASfile)) #"/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_noUKB/CSF/Sex/Male/ADGC_ADSP_FinnGen_Males_cc.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_weights.txt
xb = as.data.frame(xb); xc = as.data.frame(xc); yb = as.data.frame(yb); yc = as.data.frame(yc)

LOCI = fread(opt$Known_locus_file, header = T, sep = '\t'); LOCI <- as.data.frame(LOCI)
# ANNO = fread(opt$Annotation_file, header = T, sep = '\t'); ANNO <- as.data.frame(ANNO)


# create dataframe of female-specific genes
a = subset(x, topSNP != 0, select = c(Analyte, Gene, CHR, Tissue))
a = a[order(a$CHR), ]

# create dataframe of male-specific genes and bind with female-specific genes df
b = subset(y, topSNP != 0, select = c(Analyte, Gene, CHR, Tissue))
b = b[order(b$CHR), ]
dat = rbind(a,b)

# Function to update dat with correct values from source_df (x for Female, y for Male)
update_dat_column <- function(dat, source_df, sex_prefix) {
  # Extract relevant columns
  source_filtered <- source_df %>%
    dplyr::select(Analyte, Gene, CHR, Tissue, Discovery, TWAS.Z, TWAS.P, fdr_p)
  
  # Merge data for Primary Discovery
  dat <- dat %>%
    left_join(source_filtered %>% filter(Discovery == "Primary") %>%
       rename_with(~ paste0(sex_prefix, ".PD.", c("Z", "P", "fdr")), .cols = c(TWAS.Z, TWAS.P, fdr_p)),  # nolint: line_length_linter.
       by = c("Analyte", "Gene", "CHR", "Tissue"))

  # Merge data for Secondary Discovery
  dat <- dat %>%
    left_join(source_filtered %>% filter(Discovery == "Secondary") %>% # nolint
                rename_with(~ paste0(sex_prefix, ".SD.", c("Z", "P", "fdr")), .cols = c(TWAS.Z, TWAS.P, fdr_p)),  # nolint: line_length_linter.
              by = c("Analyte", "Gene", "CHR", "Tissue"))

  return(dat)
}

# Apply function separately for females (x) and males (y)
dat <- update_dat_column(dat, x, "Female")
dat <- update_dat_column(dat, y, "Male")

# clean-up dataframe
dat = dplyr::select(dat, -Discovery.x, -Discovery.y, -Discovery.x.x, -Discovery.y.y)

# get best PWAS results for each sex
dat$Female.Z <- ifelse(is.na(dat$Female.PD.P), dat$Female.SD.Z,
                 ifelse(is.na(dat$Female.SD.P), dat$Female.PD.Z,
                        ifelse(dat$Female.PD.P < dat$Female.SD.P, dat$Female.PD.Z, dat$Female.SD.Z))
)
dat$Female.P <- ifelse(is.na(dat$Female.PD.P), dat$Female.SD.P,
                 ifelse(is.na(dat$Female.SD.P), dat$Female.PD.P,
                        ifelse(dat$Female.PD.P < dat$Female.SD.P, dat$Female.PD.P, dat$Female.SD.P))
)
dat$Female.fdr <- ifelse(is.na(dat$Female.PD.P), dat$Female.SD.fdr,
                   ifelse(is.na(dat$Female.SD.P), dat$Female.PD.fdr,
                          ifelse(dat$Female.PD.P < dat$Female.SD.P, dat$Female.PD.fdr,
                                 ifelse(dat$Female.PD.P > dat$Female.SD.P, dat$Female.SD.fdr,
                                        pmin(dat$Female.PD.fdr, dat$Female.SD.fdr, na.rm = TRUE))))
)

dat$Male.Z <- ifelse(is.na(dat$Male.PD.P), dat$Male.SD.Z,
                 ifelse(is.na(dat$Male.SD.P), dat$Male.PD.Z,
                        ifelse(dat$Male.PD.P < dat$Male.SD.P, dat$Male.PD.Z, dat$Male.SD.Z))
)
dat$Male.P <- ifelse(is.na(dat$Male.PD.P), dat$Male.SD.P,
                 ifelse(is.na(dat$Male.SD.P), dat$Male.PD.P,
                        ifelse(dat$Male.PD.P < dat$Male.SD.P, dat$Male.PD.P, dat$Male.SD.P))
)
dat$Male.fdr <- ifelse(is.na(dat$Male.PD.P), dat$Male.SD.fdr,
                   ifelse(is.na(dat$Male.SD.P), dat$Male.PD.fdr,
                          ifelse(dat$Male.PD.P < dat$Male.SD.P, dat$Male.PD.fdr,
                                 ifelse(dat$Male.PD.P > dat$Male.SD.P, dat$Male.SD.fdr,
                                        pmin(dat$Male.PD.fdr, dat$Male.SD.fdr, na.rm = TRUE))))
)

# move columns
col_to_move = c("Female.Z", "Female.P", "Female.fdr")
dat = dat %>% relocate(all_of(col_to_move), .after = "Female.SD.fdr")

# no UKB
# load sensitivity check meta-analysis PWAS results for both discoveries across sexes
xb <- fread("ADGC_ADSP_FinnGen_Females_cc.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-strat_sex-strat_W23_weights.txt")
xc <- fread("ADGC_ADSP_FinnGen_Females_cc.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_weights.txt")
yb <- fread("ADGC_ADSP_FinnGen_Males_cc.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-strat_sex-strat_W23_weights.txt")
yc <- fread("ADGC_ADSP_FinnGen_Males_cc.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_weights.txt")
xb = as.data.frame(xb); xc = as.data.frame(xc); yb = as.data.frame(yb); yc = as.data.frame(yc)

# update sensitivity meta-analysis PWAS dataframes to include the following columns: Analyte, Tissue
xb$Analyte <- NA
yb$Analyte <- NA

xb$Tissue <- "B"
xc$Tissue <- "C"
yb$Tissue <- "B"
yc$Tissue <- "C"

# reduce sensitivity meta-analysis dataframes to necessary columns
xb <- subset(xb, select = c(Analyte, Gene, CHR, Tissue, TWAS.Z, TWAS.P, fdr_p, Discovery))
yb <- subset(yb, select = c(Analyte, Gene, CHR, Tissue, TWAS.Z, TWAS.P, fdr_p, Discovery))

xc <- subset(xc, select = c(ID, Gene, CHR, Tissue, TWAS.Z, TWAS.P, fdr_p, Discovery))
yc <- subset(yc, select = c(ID, Gene, CHR, Tissue, TWAS.Z, TWAS.P, fdr_p, Discovery))
names(xc)[1] <- "Analyte"
names(yc)[1] <- "Analyte"

# rbind sensitivity meta-analysis PWAS results per sex
f_sen <- rbind(xb, xc)
m_sen <- rbind(yb, yc)

add_sensitivity <- function(dat, df, prefix) {
  # merge data for primary discovery
  dat <- dat %>%
    left_join(df %>% filter(Discovery == "Primary") %>%
                rename_with(~ paste0(prefix, ".PD.", c("Z", "P", "fdr")), .cols = c(TWAS.Z, TWAS.P, fdr_p)), # nolint
              by = c("Analyte", "Gene", "CHR", "Tissue"))
  
  # merge data for secondary discovery
  dat <- dat %>%
    left_join(df %>% filter(Discovery == "Secondary") %>%
                rename_with(~ paste0(prefix, ".SD.", c("Z", "P", "fdr")), .cols = c(TWAS.Z, TWAS.P, fdr_p)),
              by = c("Analyte", "Gene", "CHR", "Tissue"))

  return(dat)
}

# apply function separately for females (f_sen) and males (m_sen)
dat <- add_sensitivity(dat, f_sen, "noUKB.F")
dat <- add_sensitivity(dat, m_sen, "noUKB.M")

# clean-up dataframe
dat = dplyr::select(dat, -Discovery.x, -Discovery.y, -Discovery.x.x, -Discovery.y.y)

# get best PWAS results per sex across the sensitivity meta
dat$noUKB.Female.Z <- ifelse(is.na(dat$noUKB.F.PD.P), dat$noUKB.F.SD.Z,
                 ifelse(is.na(dat$noUKB.F.SD.P), dat$noUKB.F.PD.Z,
                        ifelse(dat$noUKB.F.PD.P < dat$noUKB.F.SD.P, dat$noUKB.F.PD.Z, dat$noUKB.F.SD.Z))
)
dat$noUKB.Female.P <- ifelse(is.na(dat$noUKB.F.PD.P), dat$noUKB.F.SD.P,
                 ifelse(is.na(dat$noUKB.F.SD.P), dat$noUKB.F.PD.P,
                        ifelse(dat$noUKB.F.PD.P < dat$noUKB.F.SD.P, dat$noUKB.F.PD.P, dat$noUKB.F.SD.P))
)
dat$noUKB.Female.fdr <- ifelse(is.na(dat$noUKB.F.PD.P), dat$noUKB.F.SD.fdr,
                   ifelse(is.na(dat$noUKB.F.SD.P), dat$noUKB.F.PD.fdr,
                          ifelse(dat$noUKB.F.PD.P < dat$noUKB.F.SD.P, dat$noUKB.F.PD.fdr,
                                 ifelse(dat$noUKB.F.PD.P > dat$noUKB.F.SD.P, dat$noUKB.F.SD.fdr,
                                        pmin(dat$noUKB.F.PD.fdr, dat$noUKB.F.SD.fdr, na.rm = TRUE))))
)

dat$noUKB.Male.Z <- ifelse(is.na(dat$noUKB.M.PD.P), dat$noUKB.M.SD.Z,
                 ifelse(is.na(dat$noUKB.M.SD.P), dat$noUKB.M.PD.Z,
                        ifelse(dat$noUKB.M.PD.P < dat$noUKB.M.SD.P, dat$noUKB.M.PD.Z, dat$noUKB.M.SD.Z))
)
dat$noUKB.Male.P <- ifelse(is.na(dat$noUKB.M.PD.P), dat$noUKB.M.SD.P,
                 ifelse(is.na(dat$noUKB.M.SD.P), dat$noUKB.M.PD.P,
                        ifelse(dat$noUKB.M.PD.P < dat$noUKB.M.SD.P, dat$noUKB.M.PD.P, dat$noUKB.M.SD.P))
)
dat$noUKB.Male.fdr <- ifelse(is.na(dat$noUKB.M.PD.P), dat$noUKB.M.SD.fdr,
                   ifelse(is.na(dat$noUKB.M.SD.P), dat$noUKB.M.PD.fdr,
                          ifelse(dat$noUKB.M.PD.P < dat$noUKB.M.SD.P, dat$noUKB.M.PD.fdr,
                                 ifelse(dat$noUKB.M.PD.P > dat$noUKB.M.SD.P, dat$noUKB.M.SD.fdr,
                                        pmin(dat$noUKB.M.PD.fdr, dat$noUKB.M.SD.fdr, na.rm = TRUE))))
)

# move columns
col_to_move = c("noUKB.Female.Z", "noUKB.Female.P", "noUKB.Female.fdr")
dat = dat %>% relocate(all_of(col_to_move), .after = "noUKB.F.SD.fdr")

## get start and end for each gene. Create variables for using ensembl with human genes for build 38 Ensembl Release 113
mart = useMart("ensembl")
mart = useDataset("hsapiens_gene_ensembl", mart)
ensembl = useEnsembl(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')

#set variable for filter based on input
filters = c("hgnc_symbol")

#set variable for attributes that I would like to receive from the search
attributes = c("hgnc_symbol", "external_gene_name", "entrezgene_accession", "chromosome_name", "start_position","end_position")

genes = dat$Gene

#create variable to run Ensembl search
g1 = getBM(attributes=attributes, filters=filters, values=genes, mart = ensembl)

#create subset where hgnc, external gene name and entrez gene match
g <- g1 %>%
  filter(hgnc_symbol == external_gene_name & hgnc_symbol == entrezgene_accession
); nrow(g)  #39
length(unique(g$hgnc_symbol))   #37 meaning there is a duplicate

# remove rows that have chromosome_name not 1:22
g <- subset(g, chromosome_name %in% c(1:22))

# add start and end columns to dat
dat$START <- g$start_position[match(dat$Gene, g$hgnc_symbol)]
dat$END <- g$end_position[match(dat$Gene, g$hgnc_symbol)]


# Code for Discovery details
dat <- dat %>%
  mutate(Discovery = case_when(
    # FEMALE LOGIC — prioritize if any Female FDR is significant
    (!is.na(Female.PD.fdr) & Female.PD.fdr < 0.05) |
    (!is.na(Female.SD.fdr) & Female.SD.fdr < 0.05) ~ case_when(
      
      !is.na(Female.PD.fdr) & Female.PD.fdr < 0.05 &
        !is.na(Female.SD.fdr) & Female.SD.fdr < 0.05 ~ "Prim+Sec",

      !is.na(Female.PD.fdr) & Female.PD.fdr < 0.05 &
        (is.na(Female.SD.fdr) | Female.SD.fdr >= 0.05) ~ "Primary",

      !is.na(Female.SD.fdr) & Female.SD.fdr < 0.05 &
        (is.na(Female.PD.fdr) | Female.PD.fdr >= 0.05) ~ "Secondary",

      TRUE ~ ""
    ),

    # MALE LOGIC — only if both Female FDRs are non-significant or missing
    (is.na(Female.PD.fdr) | Female.PD.fdr >= 0.05) &
    (is.na(Female.SD.fdr) | Female.SD.fdr >= 0.05) ~ case_when(

      !is.na(Male.PD.fdr) & Male.PD.fdr < 0.05 &
        !is.na(Male.SD.fdr) & Male.SD.fdr < 0.05 ~ "Prim+Sec",

      !is.na(Male.PD.fdr) & Male.PD.fdr < 0.05 &
        (is.na(Male.SD.fdr) | Male.SD.fdr >= 0.05) ~ "Primary",

      !is.na(Male.SD.fdr) & Male.SD.fdr < 0.05 &
        (is.na(Male.PD.fdr) | Male.PD.fdr >= 0.05) ~ "Secondary",

      TRUE ~ ""
    ),

    # FALLBACK: If P-values are NA
    is.na(Female.PD.P) ~ "Secondary",
    is.na(Female.SD.P) ~ "Primary",

    TRUE ~ ""
  ))

# Code for Sensitivity details 
dat <- dat %>%
  mutate(Sensitivity = case_when(
    # Primary: Female-driven with fallback to Male
    Discovery == "Primary" & !is.na(Female.PD.fdr) & Female.PD.fdr < 0.05 &
      !is.na(Female.PD.Z) & !is.na(noUKB.Female.Z) & !is.na(noUKB.Female.P) ~
      ifelse(sign(Female.PD.Z) == sign(noUKB.Female.Z) & noUKB.Female.P < 0.05 &
               !(sign(Female.PD.Z) == sign(noUKB.Male.Z) & noUKB.Male.P < 0.05), "Yes", "No"),

    Discovery == "Primary" & (is.na(Female.PD.fdr) | Female.PD.fdr >= 0.05) &
      !is.na(Male.PD.Z) & !is.na(noUKB.Male.Z) & !is.na(noUKB.Male.P) ~
      ifelse(sign(Male.PD.Z) == sign(noUKB.Male.Z) & noUKB.Male.P < 0.05 &
               !(sign(Male.PD.Z) == sign(noUKB.Female.Z) & noUKB.Female.P < 0.05), "Yes", "No"),

    # Secondary: Female-driven with fallback to Male
    Discovery == "Secondary" & !is.na(Female.SD.fdr) & Female.SD.fdr < 0.05 &
      !is.na(Female.SD.Z) & !is.na(noUKB.Female.Z) & !is.na(noUKB.Female.P) ~
      ifelse(sign(Female.SD.Z) == sign(noUKB.Female.Z) & noUKB.Female.P < 0.05 &
               !(sign(Female.SD.Z) == sign(noUKB.Male.Z) & noUKB.Male.P < 0.05), "Yes", "No"),

    Discovery == "Secondary" & (is.na(Female.SD.fdr) | Female.SD.fdr >= 0.05) &
      !is.na(Male.SD.Z) & !is.na(noUKB.Male.Z) & !is.na(noUKB.Male.P) ~
      ifelse(sign(Male.SD.Z) == sign(noUKB.Male.Z) & noUKB.Male.P < 0.05 &
               !(sign(Male.SD.Z) == sign(noUKB.Female.Z) & noUKB.Female.P < 0.05), "Yes", "No"),

    # Prim+Sec: Female if both fdrs < 0.05; fallback to Male
    Discovery == "Prim+Sec" & !is.na(Female.PD.fdr) & Female.PD.fdr < 0.05 &
      !is.na(Female.SD.fdr) & Female.SD.fdr < 0.05 &
      !is.na(Female.PD.Z) & !is.na(noUKB.Female.Z) & !is.na(noUKB.Female.P) ~
      ifelse(sign(Female.PD.Z) == sign(noUKB.Female.Z) & noUKB.Female.P < 0.05 &
               !(sign(Female.PD.Z) == sign(noUKB.Male.Z) & noUKB.Male.P < 0.05), "Yes", "No"),

    Discovery == "Prim+Sec" & 
      (is.na(Female.PD.fdr) | Female.PD.fdr >= 0.05) &
      (is.na(Female.SD.fdr) | Female.SD.fdr >= 0.05) &
      !is.na(Male.PD.Z) & !is.na(noUKB.Male.Z) & !is.na(noUKB.Male.P) ~
      ifelse(sign(Male.PD.Z) == sign(noUKB.Male.Z) & noUKB.Male.P < 0.05 &
               !(sign(Male.PD.Z) == sign(noUKB.Female.Z) & noUKB.Female.P < 0.05), "Yes", "No"),

    # All others
    TRUE ~ ""
  ))


# Make sure LOCI$GRCh38_POS is numeric
LOCI$GRCh38_POS <- as.numeric(LOCI$GRCh38_POS)

# Function to check if any LOCI position falls within the ±1Mb window for a row in dat
check_novelty <- function(chr, start, end) {
  window_start <- start - 1e6
  window_end   <- end + 1e6
  any(LOCI$CHR == chr & LOCI$GRCh38_POS >= window_start & LOCI$GRCh38_POS <= window_end)
}

# Apply the function row-wise
dat$Novelty <- mapply(function(chr, start, end) {
  if (check_novelty(chr, start, end)) {
    "Known locus"
  } else {
    "Novel gene"
  }
}, dat$CHR, dat$START, dat$END)

fwrite(dat, file = "ADGC_ADSP_UKB_FinnGen_cc_rb.gen090.noAPOE.shared_var.W23nCSF.PWAS_noMHC_ext2Mb_top-sex-specific-genes_sensitivity-check.txt", col.names = T, row.names = F, quote = F, na = NA, sep = '\t')
