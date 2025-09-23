#### create table S23
## table 2: AD PWAS genes + ancestry consistency + SMR + COLOC with tailored pairs
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
  make_option("--dir", action="store", default=NA, type='character',
              help="Path for sex-specific top gene [required]"),
  make_option("--female_top_gene", action="store", default=NA, type='character',
              help="top female  specific  genes file [required]"),
  make_option("--male_top_gene", action="store", default=NA, type='character',
              help="top male  specific  genes file [required]"),
			  
  make_option("--PD_dir", action="store", default=NA, type='character',
              help="Path for EU_all PWAS results location [required]"),
			  
  make_option("--fb_combZ", action="store", default=NA, type='character',
              help="Combine PWAS Z-scores for AFR & EU_all female brain [required]"),
  make_option("--mb_combZ", action="store", default=NA, type='character',
              help="Combine PWAS Z-scores for AFR & EU_all male brain [required]"),
  make_option("--fc_combZ", action="store", default=NA, type='character',
              help="Combine PWAS Z-scores for AFR & EU_all female CSF [required]"),
  make_option("--mc_combZ", action="store", default=NA, type='character',
              help="Combine PWAS Z-scores for AFR & EU_all male CSF [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

#############################################################

dir = opt$dir # "/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/Figures_Tables/"
x = fread(paste0(dir, opt$female_top_genel), header = T, sep = '\t') 
y = fread(paste0(dir, opt$male_top_gene), header = T, sep = '\t') 
x = as.data.frame(x); y = as.data.frame(y)
fb_dir = "Brain/Sex/Female/"; fc_dir = "CSF/Sex/Female/"
mb_dir = "Brain/Sex/Male/"; mc_dir = "CSF/Sex/Male/"

#AFR primary analysis PWAS results from both brain and CSF

w_AFR = fread(paste0(opt$PD_dir, fb_dir, opt$fb_combZ))
x_AFR = fread(paste0(opt$PD_dir, mb_dir, opt$mb_combZ))
y_AFR = fread(paste0(opt$PD_dir, fc_dir, opt$fc_combZ))
z_AFR = fread(paste0(opt$PD_dir, mc_dir, opt$mc_combZ))
w_AFR = as.data.frame(w_AFR); x_AFR = as.data.frame(x_AFR); y_AFR = as.data.frame(y_AFR); z_AFR = as.data.frame(z_AFR)

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
                rename_with(~ paste0(sex_prefix, ".PD.", c("Z", "P", "fdr")), .cols = c(TWAS.Z, TWAS.P, fdr_p)), 
              by = c("Analyte", "Gene", "CHR", "Tissue"))

  # Merge data for Secondary Discovery
  dat <- dat %>%
    left_join(source_filtered %>% filter(Discovery == "Secondary") %>% 
                rename_with(~ paste0(sex_prefix, ".SD.", c("Z", "P", "fdr")), .cols = c(TWAS.Z, TWAS.P, fdr_p)), 
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

## AFR
replace_genes <- function(df, column, replacements) {
  # Update replacements to use word boundaries for regex
  replacements_with_boundaries <- setNames(
    replacements,
    paste0("\\b", names(replacements), "\\b")
  )
  
  df %>% 
    mutate(!!sym(column) := str_replace_all(.data[[column]], replacements_with_boundaries))
}

# Define the replacement rules as a named vector
source("/storage2/fs1/belloy2/Active/04_Code/username/PWAS/EU_all/Wingo23_NGI-CSFc_update_gene_names.R")

dfs <- list(w_AFR = w_AFR, x_AFR = x_AFR, y_AFR = y_AFR, z_AFR = z_AFR)
col <- c("Gene")

# Apply the function to each data frame and each column if it exists
for (df_name in names(dfs)) {
  for (col in col) {
    if (col %in% colnames(dfs[[df_name]])) {
      dfs[[df_name]] <- replace_genes(dfs[[df_name]], col, replacements)
    }
  }
}

# Unpack modified data frames back to individual variables
list2env(dfs, .GlobalEnv)

# subset the AFR data to the PWAS z-score and p-value per discovery and sex
w_AFR = subset(w_AFR, select = c(Gene, AFR.Female.PD.PWAS.Z:AFR.Male.SD.PWAS.P, comb.Male.PD.PWAS.Z:comb.Male.SD.PWAS.P))
x_AFR = subset(x_AFR, select = c(Gene, AFR.Female.PD.PWAS.Z:AFR.Female.SD.PWAS.P, comb.Female.PD.PWAS.Z:comb.Female.SD.PWAS.P, AFR.Male.PD.PWAS.Z:comb.Male.SD.PWAS.P))
w_AFR$Analyte = NA; x_AFR$Analyte = NA
y_AFR = subset(y_AFR, select = c(ID, Gene, AFR.Female.PD.PWAS.Z:AFR.Male.SD.PWAS.P, comb.Male.PD.PWAS.Z:comb.Male.SD.PWAS.P))
z_AFR = subset(z_AFR, select = c(ID, Gene, AFR.Female.PD.PWAS.Z:AFR.Female.SD.PWAS.P, comb.Female.PD.PWAS.Z:comb.Female.SD.PWAS.P, AFR.Male.PD.PWAS.Z:comb.Male.SD.PWAS.P))
names(y_AFR)[1] = "Analyte"; names(z_AFR)[1] = "Analyte"

d = rbind(w_AFR, x_AFR, y_AFR, z_AFR)

d = subset(d, !(Gene %in% c("CEBPZOS", "PACSIN1", "CLIC1")))

dat2 = merge(dat, d, by = c("Analyte", "Gene"))
dat2 = dat2[order(dat2$CHR), ]

# obtain best z and p for AFR per sex
dat2$AFR.Female.Z <- ifelse(is.na(dat2$AFR.Female.PD.PWAS.P), dat2$AFR.Female.SD.PWAS.Z,
                 ifelse(is.na(dat2$AFR.Female.SD.PWAS.P), dat2$AFR.Female.PD.PWAS.Z,
                        ifelse(dat2$AFR.Female.PD.PWAS.P < dat2$AFR.Female.SD.PWAS.P, dat2$AFR.Female.PD.PWAS.Z, dat2$AFR.Female.SD.PWAS.Z))
)
dat2$AFR.Female.P <- ifelse(is.na(dat2$AFR.Female.PD.PWAS.P), dat2$AFR.Female.SD.PWAS.P,
                 ifelse(is.na(dat2$AFR.Female.SD.PWAS.P), dat2$AFR.Female.PD.PWAS.P,
                        ifelse(dat2$AFR.Female.PD.PWAS.P < dat2$AFR.Female.SD.PWAS.P, dat2$AFR.Female.PD.PWAS.P, dat2$AFR.Female.SD.PWAS.P))
)

dat2$AFR.Male.Z <- ifelse(is.na(dat2$AFR.Male.PD.PWAS.P), dat2$AFR.Male.SD.PWAS.Z,
                 ifelse(is.na(dat2$AFR.Male.SD.PWAS.P), dat2$AFR.Male.PD.PWAS.Z,
                        ifelse(dat2$AFR.Male.PD.PWAS.P < dat2$AFR.Male.SD.PWAS.P, dat2$AFR.Male.PD.PWAS.Z, dat2$AFR.Male.SD.PWAS.Z))
)
dat2$AFR.Male.P <- ifelse(is.na(dat2$AFR.Male.PD.PWAS.P), dat2$AFR.Male.SD.PWAS.P,
                 ifelse(is.na(dat2$AFR.Male.SD.PWAS.P), dat2$AFR.Male.PD.PWAS.P,
                        ifelse(dat2$AFR.Male.PD.PWAS.P < dat2$AFR.Male.SD.PWAS.P, dat2$AFR.Male.PD.PWAS.P, dat2$AFR.Male.SD.PWAS.P))
)

# move columns
col_to_move = c("AFR.Female.Z", "AFR.Female.P")
dat2 = dat2 %>% relocate(all_of(col_to_move), .after = "AFR.Female.SD.PWAS.P")

# obtain best comb p and z per sex
dat2$comb.Female.Z <- ifelse(is.na(dat2$comb.Female.PD.PWAS.P), dat2$comb.Female.SD.PWAS.Z,
                  ifelse(is.na(dat2$comb.Female.SD.PWAS.P), dat2$comb.Female.PD.PWAS.Z,
                          ifelse(abs(dat2$comb.Female.PD.PWAS.Z) > abs(dat2$comb.Female.SD.PWAS.Z), dat2$comb.Female.PD.PWAS.Z, dat2$comb.Female.SD.PWAS.Z))
)

dat2$comb.Female.P <- ifelse(is.na(dat2$comb.Female.PD.PWAS.P), dat2$comb.Female.SD.PWAS.P,
                  ifelse(is.na(dat2$comb.Female.SD.PWAS.P), dat2$comb.Female.PD.PWAS.P,
                          ifelse(abs(dat2$comb.Female.PD.PWAS.Z) > abs(dat2$comb.Female.SD.PWAS.Z), dat2$comb.Female.PD.PWAS.P, dat2$comb.Female.SD.PWAS.P))
)

dat2$comb.Male.Z <- ifelse(is.na(dat2$comb.Male.PD.PWAS.P), dat2$comb.Male.SD.PWAS.Z,
                  ifelse(is.na(dat2$comb.Male.SD.PWAS.P), dat2$comb.Male.PD.PWAS.Z,
                          ifelse(abs(dat2$comb.Male.PD.PWAS.Z) > abs(dat2$comb.Male.SD.PWAS.Z), dat2$comb.Male.PD.PWAS.Z, dat2$comb.Male.SD.PWAS.Z))
)

dat2$comb.Male.P <- ifelse(is.na(dat2$comb.Male.PD.PWAS.P), dat2$comb.Male.SD.PWAS.P,
                  ifelse(is.na(dat2$comb.Male.SD.PWAS.P), dat2$comb.Male.PD.PWAS.P,
                          ifelse(abs(dat2$comb.Male.PD.PWAS.Z) < abs(dat2$comb.Male.SD.PWAS.Z), dat2$comb.Male.PD.PWAS.P, dat2$comb.Male.SD.PWAS.P))
)

# move comb columns
mov_cols = c("comb.Female.Z", "comb.Female.P")
dat2 = dat2 %>% relocate(all_of(mov_cols), .after = "comb.Female.SD.PWAS.P")

## get start and end for each gene. Create variables for using ensembl with human genes for build 38 Ensembl Release 113
mart = useMart("ensembl")
mart = useDataset("hsapiens_gene_ensembl", mart)
ensembl = useEnsembl(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')

#set variable for filter based on input
filters = c("hgnc_symbol")

#set variable for attributes that I would like to receive from the search
attributes = c("hgnc_symbol", "external_gene_name", "entrezgene_accession", "chromosome_name", "start_position","end_position")

genes = dat2$Gene

#create variable to run Ensembl search
g1 = getBM(attributes=attributes, filters=filters, values=genes, mart = ensembl)

#create subset where hgnc, external gene name and entrez gene match
g <- g1 %>%
  filter(hgnc_symbol == external_gene_name & hgnc_symbol == entrezgene_accession
); nrow(g)  #38
length(unique(g$hgnc_symbol))   #36; meaning there is a duplicate

# remove rows that have chromosome_name not 1:22
g <- subset(g, chromosome_name %in% c(1:22))

# add start and end columns to dat2
dat2$START <- g$start_position[match(dat2$Gene, g$hgnc_symbol)]
dat2$END <- g$end_position[match(dat2$Gene, g$hgnc_symbol)]

# This code adds a new column Sensitivity with "Yes" or "No" according to your two-part logic:
#Part 1: for females when FDR and P are significant and the combined Z-score is stronger.
#Part 2: fallback to males when female FDR isn't significant but male values meet the threshold and combined effect is stronger.
dat2 <- dat2 %>%
  mutate(Sensitivity = case_when(
    # Female: FDR < 0.05 and combined P < 0.05
    !is.na(Female.fdr) & Female.fdr < 0.05 &
      !is.na(comb.Female.P) & comb.Female.P < 0.05 &
      !is.na(Female.Z) & !is.na(comb.Female.Z) &
      (
        (Female.Z > 0 & comb.Female.Z > Female.Z) |
        (Female.Z < 0 & comb.Female.Z < Female.Z)
      ) ~ "Yes",

    # Male: Female FDR not significant, Male FDR < 0.05 and combined P < 0.05
    (!is.na(Female.fdr) & Female.fdr >= 0.05 | is.na(Female.fdr)) &
      !is.na(Male.fdr) & Male.fdr < 0.05 &
      !is.na(comb.Male.P) & comb.Male.P < 0.05 &
      !is.na(Male.Z) & !is.na(comb.Male.Z) &
      (
        (Male.Z > 0 & comb.Male.Z > Male.Z) |
        (Male.Z < 0 & comb.Male.Z < Male.Z)
      ) ~ "Yes",

    # All other cases
    TRUE ~ "No"
  ))

fwrite(dat2, paste0(dir, "ADGC_ADSP_UKB_FinnGen_cc_rb.gen090.noAPOE.shared_var.W23nCSF.PWAS_noMHC_ext2Mb_top-sex-specific-genes_AFRad_full-PWAS.txt"), col.names = T, row.names = F, quote = F, na = NA, sep = '\t')
