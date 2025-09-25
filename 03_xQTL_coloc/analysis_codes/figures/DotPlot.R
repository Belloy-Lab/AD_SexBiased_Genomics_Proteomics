#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(haven)
  library(gridExtra)
  library(ggpattern)
  library(UpSetR)
  library(stringr)
  library(biomaRt)
  library(CMplot)
  library(readxl)
  library(ggpubr)
  library(RColorBrewer)
  library(scales)
  library(patchwork)
  library(forcats)
  library(optparse)
})

# Define command line options
option_list = list(
  make_option("--dir", type = "character", default = NULL,
              help = "Directory containing input files", metavar = "character"),
  make_option("--data", type = "character", default = NULL,
              help = "Female top sex-specific genes from brain/CSF PWAS with dotplot matrix ", metavar = "character"),
  make_option("--f_ss_genes", type = "character", default = NULL,
              help = "Female top sex-specific genes from both brain/CSF PWAS", metavar = "character"),
  make_option("--m_ss_genes", type = "character", default = NULL,
              help = "Male top sex-specific genes from both brain/CSF PWAS", metavar = "character")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Assign arguments to variables
dir <- opt$dir
setwd(dir)
# load data to extract PWAS results
f_dat = fread(paste0(dir, opt$data))
m_dat = fread(paste0(dir, opt$data))
f = fread(paste0(dir, opt$f_ss_genes))
m = fread(paste0(dir, opt$f_ss_genes)) 
f_dat= as.data.frame(f_dat); m_dat= as.data.frame(m_dat); m = as.data.frame(m); f = as.data.frame(f)

f_dat$COLOC_PP4 = as.numeric(f_dat$COLOC_PP4); f_dat$SMR_pFDR = as.numeric(f_dat$SMR_pFDR); f_dat$HEIDI_P = as.numeric(f_dat$HEIDI_P)
f_dat$best_pQTL_COLOC_PP4 = as.numeric(f_dat$best_pQTL_COLOC_PP4)
f_dat<- dplyr::select(f_dat, -Novelty, -Female_PWAS.P, -Female_PWAS.Z, -Male_PWAS.P, -Male_PWAS.Z)

# change name of column
names(f_dat)[6] <- "Prim.v.Sec"

# initiate columns to add female data across tissue and discovery
f_dat$Female.B.PD.Z <- NA
f_dat$Female.B.PD.P <- NA
f_dat$Female.B.PD.pFDR <- NA

f_dat$Female.B.SD.Z <- NA
f_dat$Female.B.SD.P <- NA
f_dat$Female.B.SD.pFDR <- NA

f_dat$Female.C.PD.Z <- NA
f_dat$Female.C.PD.P <- NA
f_dat$Female.C.PD.pFDR <- NA

f_dat$Female.C.SD.Z <- NA
f_dat$Female.C.SD.P <- NA
f_dat$Female.C.SD.pFDR <- NA

# Processing top gene sex-spec data to extract PWAS results (male and female)
f$TWAS.Z = as.numeric(f$TWAS.Z); f$TWAS.P = as.numeric(f$TWAS.P); f$CHR = as.numeric(f$CHR)
m$TWAS.Z = as.numeric(m$TWAS.Z); m$TWAS.P = as.numeric(m$TWAS.P); m$CHR = as.numeric(m$CHR)

assign_brain_data <- function(f_dat, f, col_prefix, discovery, tissue) {
  # Filter the f dataframe
  f_filtered <- f %>% filter(Discovery == discovery & Tissue == tissue)
  
  # Merge with the f_datdataframe
  merged_data <- f_dat%>%
    left_join(f_filtered %>% dplyr::select(Gene, TWAS.Z, TWAS.P, fdr_p), by = "Gene")
  
  # Update columns in f_datwith matching gene names
  f_dat[[paste0(col_prefix, "Z")]] <- merged_data$TWAS.Z
  f_dat[[paste0(col_prefix, "P")]] <- merged_data$TWAS.P
  f_dat[[paste0(col_prefix, "pFDR")]] <- merged_data$fdr_p
  
  return(f_dat)
}

# Update f_datwith filtered data for each appropriate column
f_dat<- assign_brain_data(f_dat, f, "Female.B.PD.", "Primary", "B")
f_dat<- assign_brain_data(f_dat, f, "Female.B.SD.", "Secondary", "B")

# Function specifically for processing CSF data
assign_csf_data <- function(f_dat, f, col_prefix, discovery, tissue) {
  # Filter the f dataframe based on discovery and tissue
  f_filtered <- f %>% filter(Discovery == discovery & Tissue == tissue)
  
  # Split f_datdataframe into with and without Analyte scenarios
  dat_with_analyte <- f_dat%>% filter(!is.na(Analyte))
  dat_no_analyte <- f_dat%>% filter(is.na(Analyte))
  
  # Handle cases with Analyte
  if (nrow(dat_with_analyte) > 0) {
    merged_with_analyte <- dat_with_analyte %>%
      left_join(f_filtered %>% dplyr::select(Analyte, Gene, TWAS.Z, TWAS.P, fdr_p), by = c("Analyte", "Gene"))
      
    dat_with_analyte[[paste0(col_prefix, "Z")]] <- merged_with_analyte$TWAS.Z
    dat_with_analyte[[paste0(col_prefix, "P")]] <- merged_with_analyte$TWAS.P
    dat_with_analyte[[paste0(col_prefix, "pFDR")]] <- merged_with_analyte$fdr_p
  }

  # Handle cases without Analyte (only by Gene)
  if (nrow(dat_no_analyte) > 0) {
    merged_no_analyte <- dat_no_analyte %>%
      left_join(f_filtered %>% dplyr::select(Gene, TWAS.Z, TWAS.P, fdr_p), by = "Gene")
      
    dat_no_analyte[[paste0(col_prefix, "Z")]] <- merged_no_analyte$TWAS.Z
    dat_no_analyte[[paste0(col_prefix, "P")]] <- merged_no_analyte$TWAS.P
    dat_no_analyte[[paste0(col_prefix, "pFDR")]] <- merged_no_analyte$fdr_p
  }
  
  # Combine the dataframes back together
  f_dat<- bind_rows(dat_with_analyte, dat_no_analyte)
  
  return(f_dat)
}

# Update f_dat with the appropriate columns
f_dat<- assign_csf_data(f_dat, f, "Female.C.PD.", "Primary", "C")
f_dat<- assign_csf_data(f_dat, f, "Female.C.SD.", "Secondary", "C")

# add info in f_datto indicate which Tissue the finding originated from
fs <- subset(f, topSNP > 0)
ms <- subset(m, topSNP > 0)
fs$Sex = "F"
ms$Sex = "M"
ss <- rbind(fs, ms)

f_dat<- f_dat%>%
    left_join(ss %>% dplyr::select(Analyte, Gene, Tissue, Sex), by = c("Analyte", "Gene")
)

# order f_dat based on Chromosome, Gene_start, Gene_end
f_dat<- f_dat[order(f_dat$Chromosome, f_dat$Gene_start, f_dat$Gene_end), ]

# Add columns to indicate significance labels and SMR support 
f_dat <- f_dat %>%
  mutate(
    Brain_P_SMR = ifelse(Sex == "F" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "B" & Prim.v.Sec %in% c("Primary", "Both"), "black", "grey60"),
    Brain_S_SMR = ifelse(Sex == "F" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "B" & Prim.v.Sec %in% c("Secondary", "Both"), "black", "grey60"),
    
    CSF_P_SMR = ifelse(Sex == "F" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "C" & Prim.v.Sec %in% c("Primary", "Both"), "black", "grey60"),
    CSF_S_SMR = ifelse(Sex == "F" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "C" & Prim.v.Sec %in% c("Secondary", "Both"), "black", "grey60"),

    Brain_P_stroke = ifelse(Sex == "F" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "B" & Prim.v.Sec %in% c("Primary", "Both"), 0.5, 0.1),
    Brain_S_stroke = ifelse(Sex == "F" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "B" & Prim.v.Sec %in% c("Secondary", "Both"), 0.5, 0.1),
    
    CSF_P_stroke = ifelse(Sex == "F" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "C" & Prim.v.Sec %in% c("Primary", "Both"), 0.5, 0.1),
    CSF_S_stroke = ifelse(Sex == "F" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "C" & Prim.v.Sec %in% c("Secondary", "Both"), 0.5, 0.1),

    label_BP = ifelse(Female.B.PD.pFDR < 0.05, "**", ifelse(Female.B.PD.P < 0.05, "*", NA)),
    label_BS = ifelse(Female.B.SD.pFDR < 0.05, "**", ifelse(Female.B.SD.P < 0.05, "*", NA)),
    label_CP = ifelse(Female.C.PD.pFDR < 0.05, "**", ifelse(Female.C.PD.P < 0.05, "*", NA)),
    label_CS = ifelse(Female.C.SD.pFDR < 0.05, "**", ifelse(Female.C.SD.P < 0.05, "*", NA))
)

# Male analysis
m_dat$COLOC_PP4 = as.numeric(m_dat$COLOC_PP4); m_dat$SMR_pFDR = as.numeric(m_dat$SMR_pFDR); m_dat$HEIDI_P = as.numeric(m_dat$HEIDI_P)
m_dat$best_pQTL_COLOC_PP4 = as.numeric(m_dat$best_pQTL_COLOC_PP4)
m_dat<- dplyr::select(m_dat, -Novelty, -Female_PWAS.P, -Female_PWAS.Z, -Male_PWAS.P, -Male_PWAS.Z)

# change name of column
names(m_dat)[6] <- "Prim.v.Sec"

# initiate columns to add female data across tissue and discovery
m_dat$Male.B.PD.Z <- NA
m_dat$Male.B.PD.P <- NA
m_dat$Male.B.PD.pFDR <- NA

m_dat$Male.B.SD.Z <- NA
m_dat$Male.B.SD.P <- NA
m_dat$Male.B.SD.pFDR <- NA

m_dat$Male.C.PD.Z <- NA
m_dat$Male.C.PD.P <- NA
m_dat$Male.C.PD.pFDR <- NA

m_dat$Male.C.SD.Z <- NA
m_dat$Male.C.SD.P <- NA
m_dat$Male.C.SD.pFDR <- NA

assign_brain_data <- function(m_dat, m, col_prefix, discovery, tissue) {
  # Filter the m dataframe
  m_filtered <- m %>% filter(Discovery == discovery & Tissue == tissue)
  
  # Merge with the m_dat dataframe
  merged_data <- m_dat%>%
    left_join(m_filtered %>% dplyr::select(Gene, TWAS.Z, TWAS.P, fdr_p), by = "Gene")
  
  # Update columns in m_dat with matching gene names
  m_dat[[paste0(col_prefix, "Z")]] <- merged_data$TWAS.Z
  m_dat[[paste0(col_prefix, "P")]] <- merged_data$TWAS.P
  m_dat[[paste0(col_prefix, "pFDR")]] <- merged_data$fdr_p
  
  return(m_dat)
}

# Update m_dat with filtered data for each appropriate column
m_dat<- assign_brain_data(m_dat, m, "Male.B.PD.", "Primary", "B")
m_dat<- assign_brain_data(m_dat, m, "Male.B.SD.", "Secondary", "B")

# Function specifically for processing CSF data
assign_csf_data <- function(m_dat, m, col_prefix, discovery, tissue) {
  # Filter the m dataframe based on discovery and tissue
  m_filtered <- m %>% filter(Discovery == discovery & Tissue == tissue)
  
  # Split m_dat dataframe into with and without Analyte scenarios
  dat_with_analyte <- m_dat%>% filter(!is.na(Analyte))
  dat_no_analyte <- m_dat%>% filter(is.na(Analyte))
  
  # Handle cases with Analyte
  if (nrow(dat_with_analyte) > 0) {
    merged_with_analyte <- dat_with_analyte %>%
      left_join(m_filtered %>% dplyr::select(Analyte, Gene, TWAS.Z, TWAS.P, fdr_p), by = c("Analyte", "Gene"))
      
    dat_with_analyte[[paste0(col_prefix, "Z")]] <- merged_with_analyte$TWAS.Z
    dat_with_analyte[[paste0(col_prefix, "P")]] <- merged_with_analyte$TWAS.P
    dat_with_analyte[[paste0(col_prefix, "pFDR")]] <- merged_with_analyte$fdr_p
  }

  # Handle cases without Analyte (only by Gene)
  if (nrow(dat_no_analyte) > 0) {
    merged_no_analyte <- dat_no_analyte %>%
      left_join(m_filtered %>% dplyr::select(Gene, TWAS.Z, TWAS.P, fdr_p), by = "Gene")
      
    dat_no_analyte[[paste0(col_prefix, "Z")]] <- merged_no_analyte$TWAS.Z
    dat_no_analyte[[paste0(col_prefix, "P")]] <- merged_no_analyte$TWAS.P
    dat_no_analyte[[paste0(col_prefix, "pFDR")]] <- merged_no_analyte$fdr_p
  }
  
  # Combine the dataframes back together
  m_dat<- bind_rows(dat_with_analyte, dat_no_analyte)
  
  return(m_dat)
}

# Update m_dat with the appropriate columns
m_dat<- assign_csf_data(m_dat, m, "Male.C.PD.", "Primary", "C")
m_dat<- assign_csf_data(m_dat, m, "Male.C.SD.", "Secondary", "C")

# add info in m_dat to indicate which Tissue the finding originated from
fs <- subset(f, topSNP > 0)
ms <- subset(m, topSNP > 0)
fs$Sex = "F"
ms$Sex = "M"
ss <- rbind(fs, ms)

m_dat<- m_dat%>%
    left_join(ss %>% dplyr::select(Analyte, Gene, Tissue, Sex), by = c("Analyte", "Gene")
)

# order m_dat based on Chromosome, Gene_start, Gene_end
m_dat<- m_dat[order(m_dat$Chromosome, m_dat$Gene_start, m_dat$Gene_end), ]

# Add columns to indicate significance labels and SMR support 
m_dat <- m_dat %>%
  mutate(
    Brain_P_SMR = ifelse(Sex == "M" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "B" & Prim.v.Sec %in% c("Primary", "Both"), "black", "grey60"),
    Brain_S_SMR = ifelse(Sex == "M" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "B" & Prim.v.Sec %in% c("Secondary", "Both"), "black", "grey60"),
    
    CSF_P_SMR = ifelse(Sex == "M" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "C" & Prim.v.Sec %in% c("Primary", "Both"), "black", "grey60"),
    CSF_S_SMR = ifelse(Sex == "M" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "C" & Prim.v.Sec %in% c("Secondary", "Both"), "black", "grey60"),

    Brain_P_stroke = ifelse(Sex == "M" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "B" & Prim.v.Sec %in% c("Primary", "Both"), 0.5, 0.1),
    Brain_S_stroke = ifelse(Sex == "M" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "B" & Prim.v.Sec %in% c("Secondary", "Both"), 0.5, 0.1),
    
    CSF_P_stroke = ifelse(Sex == "M" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "C" & Prim.v.Sec %in% c("Primary", "Both"), 0.5, 0.1),
    CSF_S_stroke = ifelse(Sex == "M" & SMR_pFDR < 0.05 & (HEIDI_P > 0.05 | is.na(HEIDI_P)) & Tissue == "C" & Prim.v.Sec %in% c("Secondary", "Both"), 0.5, 0.1),

    label_BP = ifelse(Male.B.PD.pFDR < 0.05, "**", ifelse(Male.B.PD.P < 0.05, "*", "")),
    label_BS = ifelse(Male.B.SD.pFDR < 0.05, "**", ifelse(Male.B.SD.P < 0.05, "*", "")),
    label_CP = ifelse(Male.C.PD.pFDR < 0.05, "**", ifelse(Male.C.PD.P < 0.05, "*", "")),
    label_CS = ifelse(Male.C.SD.pFDR < 0.05, "**", ifelse(Male.C.SD.P < 0.05, "*", ""))
)

# sort on CHR/BP and by sex making males on bottom
f_dat<- f_dat[order(f_dat$Sex, f_dat$Chromosome, f_dat$Gene_start, f_dat$Gene_end), ]
m_dat<- m_dat[order(m_dat$Sex, m_dat$Chromosome, m_dat$Gene_start, m_dat$Gene_end), ]

fwrite(f_dat, paste0(dir, "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_top-fsg_dotplot_matrix11.txt"), col.names = T, row.names = F, quote = F, sep = '\t', na = NA)
fwrite(m_dat, paste0(dir, "ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_top-msg_dotplot_matrix11.txt"), col.names = T, row.names = F, quote = F, sep = '\t', na = NA)

gene_order <- f_dat$Gene  # assuming f_dat and m_dat are in the same order

f_dat$Gene <- factor(f_dat$Gene, levels = rev(gene_order))
m_dat$Gene <- factor(m_dat$Gene, levels = rev(gene_order))

# determine the range of z-scores between male and female data for color scale
max_z = ifelse(max(c(f_dat$Female.B.PD.Z, f_dat$Female.B.SD.Z, f_dat$Female.C.PD.Z, f_dat$Female.C.SD.Z), na.rm = TRUE) > max(c(m_dat$Male.B.PD.Z, m_dat$Male.B.SD.Z, m_dat$Male.C.PD.Z, m_dat$Male.C.SD.Z), na.rm = TRUE),
                max(c(f_dat$Female.B.PD.Z, f_dat$Female.B.SD.Z, f_dat$Female.C.PD.Z, f_dat$Female.C.SD.Z), na.rm = TRUE),
                max(c(m_dat$Male.B.PD.Z, m_dat$Male.B.SD.Z, m_dat$Male.C.PD.Z, m_dat$Male.C.SD.Z), na.rm = TRUE)
)
max_z <- ceiling(max_z) + 1
min_z <- -(max_z)

# determine smallest abs(z-score) between male and female for size lower limit
min_abs_z_female <- min(abs(c(
  f_dat$Female.B.PD.Z,
  f_dat$Female.B.SD.Z,
  f_dat$Female.C.PD.Z,
  f_dat$Female.C.SD.Z
)), na.rm = TRUE)

min_abs_z_male <- min(abs(c(
  m_dat$Male.B.PD.Z,
  m_dat$Male.B.SD.Z,
  m_dat$Male.C.PD.Z,
  m_dat$Male.C.SD.Z
)), na.rm = TRUE)

min_size = ifelse(min_abs_z_female < min_abs_z_male, min_abs_z_female, min_abs_z_male)

max_size = ifelse(max(abs(c(
  f_dat$Female.B.PD.Z,
  f_dat$Female.B.SD.Z,
  f_dat$Female.C.PD.Z,
  f_dat$Female.C.SD.Z
)), na.rm = TRUE) > max(abs(c(
  m_dat$Male.B.PD.Z,
  m_dat$Male.B.SD.Z,
  m_dat$Male.C.PD.Z,
  m_dat$Male.C.SD.Z
)), na.rm = TRUE), max(abs(c(
  f_dat$Female.B.PD.Z,
  f_dat$Female.B.SD.Z,
  f_dat$Female.C.PD.Z,
  f_dat$Female.C.SD.Z
)), na.rm = TRUE), max(abs(c(
  m_dat$Male.B.PD.Z,
  m_dat$Male.B.SD.Z,
  m_dat$Male.C.PD.Z,
  m_dat$Male.C.SD.Z
)), na.rm = TRUE))

max_size <- ceiling(max_size)

f_dat <- f_dat %>%
    dplyr::mutate(COLOC_color = ifelse(f_dat$Sex == "F", "#700000", "#3c7899")
)


## smaller plot
f_dot <- ggplot() +
    geom_point(data = f_dat, aes(x = "Brain Primary", y = Gene, size = abs(Female.B.PD.Z), fill = Female.B.PD.Z), shape = 21, color = f_dat$Brain_P_SMR, stroke = f_dat$Brain_P_stroke * 0.6) +
    geom_text(data = f_dat, aes(x = "Brain Primary", y = Gene, label = label_BP), vjust = 0.8, size = 2.5) +

    geom_point(data = f_dat, aes(x = "Brain Secondary", y = Gene, size = abs(Female.B.SD.Z), fill = Female.B.SD.Z), shape = 21, color = f_dat$Brain_S_SMR, stroke = f_dat$Brain_S_stroke * 0.6) +
    geom_text(data = f_dat, aes(x = "Brain Secondary", y = Gene, label = label_BS), vjust = 0.8, size = 2.5) +

    geom_point(data = f_dat, aes(x = "CSF Primary", y = Gene, size = abs(Female.C.PD.Z), fill = Female.C.PD.Z), shape = 21, color = f_dat$CSF_P_SMR, stroke = f_dat$CSF_P_stroke * 0.6) +
    geom_text(data = f_dat, aes(x = "CSF Primary", y = Gene, label = label_CP), vjust = 0.8, size = 2.5) +

    geom_point(data = f_dat, aes(x = "CSF Secondary", y = Gene, size = abs(Female.C.SD.Z), fill = Female.C.SD.Z), shape = 21, color = f_dat$CSF_S_SMR, stroke = f_dat$CSF_S_stroke * 0.6) +
    geom_text(data = f_dat, aes(x = "CSF Secondary", y = Gene, label = label_CS), vjust = 0.8, size = 2.5) +

    scale_size_continuous(
        limits = c(min_size, max_size),
        breaks = seq(1, max_size, by = 1),
        range = c(1.2, 3.8),
        name = NULL
    ) +
    scale_fill_gradient2(
        low = "#8F65BA", high = "#B8E868",
        limits = c(min_z, max_z),
        breaks = seq(min_z, max_z, by = 1),
        labels = seq(min_z, max_z, by = 1),
        name = NULL
    ) +
    theme_bw() +
    theme(
        plot.margin = margin(0, 0, 0, 0, "points"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.ticks.y.left = element_blank(),
        axis.text.x = element_text(size = 6, angle = 35, hjust = 0.9),
        axis.title.x = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        legend.text = element_text(size = 5.36, color = "grey50"),
        legend.title = element_text(size = 5.36, color = "grey50", margin = margin(b  = -6)),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.spacing.y = unit(0.1, 'mm'),
        legend.box.margin = margin(0, 0, -40, 0, unit = "pt"),
        legend.margin = margin(0,0,0,0),
        panel.grid = element_line(linewidth = 0.1)
    ) +
    guides(
        size = guide_legend(
            title = NULL,
            direction = "horizontal",
            label.position = NULL,
            override.aes = list(shape = 1),
            nrow = 1,
            order = 1
        ),
        fill = guide_colorbar(
            title = NULL,
            barwidth = 8,
            barheight = 0.67,
            title.position = "left",
            direction = "horizontal",
            label.position  = "bottom",
            order = 2
        )
)

m_dot <- ggplot() +
    geom_point(data = m_dat, aes(x = "Brain Primary", y = Gene, size = abs(Male.B.PD.Z), fill = Male.B.PD.Z), shape = 21, color = m_dat$Brain_P_SMR, stroke = m_dat$Brain_P_stroke * 0.6, show.legend = FALSE) +
    geom_text(data = m_dat, aes(x = "Brain Primary", y = Gene, label = label_BP), vjust = 0.8, size = 2.5, show.legend = FALSE) +

    geom_point(data = m_dat, aes(x = "Brain Secondary", y = Gene, size = abs(Male.B.SD.Z), fill = Male.B.SD.Z), shape = 21, color = m_dat$Brain_S_SMR, stroke = m_dat$Brain_S_stroke * 0.6, show.legend = FALSE) +
    geom_text(data = m_dat, aes(x = "Brain Secondary", y = Gene, label = label_BS), vjust = 0.8, size = 2.5, show.legend = FALSE) +

    geom_point(data = m_dat, aes(x = "CSF Primary", y = Gene, size = abs(Male.C.PD.Z), fill = Male.C.PD.Z), shape = 21, color = m_dat$CSF_P_SMR, stroke = m_dat$CSF_P_stroke * 0.6, show.legend = FALSE) +
    geom_text(data = m_dat, aes(x = "CSF Primary", y = Gene, label = label_CP), vjust = 0.8, size = 2.5, show.legend = FALSE) +

    geom_point(data = m_dat, aes(x = "CSF Secondary", y = Gene, size = abs(Male.C.SD.Z), fill = Male.C.SD.Z), shape = 21, color = m_dat$CSF_S_SMR, stroke = m_dat$CSF_S_stroke * 0.6, show.legend = FALSE) +
    geom_text(data = m_dat, aes(x = "CSF Secondary", y = Gene, label = label_CS), vjust = 0.8, size = 2.5, show.legend = FALSE) +

    scale_size_continuous(
        limits = c(min_size, max_size),
        breaks = seq(1, max_size, by = 1),
        range = c(1.2, 3.8),
        name = NULL
    ) +
    scale_fill_gradient2(
        low = "#8F65BA", high = "#B8E868",
        limits = c(min_z, max_z),
        breaks = seq(min_z, max_z, by = 1),
        labels = seq(min_z, max_z, by = 1),
        name = NULL
    ) +
    theme_bw() +
    theme(
        plot.margin = margin(0, -5, 0, 0, "points"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.x = element_text(size = 6, angle = 35, hjust = 0.9),
        axis.title.x = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        legend.position = "none",
        panel.grid = element_line(linewidth = 0.1)
)

## connected scatterplot
coloc <- ggplot(f_dat, aes(x = best_pQTL_COLOC_PP4, y = Gene)) +
  geom_path(group = 1, color = "black", linewidth = 0.15) +
  geom_point(aes(x = best_pQTL_COLOC_PP4), shape = 21, fill = f_dat$COLOC_color, stroke = 0, size = 2.75) +
  scale_x_continuous(
    limits = c(0, 1), 
    breaks = c(0, 0.4, 0.7, 1),
    labels = function(x) ifelse(x == 0, "0", format(x, nsmall = 1))) +
  geom_vline(xintercept = 0.4, color = "black", linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = 0.7, color = "black", linewidth = 0.3, linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.x.bottom = element_line(linewidth = 0.2),
    axis.line.x.bottom = element_line(linewidth = 0.2),
    plot.margin = margin(0, 0, 0, -2, "points")
)

locus_dat <- fread(paste0(dir, "top_PP4_per_PWAS_locus.txt"), header = T, sep = '\t'); locus_dat = as.data.frame(locus_dat)

locus_dat <- locus_dat %>% 
    mutate(Gene = locus_label) %>%
    separate_rows(locus_label, sep = '/'
)

locus_dat = as.data.frame(locus_dat)

colnames(locus_dat) = c("Gene", "qtl_type", "best_PP4", "locus_label")

locus_dat <- locus_dat %>%
  dplyr::distinct(Gene, .keep_all = TRUE
)

# update locus label for ENO3 to ENO3/RABEP1, then duplicate the row and change the Gene label from ENO3 to RABEP1
locus_dat <- locus_dat %>%
  mutate(locus_label = ifelse(locus_label == "ENO3", "ENO3/RABEP1", locus_label)
)
duplicated_row <- locus_dat %>%
  filter(locus_label == "ENO3/RABEP1" & Gene == "ENO3") %>%
  mutate(Gene = "RABEP1"
)
locus_dat <- bind_rows(locus_dat, duplicated_row)

locus_dat$Gene <- factor(locus_dat$Gene, levels = gene_order)
locus_dat <- locus_dat[order(locus_dat$Gene), ]
locus_dat$Gene <- factor(locus_dat$Gene, levels = rev(gene_order))

# Map gene to numeric y-position (from factor level order)
gene_ypos <- data.frame(
  Gene = rev(gene_order),
  y = seq_along(gene_order)
)

# Join y-positions with locus_dat
locus_center <- locus_dat %>%
left_join(gene_ypos, by = "Gene") %>%
  group_by(locus_label) %>%
  mutate(
    center_y = mean(y),
    center_PP4 = mean(best_PP4)
    ) %>%
    ungroup(
)

locus_center = as.data.frame(locus_center)

coloc <- coloc + 
    geom_path(data = locus_center, aes(x = best_PP4, y = center_y, group = 1),
                color = "black", linetype = "dashed", linewidth = 0.15, inherit.aes = FALSE) +
    geom_point(data = locus_center, aes(x = best_PP4, y = center_y),
                shape = 23, fill = "grey30", stroke = 0, size = 2.75, inherit.aes = FALSE
)
p <- f_dot + m_dot + coloc + plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "bottom", legend.box = "vertical") & coord_cartesian(clip = "off")
pdf(paste0(dir, "ADGC_ADSP_UKB_FinnGen_cc_rb.gen090.noAPOE.shared_var.ssg_W23nCSFc_dotplot_matrix11.pdf"), width = 4, height = 6)
plot(p)
dev.off()
