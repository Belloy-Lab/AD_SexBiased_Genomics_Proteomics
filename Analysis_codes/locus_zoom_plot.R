#package requirements
library(locuszoomr)
library(data.table)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(LDlinkR)
library(ggplot2)
library(locuscomparer)
library(htmlwidgets)
library(readxl)
library(patchwork)
library(cowplot)
source("Analysis_codes/functions_locus_zoom_plot.R")

#data requirements
##AD GWAS data
AD_GWAS_female_EUROFINN_meta <- fread(paste0("/Path/to/Female/GWAS/summstat/file"))

AD_GWAS_male_EUROFINN_meta <-fread(paste0("/Path/to/Male/GWAS/summstat/file"))

##wingo pQTL data
w23_pQTL_female <- fread(paste0("/Path/to/Female/Brain/pQTL/file"))

w23_pQTL_male <- fread(paste0("/Path/to/Male/Brain/pQTL/file"))

w23_pQTL_combined <- fread(paste0("/Path/to/combined/Brain/pQTL/file"))

# Set path to store outputs
out_dir <- "/Path/to/store/result/plots/"

##CSF data
csf_data <- read_in_csf("X3642.4")
csf_female <- csf_data[["csf_female"]]
csf_male <- csf_data[["csf_male"]]
csf_combined <- csf_data[["csf_combined"]]

#data preprocessing
AD_GWAS_male_EUROFINN_meta <- as.data.frame(AD_GWAS_male_EUROFINN_meta)

AD_GWAS_female_EUROFINN_meta <- AD_GWAS_female_EUROFINN_meta %>%
  dplyr::rename(A1 = ALLELE1, A2 = ALLELE0, MAF = A1FREQ) %>%
  dplyr::mutate(
    id12 = paste(CHR, BP, A1, A2, sep = ":"),
    id21 = paste(CHR, BP, A2, A1, sep = ":")
  )

AD_GWAS_male_EUROFINN_meta <- AD_GWAS_male_EUROFINN_meta %>%
  dplyr::rename(A1 = `ALLELE1`, A2 = `ALLELE0`, MAF = `A1FREQ`) %>%
  dplyr::mutate(
    id12 = paste(CHR, BP, A1, A2, sep = ":"),
    id21 = paste(CHR, BP, A2, A1, sep = ":")
  )

w23_pQTL_female <- as.data.frame(w23_pQTL_female)
w23_pQTL_male <- as.data.frame(w23_pQTL_male)
w23_pQTL_combined <- as.data.frame(w23_pQTL_combined)

w23_pQTL_female <- w23_pQTL_female %>%
  dplyr::mutate(
    id12 = paste(CHR, BP, A1, A2, sep = ":"),
    id21 = paste(CHR, BP, A2, A1, sep = ":")
  )

w23_pQTL_male <- w23_pQTL_male %>%
  dplyr::mutate(
    id12 = paste(CHR, BP, A1, A2, sep = ":"),
    id21 = paste(CHR, BP, A2, A1, sep = ":")
  )

w23_pQTL_combined <- w23_pQTL_combined %>% 
  dplyr::mutate(
    id12 = paste(CHR, BP, A1, A2, sep = ":"),
    id21 = paste(CHR, BP, A2, A1, sep = ":")
  )

##wingo pQTL data
data_list <- list() 
data_list[["AD_female"]] <- AD_GWAS_female_EUROFINN_meta
data_list[["AD_male"]] <- AD_GWAS_male_EUROFINN_meta
data_list[["pQTL_female"]] <- w23_pQTL_female
data_list[["pQTL_male"]] <- w23_pQTL_male
data_list[["pQTL_combined"]] <- w23_pQTL_combined
## CSF data
data_list <- list() 
data_list[["AD_female"]] <- AD_GWAS_female_EUROFINN_meta
data_list[["AD_male"]] <- AD_GWAS_male_EUROFINN_meta
data_list[["pQTL_female"]] <- csf_female
data_list[["pQTL_male"]] <- csf_male
data_list[["pQTL_combined"]] <- csf_combined

source("Analysis_codes/functions_locus_zoom_plot.R")
#input requirements
##wingo example 
#   Protein CH Gene start  Gene end Novelty status Discovery Tissue Z-score
# 1   ACSL6  5  131949973 132012243     Novel gene   Primary  Brain    3.58
#   P-value  FDR-P PWAS_Sex            ENSG Start_hg38  End_hg38  Gene
# 1 0.00034 0.0206        F ENSG00000164398  131313301 132512200 ACSL6
# gene_temp <- table2_updated[i,]
# bp_start <- gene_temp$Start_hg38
# bp_end <- gene_temp$End_hg38
# chr <- gene_temp$CH
# Gene_id <- gene_temp$ENSG
# Gene_name1 <- gene_temp$Gene
# Gene_name2 <- gene_temp$Protein
# level <- gene_temp$Discovery
#----example--
bp_start <- 131549973
bp_end <- 132512243
chr <- 5
Gene_id <- "ENSG00000164398"
Gene_name1 <- "ACSL6"
Gene_name2 <- "ACSL6"
level <- "Primary"
# defalut position <- "topright" other position <- "bottomright" "topleft"
position <- "topleft"
#max rows to show on the gene track
maxrows_num <- 10
#default AD_sex_code = "AD_female" pQTL_sex_code = "pQTL_female"
AD_sex_code <- "AD_female"
pQTL_sex_code <- "pQTL_female"

index_num <- 1
if (level == "Primary"){
    stack_plot = S14_plots(w23_pQTL_female,w23_pQTL_male,w23_pQTL_combined,
                    AD_GWAS_female_EUROFINN_meta,AD_GWAS_male_EUROFINN_meta,
                    data_list,AD_sex_code,pQTL_sex_code,bp_start,bp_end,chr,Gene_id,Gene_name1,Gene_name2,1e6,index_num,level,maxrows_num,position,ld_table,out_dir)
} else if (level == "Secondary"){
    stack_plot = S14_plots(w23_pQTL_female,w23_pQTL_male,w23_pQTL_combined,
                    AD_GWAS_female_EUROFINN_meta,AD_GWAS_male_EUROFINN_meta,
                    data_list,AD_sex_code,"pQTL_combined",bp_start,bp_end,chr,Gene_id,Gene_name1,Gene_name2,1e6,index_num,level,maxrows_num,position,ld_table,out_dir)
}

##CSF example
#    ENSG_ID   ENSG Analyte updated_Gene   Gene Tissue   CHR Start_hg19 End_hg19
#     <char> <char>  <char>       <char> <char> <char> <int>      <int>    <int>
# 1:    <NA>   <NA> X3642.4         CD84   CD84    CSF     1         NA       NA
#    Start_hg38  End_hg38 hg19_range hg38_range PWAS_Discovery
#         <int>     <int>      <int>      <int>         <char>
# 1:  159526132 161525673         NA    1999541   Secondary
# gene_temp <- gene_f_X3642.4[1,]
# bp_start <- gene_temp$Start_hg38
# bp_end <- gene_temp$End_hg38
# chr <- gene_temp$CHR
# Gene_name1 <- gene_temp$Gene
# Gene_name2 <- gene_temp$updated_Gene
# Gene_id <- gene_temp$Analyte
# level <- gene_temp$PWAS_Discovery
#----example--
bp_start <- 159526132
bp_end <- 161525673
chr <- 1
Gene_id <- "X3642.4"
Gene_name1 <- "CD84"
Gene_name2 <- "CD84"
level <- "Secondary"
# defalut position <- "topright" other position <- "bottomright" "topleft"
position <- "topleft"
#max rows to show on the gene track
maxrows_num <- 6
#default AD_sex_code = "AD_female" pQTL_sex_code = "pQTL_female"
AD_sex_code <- "AD_female"
pQTL_sex_code <- "pQTL_female"

index_num <- 2
if (level == "Primary"){
    stack_plot = S14_plots(csf_female,csf_male,csf_combined,
                    AD_GWAS_female_EUROFINN_meta,AD_GWAS_male_EUROFINN_meta,
                    data_list,AD_sex_code,pQTL_sex_code,bp_start,bp_end,chr,Gene_id,Gene_name1,Gene_name2,2e6,index_num,level,maxrows_num,position,ld_table,out_dir)
} else if (level == "Secondary"){
    stack_plot = S14_plots(csf_female,csf_male,csf_combined,
                    AD_GWAS_female_EUROFINN_meta,AD_GWAS_male_EUROFINN_meta,
                    data_list,AD_sex_code,pQTL_sex_code,bp_start,bp_end,chr,Gene_id,Gene_name1,Gene_name2,2e6,index_num,level,maxrows_num,position,ld_table,out_dir)
}
