library(data.table)
library(plyr)
library(dplyr)
library(locuszoomr)
library(AnnotationHub)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(RColorBrewer)
library(rtracklayer)
library(patchwork)
library(grid)

# souece the function script
source("Analysis_codes/functions_overlay_plot.r")

female_90 <- fread("/Path/to/Female/GWAS/summstat/file")

male_90 <- fread("/Path/to/Male/GWAS/summstat/file")

w23_pQTL_female <- fread(paste0("/Path/to/Female/Brain/pQTL/file"))

w23_pQTL_male <- fread(paste0("/Path/to/Male/Brain/pQTL/file"))

w23_pQTL_combined <- fread(paste0("/Path/to/Combined/Brain/pQTL/file"))

brain_meta_DLPFC <- fread("/Path/to/DLPFC/brain/meta/mQTL/file")

wingo_female_pQTL <- fread("/Path/to/DLPFC/female/brain/pQTL/data.txt")

eQTLgen <- fread("/Path/to/blood/eQTL/blood_eQTL.csv")

eQTL_catalogue_nedelec_naive <- fread("/Path/to/nedelec_naive/eQTL/data.tsv.gz")

gtf <- rtracklayer::import("/Path/to/gencode.v47.basic.annotation.gtf.gz")
ref_seq <- as.data.table(gtf)

gene_f <- fread(paste0("/Path/to/Female/prioratized/genes/from/Brain/CSF/PWAA/results/list.txt"))

#define output directory
out_dir <- paste0("/Path/to/store/output/Overlay_plots/")
fold_w <- paste0("/Path/to/store/temporary/LD_temp_save/")

file_snp <- "_common_id.txt"
pl <- "plink1.9"
in_fold <- "/Path/to/EU_TOPMed/data/in/PLINK/format/separeted/byCHR/"
in_file <- "TOPMed_ch"
keep_file <- "/Path/to/EU_nHISP_TOPMed_list_Unrel.txt"

esngid <- esngid <- sub("\\..*", "", ref_seq[ref_seq$gene_name == "PSEN1", "gene_id"][1])
w23_pQTL_female_gene <- w23_pQTL_female %>% 
    dplyr::filter(gene_id == esngid)

esngid <- esngid <- sub("\\..*", "", ref_seq[ref_seq$gene_name == "PAPLN", "gene_id"][1])
wingo_female_pQTL_gene <- wingo_female_pQTL  %>% 
    dplyr::filter(gene_id == esngid) %>%
    dplyr::mutate(P = pvalue)

esngid <- esngid <- sub("\\..*", "", ref_seq[ref_seq$gene_name == "HEATR4", "gene_id"][1])
eQTLgen_gene <- eQTLgen %>% 
    dplyr::filter(gene_id == esngid) %>% 
    dplyr::mutate(P = pvalue)

brain_meta_DLPFC_gene <- brain_meta_DLPFC %>% 
    dplyr::filter(Probe == "cg24036523") %>% 
    dplyr::mutate(P = p, A1 = ALLELE1, A2 = ALLELE0)


esngid <- esngid <- sub("\\..*", "", ref_seq[ref_seq$gene_name == "PSEN1", "gene_id"][1])
eQTL_catalogue_nedelec_naive_gene <- eQTL_catalogue_nedelec_naive %>% 
    dplyr::filter(gene_id == esngid) %>% 
    dplyr::mutate(P = pvalue)

rsid <- "rs11159021"
rsid_2 <- "rs3742825"
rsid_3 <- "rs74986264" 
chr <- 14
bp <- 73238768

merged_data_list <- select_data_merge_LD(rsid,rsid_2,rsid_3,chr,bp,female_90,male_90,
                                        w23_pQTL_female_gene,
                                        wingo_female_pQTL_gene,
                                        eQTLgen_gene,
                                        brain_meta_DLPFC_gene,
                                        eQTL_catalogue_nedelec_naive_gene,
                                        fold_w,file_snp,in_file,pl,in_fold,keep_file,7e5)
merged_female_withld <- merged_data_list[["female"]]
merged_male_withld <- merged_data_list[["male"]]
merged_w23_pQTL_female_gene_with_ld <- merged_data_list[["1"]]
merged_wingo_female_pQTL_gene_with_ld <- merged_data_list[["2"]]
merged_eQTLgen_gene_with_ld <- merged_data_list[["3"]]
merged_brain_meta_DLPFC_gene_with_ld <- merged_data_list[["4"]]
merged_eQTL_catalogue_nedelec_naive_gene_with_ld <- merged_data_list[["5"]]

max_ylim_female <- find_max_ylim_AD(merged_female_withld, rsid)
max_ylim_male <- find_max_ylim_AD(merged_male_withld, rsid)
max_ylim <- max(max_ylim_female, max_ylim_male)
# margin_y 1 7 13 
f <- get_ggplot(max_ylim,bp,merged_female_withld,merged_female_withld,rsid,rsid_2,rsid_3,"AD Female GWAS",13,out_dir)
m <- get_ggplot(max_ylim,bp,merged_female_withld,merged_male_withld,rsid,rsid_2,rsid_3,"AD Male GWAS",13,out_dir)
q1 <- get_ggplot(max_ylim,bp,merged_female_withld,merged_w23_pQTL_female_gene_with_ld,rsid,rsid_2,rsid_3,"w23_pQTL_female",13,out_dir)
q2 <- get_ggplot(max_ylim,bp,merged_female_withld,merged_wingo_female_pQTL_gene_with_ld,rsid,rsid_2,rsid_3,"wingo_female_DLPFC_pQTL_Female_PSEN1_PAPLN",13,out_dir)
q3 <- get_ggplot(max_ylim,bp,merged_female_withld,merged_eQTLgen_gene_with_ld,rsid,rsid_2,rsid_3,"eQTLgen_blood_eQTL_Female_PSEN1_HEATR4",7,out_dir)
q4 <- get_ggplot(max_ylim,bp,merged_female_withld,merged_brain_meta_DLPFC_gene_with_ld,rsid,rsid_2,rsid_3,"BrainMeta_DLPFC_mQTL_Female_PSEN1_cg24036523",7,out_dir)
q5 <- get_ggplot(max_ylim,bp,merged_female_withld,merged_eQTL_catalogue_nedelec_naive_gene_with_ld,rsid,rsid_2,rsid_3,"eQTL_catalogue_nedelec_naive_eQTL_Female_PSEN1_PSEN1",13,out_dir)
# margin_y 0 2 6
c1 <- get_ggplot_compare(merged_female_withld,merged_w23_pQTL_female_gene_with_ld,"middle",rsid_2,"w23_pQTL_female",
                        FALSE,11.5,c(0,2,4,6,8),c("0","2","4","6","8"),out_dir)
c2 <- get_ggplot_compare(merged_female_withld,merged_wingo_female_pQTL_gene_with_ld,"left",rsid_3,"wingo_female_DLPFC",
                        FALSE,11.5,c(0,1,2,3,4,5,6),c("0","1","2","3","4","5","6"),out_dir)
c3 <- get_ggplot_compare(merged_female_withld,merged_eQTLgen_gene_with_ld,"right",rsid,"eQTLgen_blood_eQTL",
                        FALSE,6,c(0,5,10,15,20),c("0","5","10","15","20"),out_dir)
c4 <- get_ggplot_compare(merged_female_withld,merged_brain_meta_DLPFC_gene_with_ld,"middle",rsid_2,"BrainMeta_DLPFC_mQTL",
                        FALSE,6,c(2,4,6,8,10),c("2","4","6","8","10"),out_dir)
c5 <- get_ggplot_compare(merged_female_withld,merged_eQTL_catalogue_nedelec_naive_gene_with_ld,"middle",rsid_2,"eQTL_catalogue_nedelec_naive_eQTL",
                        FALSE,11.5,c(0,1,2,3,4,5,6),c("0","1","2","3","4","5","6"),out_dir)
maxrows_num <- 15
g <- get_ggplot_genetrack(maxrows_num,bp,merged_female_withld,merged_female_withld,rsid,rsid_2,rsid_3,out_dir,"PSEN1","PAPLN","HEATR4","cg24036523")
























































