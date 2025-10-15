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

csf_data <- read_in_csf("X3642.4")
csf_combined_CD84 <- csf_data[["csf_combined"]]

eQTLgen <- fread("/Path/to/blood/eQTL/blood_eQTL.csv")

kosoy_microglia <- fread("/Path/to/kosoy/eQTL/microglia/microglia_eQTL.txt")

csf_data <- read_in_csf("X3311.27")
csf_combined_FCGR3B <- csf_data[["csf_combined"]]

wingo_both_DLPFC <- fread("/Path/to/non-strat/wingo/pQTL/both_DLPFC/both_DLPFC_pQTL.txt")

eQTL_catalogue_cedar_monocytes <- fread("/Path/to/eQTL/cedar_monocytes/cedar_monocytes_eQTL.tsv.gz")

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

esngid <- "X3642.4"
csf_combined_CD84_gene <- csf_combined_CD84 %>% 
    dplyr::filter(gene_id == esngid)

esngid <- esngid <- sub("\\..*", "", ref_seq[ref_seq$gene_name == "NDUFS2", "gene_id"][1])
eQTLgen_gene <- eQTLgen  %>% 
    dplyr::filter(gene_id == esngid) %>%
    dplyr::mutate(P = pvalue)

esngid <- esngid <- sub("\\..*", "", ref_seq[ref_seq$gene_name == "FCER1G", "gene_id"][1])
kosoy_microglia_gene <- kosoy_microglia %>% 
    dplyr::filter(gene_id == esngid) %>% 
    dplyr::mutate(P = pvalue)

esngid <- "X3311.27"
csf_combined_FCGR3B_gene <- csf_combined_FCGR3B %>% 
    dplyr::filter(gene_id == esngid)

esngid <- esngid <- sub("\\..*", "", ref_seq[ref_seq$gene_name == "FCER1G", "gene_id"][1])
eQTL_catalogue_cedar_monocytes_gene <- eQTL_catalogue_cedar_monocytes %>% 
    dplyr::filter(gene_id == esngid) %>% 
    dplyr::mutate(P = pvalue)

rsid <- "rs2070902"
rsid_2 <- "rs4379692"
rsid_3 <- "rs9726326"

chr <- 1
bp <- 161216523

merged_data_list <- select_data_merge_LD(rsid,rsid_2,rsid_3,chr,bp,female_90,male_90,
                                        csf_combined_CD84_gene,
                                        eQTLgen_gene,
                                        kosoy_microglia_gene,
                                        csf_combined_FCGR3B_gene,
                                        eQTL_catalogue_cedar_monocytes_gene,
                                        fold_w,file_snp,in_file,pl,in_fold,keep_file,1e6)
merged_female_withld <- merged_data_list[["female"]]
merged_male_withld <- merged_data_list[["male"]]
csf_combined_CD84_gene_with_ld <- merged_data_list[["1"]]
eQTLgen_gene_with_ld <- merged_data_list[["2"]]
kosoy_microglia_gene_with_ld <- merged_data_list[["3"]]
csf_combined_FCGR3B_gene_with_ld <- merged_data_list[["4"]]
eQTL_catalogue_cedar_monocytes_gene_with_ld <- merged_data_list[["5"]]

max_ylim_female <- find_max_ylim_AD(merged_female_withld, rsid)
max_ylim_male <- find_max_ylim_AD(merged_male_withld, rsid)
max_ylim <- max(max_ylim_female, max_ylim_male)
# margin_y 1 7 13 
f <- get_ggplot(max_ylim,bp,merged_female_withld,merged_female_withld,rsid,rsid_2,rsid_3,"AD Female GWAS",13,out_dir)
m <- get_ggplot(max_ylim,bp,merged_female_withld,merged_male_withld,rsid,rsid_2,rsid_3,"AD Male GWAS",13,out_dir)
q1 <- get_ggplot(max_ylim,bp,merged_female_withld,csf_combined_CD84_gene_with_ld,rsid,rsid_2,rsid_3,"csf_combined_CD84_gene_with_ld",7,out_dir)
q2 <- get_ggplot(max_ylim,bp,merged_female_withld,eQTLgen_gene_with_ld,rsid,rsid_2,rsid_3,"eQTLgen_gene_with_ld",1,out_dir)
q3 <- get_ggplot(max_ylim,bp,merged_female_withld,kosoy_microglia_gene_with_ld,rsid,rsid_2,rsid_3,"kosoy_microglia_gene_with_ld",7,out_dir)
q4 <- get_ggplot(max_ylim,bp,merged_female_withld,csf_combined_FCGR3B_gene_with_ld,rsid,rsid_2,rsid_3,"csf_combined_FCGR3B_gene_with_ld",7,out_dir)
q5 <- get_ggplot(max_ylim,bp,merged_female_withld,eQTL_catalogue_cedar_monocytes_gene_with_ld,rsid,rsid_2,rsid_3,"eQTL_catalogue_cedar_monocytes_gene_with_ld",7,out_dir)
# margin_y 0 2 6
c1 <- get_ggplot_compare(merged_female_withld,csf_combined_CD84_gene_with_ld,"left",rsid_3,"csf_combined_CD84_gene_with_ld",
                        FALSE,6,c(0,5,10,15,20,25),c("0","5","10","15","20","25"),out_dir)
c2 <- get_ggplot_compare(merged_female_withld,eQTLgen_gene_with_ld,"middle",rsid_2,"eQTLgen_gene_with_ld",
                        FALSE,0.5,c(0,50,100,150),c("0","50","100","150"),out_dir)
c3 <- get_ggplot_compare(merged_female_withld,kosoy_microglia_gene_with_ld,"right",rsid,"kosoy_microglia_gene_with_ld",
                        FALSE,6,c(0,2,4,6,8,10),c("0","2","4","6","8","10"),out_dir)
c4 <- get_ggplot_compare(merged_female_withld,csf_combined_FCGR3B_gene_with_ld,"middle",rsid_2,"csf_combined_FCGR3B_gene_with_ld",
                        FALSE,6,c(0,5,10,15,20,25),c("0","5","10","15","20","25"),out_dir)
c5 <- get_ggplot_compare(merged_female_withld,eQTL_catalogue_cedar_monocytes_gene_with_ld,"right",rsid,"eQTL_catalogue_cedar_monocytes_gene_with_ld",
                        FALSE,6,c(0,2,4,6,8,10,12),c("0","2","4","6","8","10","12"),out_dir)
maxrows_num <- 8
g <- get_ggplot_genetrack(maxrows_num,bp,merged_female_withld,merged_female_withld,rsid,rsid_2,rsid_3,out_dir,"CD84","NDUFS2","FCER1G","FCGR3B")
















































