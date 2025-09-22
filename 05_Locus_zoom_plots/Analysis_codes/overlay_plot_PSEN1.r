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

female_90 <- fread("/storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update/ADGC_ADSP_UKB_FinnGen_Females_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var")

male_90 <- fread("/storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update/ADGC_ADSP_UKB_FinnGen_Males_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var")

w23_pQTL_female <- fread(paste0("/storage1/fs1/belloy/Active/05_Projects/noahc/opt_QTL/liftover/lifted_data/wingo_female_pQTL_appended_lifted_to38.txt"))

w23_pQTL_male <- fread(paste0("/storage1/fs1/belloy/Active/05_Projects/noahc/opt_QTL/liftover/lifted_data/wingo_male_pQTL_appended_lifted_to38.txt"))

w23_pQTL_combined <- fread(paste0("/storage1/fs1/belloy/Active/05_Projects/noahc/opt_QTL/liftover/lifted_data/wingo_both_pQTL_appended_lifted_to38.txt"))

brain_meta_DLPFC <- fread("/storage1/fs1/belloy/Active/05_Projects/noahc/opt_QTL/liftover/lifted_data/BrainMeta_all_mQTL_appended_lifted_to38.txt")

wingo_female_pQTL <- fread("/storage1/fs1/belloy/Active/02_Data/02_Processed/QTL/wingo/pQTL/female_DLPFC/female_DLPFC_pQTL.txt")

eQTLgen <- fread("/storage1/fs1/belloy/Active/02_Data/02_Processed/QTL/eQTLgen/eQTL/blood/blood_eQTL.csv")

eQTL_catalogue_nedelec_naive <- fread("/storage1/fs1/belloy/Active/02_Data/02_Processed/QTL/eQTL_catalogue/eQTL/nedelec_naive/nedelec_naive_eQTL.tsv.gz")

gtf <- rtracklayer::import("/storage1/fs1/belloy/Active/02_Data/01_Incoming/gencode/gencode.v47.basic.annotation.gtf.gz")
ref_seq <- as.data.table(gtf)

gene_f <- fread(paste0("/storage1/fs1/belloy/Active/05_Projects/dreid/sex_strat_AD_EU_GWAS_T-PWAS_2024/",
                       "sex-strat_AD_EU_ADGC_ADSP_UKB_FinnGen_meta/",
                       "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.PWAS_non-strat_sex-strat_W23nCSFc_female-specific-genes-qtl-list.txt"))

# define your own output directory
out_dir <- paste0("/storage2/fs1/belloy2/Active/05_Projects/sivas/LocusZoom/figures/")

fold_w <- paste0("/storage2/fs1/belloy2/Active/05_Projects/sivas/LocusZoom/tmp/")
file_snp <- "_common_id.txt"
pl <- "plink1.9"
in_fold <- "/storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_reference_genotype_data/EU_TOPMed/"
in_file <- "TOPMed_ch"

in_fold <- "/storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_reference_genotype_data/EU_TOPMed/"
in_file <- "TOPMed_ch"
keep_file <- "/storage1/fs1/belloy/Active/01_References/01_Files/EU_nHISP_TOPMed_list_Unrel.txt"

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

# source the function file
source("/storage2/fs1/belloy2/Active/04_Code/sivas/LocusZoom/functions_overlay_plot.r")

rsid <- "rs3742825"
rsid_2 <- "rs74986264"
rsid_3 <- "rs11159021" 
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

source("/storage1/fs1/belloy/Active/05_Projects/chenyu.y/overlay_plot/functions_overlay_plot.r")
max_ylim_female <- find_max_ylim_AD(merged_female_withld, rsid)
max_ylim_male <- find_max_ylim_AD(merged_male_withld, rsid)
max_ylim <- max(max_ylim_female, max_ylim_male)
# margin_y 1 7 13 
f <- get_ggplot(max_ylim,bp,merged_female_withld,merged_female_withld,rsid,rsid_2,rsid_3,"AD Female GWAS",13)
m <- get_ggplot(max_ylim,bp,merged_female_withld,merged_male_withld,rsid,rsid_2,rsid_3,"AD Male GWAS",13)
q1 <- get_ggplot(max_ylim,bp,merged_female_withld,merged_w23_pQTL_female_gene_with_ld,rsid,rsid_2,rsid_3,"w23_pQTL_female",13)
q2 <- get_ggplot(max_ylim,bp,merged_female_withld,merged_wingo_female_pQTL_gene_with_ld,rsid,rsid_2,rsid_3,"wingo_female_DLPFC_pQTL_Female_PSEN1_PAPLN",13)
q3 <- get_ggplot(max_ylim,bp,merged_female_withld,merged_eQTLgen_gene_with_ld,rsid,rsid_2,rsid_3,"eQTLgen_blood_eQTL_Female_PSEN1_HEATR4",7)
q4 <- get_ggplot(max_ylim,bp,merged_female_withld,merged_brain_meta_DLPFC_gene_with_ld,rsid,rsid_2,rsid_3,"BrainMeta_DLPFC_mQTL_Female_PSEN1_cg24036523",7)
q5 <- get_ggplot(max_ylim,bp,merged_female_withld,merged_eQTL_catalogue_nedelec_naive_gene_with_ld,rsid,rsid_2,rsid_3,"eQTL_catalogue_nedelec_naive_eQTL_Female_PSEN1_PSEN1",13)
# margin_y 0 2 6
c1 <- get_ggplot_compare(merged_female_withld,merged_w23_pQTL_female_gene_with_ld,"middle",rsid,"w23_pQTL_female",
                        FALSE,11.5,c(0,2,4,6,8),c("0","2","4","6","8"))
c2 <- get_ggplot_compare(merged_female_withld,merged_wingo_female_pQTL_gene_with_ld,"left",rsid_2,"wingo_female_DLPFC",
                        FALSE,11.5,c(0,1,2,3,4,5,6),c("0","1","2","3","4","5","6"))   
c3 <- get_ggplot_compare(merged_female_withld,merged_eQTLgen_gene_with_ld,"right",rsid_3,"eQTLgen_blood_eQTL",
                        FALSE,6,c(0,5,10,15,20),c("0","5","10","15","20"))
c4 <- get_ggplot_compare(merged_female_withld,merged_brain_meta_DLPFC_gene_with_ld,"middle",rsid,"BrainMeta_DLPFC_mQTL",
                        FALSE,6,c(2,4,6,8,10),c("2","4","6","8","10"))
c5 <- get_ggplot_compare(merged_female_withld,merged_eQTL_catalogue_nedelec_naive_gene_with_ld,"middle",rsid,"eQTL_catalogue_nedelec_naive_eQTL",
                        FALSE,11.5,c(0,1,2,3,4,5,6),c("0","1","2","3","4","5","6"))   
maxrows_num <- 15
g <- get_ggplot_genetrack(maxrows_num,bp,merged_female_withld,merged_female_withld,rsid,rsid_2,rsid_3,"PSEN1","PAPLN","HEATR4","cg24036523")
























































