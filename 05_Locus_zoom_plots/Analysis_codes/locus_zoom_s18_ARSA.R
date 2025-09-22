library(data.table)
library(plyr)
library(dplyr)
library(locuszoomr)
library(AnnotationHub)
library(EnsDb.Hsapiens.v86) # package missing
library(ensembldb)
library(GenomicRanges)
library(stats4)
library(S4Vectors)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(Biobase)
library(LDlinkR)
library(openxlsx)
library(readxl)
library(patchwork)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(rtracklayer)
library(locuscomparer)

# wingo_female_pQTL_gene_1mb[
#   CHR == 22
# ][
#   which.min(P)
# ]

female_90 <- fread("/storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update/ADGC_ADSP_UKB_FinnGen_Females_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var")

male_90 <- fread("/storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update/ADGC_ADSP_UKB_FinnGen_Males_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var")

gtf <- rtracklayer::import("/storage1/fs1/belloy/Active/02_Data/01_Incoming/gencode/gencode.v47.basic.annotation.gtf.gz")
ref_seq <- as.data.table(gtf)

#please change the out_dir to your own path to save the output
out_dir <- paste0("/storage2/fs1/belloy2/Active/05_Projects/sivas/LocusZoom/figures/")

fold_w <- paste0("/storage2/fs1/belloy2/Active/05_Projects/sivas/LocusZoom/tmp/")
file_snp <- "_common_id.txt"
pl <- "plink1.9"
in_fold <- "/storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_reference_genotype_data/EU_TOPMed/"
in_file <- "TOPMed_ch"

# please change the out_dir to your own path to source the function file
source("/storage2/fs1/belloy2/Active/04_Code/sivas/LocusZoom/function_locus_zoom_s18.r")

snp_info <- data.frame(rsID = "rs2071421",CHR = 22,BP = 50625988, A1 = "T",A2 = "C",pvalue = 1,
    stringsAsFactors = FALSE
)

wingo_female_pQTL <- fread("/storage1/fs1/belloy/Active/02_Data/02_Processed/QTL/wingo/pQTL/female_DLPFC/female_DLPFC_pQTL.txt")
esngid <- sub("\\..*", "", ref_seq[ref_seq$gene_name == "ARSA", "gene_id"][1])
wingo_female_pQTL_gene <- wingo_female_pQTL  %>% 
    dplyr::filter(gene_id == esngid)

eQTLgen_blood_eQTL <- fread("/storage1/fs1/belloy/Active/02_Data/02_Processed/QTL/eQTLgen/eQTL/blood/blood_eQTL.csv")
esngid <- "ENSG00000100299"
eQTLgen_blood_eQTL_gene <- eQTLgen_blood_eQTL  %>% 
    dplyr::filter(gene_id == esngid)

xQTL_serve_DLPFC_mQTL <- fread("/storage1/fs1/belloy/Active/02_Data/01_Incoming/QTL/xQTLServe/mQTLs/xQTLServe_mQTL_appended_lifted_to38.csv")
esngid <- "cg01556979"
xQTL_serve_DLPFC_mQTL_gene <- xQTL_serve_DLPFC_mQTL  %>% 
    dplyr::filter(CpG == esngid)  %>% 
    dplyr::rename(A1=ALLELE1,A2=ALLELE0,pvalue=p)
xQTL_serve_DLPFC_mQTL_gene <- rbind(xQTL_serve_DLPFC_mQTL_gene,snp_info,fill = TRUE)

chr <- 22
bp_1 <- 50285139
rsid_1 <- "rs34407639"
female_90_1mb <- process_AD_data(female_90,chr,bp_1,1e6)
male_90_1mb <- process_AD_data(male_90,chr,bp_1,1e6)

wingo_female_pQTL_gene_1mb <- process_qtl_data(wingo_female_pQTL_gene,chr,bp_1,1e6)
merged_wingo_female_pQTL_gene_1mb <- merge_qtl(female_90_1mb,male_90_1mb,wingo_female_pQTL_gene_1mb,rsid_1)

eQTLgen_blood_eQTL_gene_1mb <- process_qtl_data(eQTLgen_blood_eQTL_gene,chr,bp_1,1e6)
merged_eQTLgen_blood_eQTL_gene_1mb <- merge_qtl(female_90_1mb,male_90_1mb,eQTLgen_blood_eQTL_gene_1mb,rsid_1)

xQTL_serve_DLPFC_mQTL_gene_1mb <- process_qtl_data(xQTL_serve_DLPFC_mQTL_gene,chr,bp_1,1e6)
merged_xQTL_serve_DLPFC_mQTL_gene_1mb <- merge_qtl(female_90_1mb,male_90_1mb,xQTL_serve_DLPFC_mQTL_gene_1mb,rsid_1)

LD_matrix <- calculate_ld(fold_w,file_snp,rsid_1, pl,in_fold,in_file,chr)

# test_snp = "rs28379706"
# test_snp %in% merged_wingo_female_pQTL_gene_1mb$SNP
# test_snp %in% merged_eQTLgen_blood_eQTL_gene_1mb$SNP
# test_snp %in% merged_xQTL_serve_DLPFC_mQTL_gene_1mb$SNP
# test_snp %in% merged_NGI_female_CSF_pQTL_gene_1mb$SNP
# test_snp %in% merged_NGI_male_CSF_pQTL_gene_1mb$SNP

female_90_1mb_withld_1 <- get_snp_ld(LD_matrix,rsid_1,female_90_1mb)
male_90_1mb_withld_1 <- get_snp_ld(LD_matrix,rsid_1,male_90_1mb)
merged_wingo_female_pQTL_gene_1mb_withld_1 <- get_snp_ld(LD_matrix,rsid_1,merged_wingo_female_pQTL_gene_1mb)
merged_eQTLgen_blood_eQTL_gene_1mb_withld_1 <- get_snp_ld(LD_matrix,rsid_1,merged_eQTLgen_blood_eQTL_gene_1mb)
merged_xQTL_serve_DLPFC_mQTL_gene_1mb_withld_1 <- get_snp_ld(LD_matrix,rsid_1,merged_xQTL_serve_DLPFC_mQTL_gene_1mb)

rsid_2 <- "rs2071421"
bp_2 <- 50625988
female_90_1mb_withld_2 <- get_snp_ld(LD_matrix,rsid_2,female_90_1mb)
male_90_1mb_withld_2 <- get_snp_ld(LD_matrix,rsid_2,male_90_1mb)
merged_wingo_female_pQTL_gene_1mb_withld_2 <- get_snp_ld(LD_matrix,rsid_2,merged_wingo_female_pQTL_gene_1mb)    
merged_eQTLgen_blood_eQTL_gene_1mb_withld_2 <- get_snp_ld(LD_matrix,rsid_2,merged_eQTLgen_blood_eQTL_gene_1mb)
merged_xQTL_serve_DLPFC_mQTL_gene_1mb_withld_2 <- get_snp_ld(LD_matrix,rsid_2,merged_xQTL_serve_DLPFC_mQTL_gene_1mb)

merged_xQTL_serve_DLPFC_mQTL_gene_1mb_withld_1 <- merged_xQTL_serve_DLPFC_mQTL_gene_1mb_withld_1 %>% dplyr::filter(!SNP %in% rsid_2)

max_ylim_female <- ceiling(max(-log10(as.numeric(female_90_1mb$P)), na.rm = TRUE))
max_ylim_male <- ceiling(max(-log10(as.numeric(male_90_1mb$P)), na.rm = TRUE))
max_ylim <- max(max_ylim_female, max_ylim_male)

p1a <- subplot_scatter_top(
    female_90_1mb_withld_1,bp_1,max_ylim,female_90_1mb_withld_1,rsid_1,
    expression(atop(bold("AD GWAS"), bold("female  -") ~ log[10](P))),
    "\n(AD GWAS female top variant)")
p2a <- subplot_scatter_top(
    female_90_1mb_withld_2,bp_2,max_ylim,female_90_1mb_withld_2,rsid_2,
    expression(atop(bold("AD GWAS"), bold("female  -") ~ log[10](P))),
    "\n(Wingo female DLPFC - ARSA pQTL top variant)")

p1b <- subplot_scatter_AD(
    female_90_1mb_withld_1,bp_1,max_ylim,male_90_1mb_withld_1,rsid_1,
    expression(atop(bold("AD GWAS"), bold("male  -") ~ log[10](P))))
p2b <- subplot_scatter_AD(
    female_90_1mb_withld_2,bp_2,max_ylim,male_90_1mb_withld_2, rsid_2,
    expression(atop(bold("AD GWAS"), bold("male  -") ~ log[10](P))))

p1c <- subplot_scatter_qtl(
    female_90_1mb_withld_1,bp_1,merged_wingo_female_pQTL_gene_1mb_withld_1,rsid_1,
    expression(atop(bold("Wingo female DLPFC"), bold("ARSA pQTL  -") ~ log[10](P))),
    "in")
p2c <- subplot_scatter_qtl(
    female_90_1mb_withld_2, bp_2, merged_wingo_female_pQTL_gene_1mb_withld_2, rsid_2,
    expression(atop(bold("Wingo female DLPFC"), bold("ARSA pQTL  -") ~ log[10](P))),
    "in")

p1d <- subplot_scatter_qtl(
    female_90_1mb_withld_1, bp_1, merged_xQTL_serve_DLPFC_mQTL_gene_1mb_withld_1, rsid_1,
    expression(atop(bold("xQTL serve DLPFC"), bold("mQTL  -") ~ log[10](P))),
    "in"
)
p2d <- subplot_scatter_qtl(
    female_90_1mb_withld_2, bp_2, merged_xQTL_serve_DLPFC_mQTL_gene_1mb_withld_2, rsid_2,
    expression(atop(bold("xQTL serve DLPFC"), bold("mQTL  -") ~ log[10](P))),
    "in" #"out"
)

p1e <- subplot_scatter_qtl(
    female_90_1mb_withld_1, bp_1, merged_eQTLgen_blood_eQTL_gene_1mb_withld_1, rsid_1,
    expression(atop(bold("eQTLgen blood"), bold("ARSA eQTL  -") ~ log[10](P))),
    "in"
)
p2e <- subplot_scatter_qtl(
    female_90_1mb_withld_2, bp_2, merged_eQTLgen_blood_eQTL_gene_1mb_withld_2, rsid_2,
    expression(atop(bold("eQTLgen blood"), bold("ARSA eQTL  -") ~ log[10](P))),
    "in"
)

g1 <- subplot_genetrack(female_90_1mb_withld_1,bp_1,rsid_1,0.4,"ARSA")
g2 <- subplot_genetrack(female_90_1mb_withld_2,bp_2,rsid_2,0.4,"ARSA")

sc1 <- subplot_compare(female_90_1mb_withld_2,merged_wingo_female_pQTL_gene_1mb_withld_2,chr,rsid_2,
    "AD GWAS female  -",
    "Wingo female DLPFC - ARSA pQTL  -",
    "topleft")
sc2 <- subplot_compare(female_90_1mb_withld_1,merged_xQTL_serve_DLPFC_mQTL_gene_1mb_withld_1,chr,rsid_1,
    "AD GWAS female  -",
    "xQTL serve DLPFC - mQTL  -",
    "topleft")
sc3 <- subplot_compare(female_90_1mb_withld_2,merged_eQTLgen_blood_eQTL_gene_1mb_withld_2,chr,rsid_2,
    "AD GWAS female  -",
    "eQTLgen blood - ARSA eQTL  -",
    "topleft")
    
#upper part
top_plot <- wrap_plots(
    p1a, p2a,
    p1b, p2b,
    p1c, p2c,
    p1d, p2d,
    p1e, p2e,
    g1,  g2,
  ncol = 2
)

#lower part
bottom_plot <- wrap_plots(
  sc1, sc2, sc3,
  ncol = 3
)

#combine the two parts
final_plot <- top_plot / bottom_plot +
  plot_layout(
    heights = c(
      sum(rep(1,5)) + 0.5,  # upper part (6 rows + 0.5 for spacing)
      2.25                  # lower part (1 row)
    )
  )


index = 2
# Create title with slightly increased top margin
title_S5 <- ggdraw() + 
    draw_label(
        paste0("S18.", index), 
        size = 18,  # Reduce font size
        fontface = "bold",
        x = 0.02, 
        hjust = 0, 
        vjust = 1
    ) +
    theme(plot.margin = margin(0.6, 0, 0.2, 0.5, unit = "inches"))  # Slightly increased top margin

title_global <- ggdraw() + 
    draw_label(
        paste0("ARSA"),
        size = 18,  
        fontface = 'bold', 
        x = 0.5, 
        hjust = 0.5,
        vjust = 0,
        lineheight = 0.8  
    ) +
    theme(plot.margin = margin(0, 0, 0, 0, unit = "inches"))  # More spacing below title
# Adjust overall plot margins to increase top and bottom space
final_plot_with_margin <- ggdraw(final_plot) + 
    theme(plot.margin = margin(0, 0.3, 0.5, 0.3, unit = "inches"))  # Increased top and bottom margins

# Combine section title, global title, and plot into a single layout
final_plot_with_title <- plot_grid(
    title_S5,        # Section title
    title_global,    # Global title
    final_plot_with_margin,  # Main plot area
    ncol = 1,
    rel_heights = c(0.08, 0.02, 0.9)  # Allocate heights: 7% for section, 5% for global title, 88% for plot
)

# Save the final output as a PDF file with updated dimensions
ggsave(
    filename = paste0(out_dir,index,"_","ARSA" , ".pdf"), 
    plot = final_plot_with_title, 
    width = 8.27 * 2, 
    height = 11.69 * 2
)
