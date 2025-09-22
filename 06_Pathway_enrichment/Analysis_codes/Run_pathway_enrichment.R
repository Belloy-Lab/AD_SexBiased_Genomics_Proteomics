#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(stringr)
  library(gprofiler2)
  library(igraph)
  library(ggtext)
  library(clusterProfiler)
  library(enrichplot)
  library(ReactomePA)
  library(org.Hs.eg.db)
})

# Define command-line options
option_list <- list(
  make_option("--work_dir", type = "character", help = "Working directory"),
  make_option("--gene_list", type = "character", help = "Path to gene list csv file"),
  make_option("--female_out", type = "character", help = "Female GO output file"),
  make_option("--male_out", type = "character", help = "Male GO output file"),
  make_option("--female_filtered", type = "character", help = "Filtered female GO output"),
  make_option("--male_filtered", type = "character", help = "Filtered male GO output")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Function to convert GeneRatio (e.g., "6/98") to numeric
convert_gene_ratio <- function(gene_ratio) {
  sapply(strsplit(gene_ratio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

## Read in Data 1 without APOE
dat1 = fread(file.path(opt$work_dir, opt$gene_list))

## Female
f1 = dat1 %>% filter(ad_sample == "Female") %>% dplyr::select(gene_name)
f1_genes_a = unique(as.character(f1$gene_name))
f1_genes_a = bitr(f1_genes_a, fromType = "SYMBOL", toType = "ENTREZID", OrgDb ="org.Hs.eg.db")
f1_genes_entrez_a = f1_genes_a$ENTREZID

go_fa = enrichGO(
  gene = f1_genes_entrez_a,
  OrgDb = org.Hs.eg.db,
  keyType = 'ENTREZID',
  ont = 'BP',
  pvalueCutoff = 1,
  pAdjustMethod = 'BH',
  qvalueCutoff = 1,
  readable = TRUE
)

go_fa = as.data.frame(go_fa) %>% arrange(p.adjust)
fwrite(go_fa, file.path(opt$work_dir, opt$female_out), col.names = TRUE)

## Male
m1 = dat1 %>% filter(ad_sample == "Male") %>% dplyr::select(gene_name)
m1_genes_a = unique(as.character(m1$gene_name))
m1_genes_a = bitr(m1_genes_a, fromType = "SYMBOL", toType = "ENTREZID", OrgDb ="org.Hs.eg.db")
m1_genes_entrez_a = m1_genes_a$ENTREZID

go_ma = enrichGO(
  gene = m1_genes_entrez_a,
  OrgDb = org.Hs.eg.db,
  keyType = 'ENTREZID',
  ont = 'BP',
  pvalueCutoff = 1,
  pAdjustMethod = 'BH',
  qvalueCutoff = 1,
  readable = TRUE
)

go_ma = as.data.frame(go_ma) %>% arrange(p.adjust)
fwrite(go_ma, file.path(opt$work_dir, opt$male_out), col.names = TRUE)

rm(list = setdiff(ls(), c("convert_gene_ratio", "opt")))

# Read filtered GO results
f1 = fread(file.path(opt$work_dir, opt$female_out)) %>% mutate(Sex = 'Female')
m1 = fread(file.path(opt$work_dir, opt$male_out)) %>% mutate(Sex = 'Male')

f1_sig <- f1 %>% filter(p.adjust < 0.05) %>% mutate(GeneRatio_num = convert_gene_ratio(GeneRatio))
m1_sig <- m1 %>% filter(p.adjust < 0.05) %>% mutate(GeneRatio_num = convert_gene_ratio(GeneRatio))

all_pathways <- unique(c(f1_sig$ID, m1_sig$ID))

desc_ref <- bind_rows(
  f1 %>% dplyr::select(ID, Description),
  m1 %>% dplyr::select(ID, Description)
) %>% distinct(ID, .keep_all = TRUE)

f1_all <- data.frame(ID = all_pathways) %>%
  left_join(desc_ref, by = "ID") %>%
  left_join(f1 %>% dplyr::select(ID, GeneRatio, p.adjust, pvalue, geneID, Count), by = "ID") %>%
  mutate(Sex = "Female",
         GeneRatio_num = convert_gene_ratio(GeneRatio),
         GeneRatio_num = ifelse(is.na(GeneRatio_num), 0, GeneRatio_num),
         GeneRatio = ifelse(is.na(GeneRatio), 0, GeneRatio),
         p.adjust = ifelse(is.na(p.adjust), 1, p.adjust),
         pvalue = ifelse(is.na(pvalue), 1, pvalue))

m1_all <- data.frame(ID = all_pathways) %>%
  left_join(desc_ref, by = "ID") %>%
  left_join(m1 %>% dplyr::select(ID, GeneRatio, p.adjust, pvalue, geneID, Count), by = "ID") %>%
  mutate(Sex = "Male",
         GeneRatio_num = convert_gene_ratio(GeneRatio),
         GeneRatio_num = ifelse(is.na(GeneRatio_num), 0, GeneRatio_num),
         GeneRatio = ifelse(is.na(GeneRatio), 0, GeneRatio),
         p.adjust = ifelse(is.na(p.adjust), 1, p.adjust),
         pvalue = ifelse(is.na(pvalue), 1, pvalue))

combined_all <- f1_all %>%
  dplyr::select(ID, Description, GeneRatio_num_f = GeneRatio_num, p.adjust_f = p.adjust, geneID_f = geneID, Count) %>%
  left_join(m1_all %>% dplyr::select(ID, GeneRatio_num_m = GeneRatio_num, p.adjust_m = p.adjust, geneID_m = geneID, Count), by = "ID")

filtered_pathways <- combined_all %>%
  filter(
    (p.adjust_f < 0.05 & p.adjust_m >= 0.05 & GeneRatio_num_f >= 1.5 * GeneRatio_num_m) |
      (p.adjust_m < 0.05 & p.adjust_f >= 0.05 & GeneRatio_num_m >= 1.5 * GeneRatio_num_f)
  ) %>%
  pull(ID)

f1_filtered <- f1_all %>% filter(ID %in% filtered_pathways)
m1_filtered <- m1_all %>% filter(ID %in% filtered_pathways)

combined_data <- bind_rows(
  f1_filtered %>% mutate(GeneRatio_plot = GeneRatio_num),
  m1_filtered %>% mutate(GeneRatio_plot = -GeneRatio_num)
)

f_filtered = combined_data %>% 
  filter(Sex == "Female" & p.adjust < 0.05) %>% 
  dplyr::select(-GeneRatio_plot)
fwrite(f_filtered, file.path(opt$work_dir, opt$female_filtered), col.names = TRUE)

m_filtered = combined_data %>% 
  filter(Sex == "Male" & p.adjust < 0.05) %>% 
  dplyr::select(-GeneRatio_plot)
fwrite(m_filtered, file.path(opt$work_dir, opt$male_filtered), col.names = TRUE)