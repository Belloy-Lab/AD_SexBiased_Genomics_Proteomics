#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(clusterProfiler)
  library(rrvgo)
  library(org.Hs.eg.db)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(ComplexHeatmap)
  library(simplifyEnrichment)
})

# Define command-line options
option_list <- list(
  make_option("--work_dir", type = "character", help = "Working directory"),
  make_option("--female_filtered", type = "character", help = "Filtered female GO file"),
  make_option("--male_filtered", type = "character", help = "Filtered male GO file"),
  make_option("--female_out_csv", type = "character", help = "Output CSV for female clusters"),
  make_option("--male_out_csv", type = "character", help = "Output CSV for male clusters"),
  make_option("--female_dendro_pdf", type = "character", help = "Output PDF dendrogram for female"),
  make_option("--male_dendro_pdf", type = "character", help = "Output PDF dendrogram for male"),
  make_option("--female_heatmap_pdf", type = "character", help = "Output PDF heatmap for female"),
  make_option("--male_heatmap_pdf", type = "character", help = "Output PDF heatmap for male")
)
opt <- parse_args(OptionParser(option_list = option_list))

## Females
female_go = fread(file.path(opt$work_dir, opt$female_filtered)) %>% 
  mutate(log10_p = -log10(p.adjust))
go_terms_f <- female_go$ID
sim_matrix_f <- GO_similarity(go_terms_f, ont = "BP", measure = "Rel")
dist_matrix_f <- as.dist(1 - sim_matrix_f)
hclust_f <- hclust(dist_matrix_f, method = "average")
clusters_f <- cutree(hclust_f, k = 7)
cluster_df_f <- data.frame(GO_ID = go_terms_f, Cluster = clusters_f)
female_go_clustered1 <- merge(female_go, cluster_df_f, by.x = "ID", by.y = "GO_ID")

representative_terms_f1a <- female_go_clustered1 %>%
  group_by(Cluster) %>%
  slice(which.min(p.adjust)) %>%
  select(Cluster, Description) %>% 
  rename(`Cluster Name` = Description)

representative_terms_f1a = merge(female_go_clustered1, representative_terms_f1a, by = "Cluster") %>% 
  rename(`Female Gene Ratio` = GeneRatio_num,
         `Female p-value` = pvalue)

fwrite(representative_terms_f1a, file.path(opt$work_dir, opt$female_out_csv))

# Plot dendrogram
pdf(file.path(opt$work_dir, opt$female_dendro_pdf), width = 10, height = 10)
plot(hclust_f, labels = FALSE, hang = -1, main = "Dendrogram of GO Terms")
rect.hclust(hclust_f, k = 7, border = 2:8)
dev.off()

# Heatmap
go_data <- female_go %>%
  select(ID, Description) %>%
  mutate(Cluster = clusters_f)
rownames(sim_matrix_f) <- go_data$Description[match(rownames(sim_matrix_f), go_data$ID)]
colnames(sim_matrix_f) <- go_data$Description[match(colnames(sim_matrix_f), go_data$ID)]

pdf(file.path(opt$work_dir, opt$female_heatmap_pdf), width = 10, height = 10)
Heatmap(
  matrix = sim_matrix_f,
  name = "Similarity",
  col = circlize::colorRamp2(c(0, 1), c("white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_split = clusters_f,
  column_split = clusters_f,
  column_title = "94 Female GO Terms Clustered",
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 4),
  heatmap_legend_param = list(title = "Similarity")
)
dev.off()

rm(list = setdiff(ls(), c("opt", "option_list")))

## Males
male_go = fread(file.path(opt$work_dir, opt$male_filtered)) %>% 
  mutate(log10_p = -log10(p.adjust))
go_terms_m <- male_go$ID
sim_matrix_m <- GO_similarity(go_terms_m, ont = "BP", measure = "Rel")
dist_matrix_m <- as.dist(1 - sim_matrix_m)
hclust_m <- hclust(dist_matrix_m, method = "average")
clusters_m <- cutree(hclust_m, k = 5)
cluster_df_m <- data.frame(GO_ID = go_terms_m, Cluster = clusters_m)
male_go_clustered1 <- merge(male_go, cluster_df_m, by.x = "ID", by.y = "GO_ID")

representative_terms_m1a <- male_go_clustered1 %>%
  group_by(Cluster) %>%
  slice(which.min(p.adjust)) %>%
  select(Cluster, Description) %>% 
  rename(`Cluster Name` = Description)

representative_terms_m1a = merge(male_go_clustered1, representative_terms_m1a, by = "Cluster") %>% 
  rename(`Male Gene Ratio` = GeneRatio_num,
         `Male p-value` = pvalue)

fwrite(representative_terms_m1a, file.path(opt$work_dir, opt$male_out_csv))

pdf(file.path(opt$work_dir, opt$male_dendro_pdf), width = 10, height = 10)
plot(hclust_m, labels = FALSE, hang = -1, main = "Dendrogram of GO Terms")
rect.hclust(hclust_m, k = 5, border = 2:8)
dev.off()

go_data <- male_go %>%
  select(ID, Description) %>%
  mutate(Cluster = clusters_m)
rownames(sim_matrix_m) <- go_data$Description[match(rownames(sim_matrix_m), go_data$ID)]
colnames(sim_matrix_m) <- go_data$Description[match(colnames(sim_matrix_m), go_data$ID)]

pdf(file.path(opt$work_dir, opt$male_heatmap_pdf), width = 8, height = 8)
Heatmap(
  matrix = sim_matrix_m,
  name = "Similarity",
  col = circlize::colorRamp2(c(0, 1), c("white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_split = clusters_m,
  column_split = clusters_m,
  column_title = "10 Male GO Terms Clustered",
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(title = "Similarity")
)
dev.off()