## Drug Enrichment Hormone-specific network - August 28, 2025

library(data.table)
library(tidyverse)

# Set WD
setwd(opt$work_dir)
library(optparse)
option_list <- list(
  make_option("--work_dir", type = "character", help = "Working directory"),
  make_option("--hormone_csv", type = "character", help = "Classified sex hormone drugs CSV"),
  make_option("--female_filter2", type = "character", help = "Female filter2 enrichment CSV"),
  make_option("--ppi_df_female", type = "character", help = "Female PPI DF CSV"),
  make_option("--female_out", type = "character", help = "Output female hormone drugs CSV"),
  make_option("--network_out", type = "character", help = "Output network edges CSV for Cytoscape"),
  make_option("--node_out", type = "character", help = "Output node types CSV for Cytoscape")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$work_dir)) stop("--work_dir is required")
dir.create(opt$work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$work_dir)


###############################################################################################################################################################
## Investigate whether sex hormone drugs are in male or female filtered dataset

# Read in Hormone drugs
hormone <- fread(file.path(opt$work_dir, opt$hormone_csv))

# Read in female drug enrichment results (filter2: significant and 1.5-fold)
f1_FDA_filter2 = fread(file.path(opt$work_dir, "results", opt$female_filter2))

# Female
drug_enrich_res <- f1_FDA_filter2 %>%
  filter(ID %in% hormone$ID & Count >= 3)

# Write out drug_enrich_results
fwrite(drug_enrich_res, file.path(opt$work_dir, "results", opt$female_out))


###############################################################################################################################################################

# Get Unique gene list from drug enrichment results
drug_enrich_res_expanded <- drug_enrich_res %>%
  separate_rows(geneID, sep = "/") %>%
  select(ID,geneID)

# Read in PPI related genes from step 1
related_gene <- fread(file.path(opt$work_dir, "network", opt$ppi_df_female))
related_gene <- select(related_gene,g1.name,g2.name)
filtered_related_gene <- related_gene %>%
  filter(g2.name %in% drug_enrich_res_expanded$geneID)

# network
network <- data.frame(Node1=c(drug_enrich_res_expanded$ID,filtered_related_gene$g1.name),
                      Node2=c(drug_enrich_res_expanded$geneID,filtered_related_gene$g2.name))

# Write out network file
fwrite(network, file.path(opt$work_dir, "network", opt$network_out), row.names = FALSE, quote = FALSE)


###############################################################################################################################################################

# Extract unique genes from the two columns
g1_genes <- filtered_related_gene$g1.name %>% unique()
g2_genes <- filtered_related_gene$g2.name %>% unique()

# Find all genes that appeared
all_genes <- union(g1_genes, g2_genes)

# Annotate categories
gene_long <- data.frame(gene = all_genes) %>%
  mutate(source = case_when(
    gene %in% g1_genes & gene %in% g2_genes ~ "both",
    gene %in% g1_genes ~ "pwas_gene",
    gene %in% g2_genes ~ "related_gene"
  ))


# Node type file
ppi_node_type <- data.frame(Node=c(drug_enrich_res$ID,gene_long$gene),
                            Type=c(rep("drug",nrow(drug_enrich_res)),gene_long$source))

# Write out node type file
fwrite(ppi_node_type, file.path(opt$work_dir, "network", opt$node_out), row.names = FALSE, quote = FALSE)


# End of script

sessionInfo()
# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Debian GNU/Linux 12 (bookworm)
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.21.so;  LAPACK version 3.11.0
# 
# locale:
#  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8
#  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8
#  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C
# [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C
# 
# time zone: Etc/UTC
# tzcode source: system (glibc)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
#  [1] lubridate_1.9.3   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4
#  [5] purrr_1.0.2       readr_2.1.5       tidyr_1.3.1       tibble_3.2.1
#  [9] ggplot2_3.5.2     tidyverse_2.0.0   data.table_1.16.4
# 
# loaded via a namespace (and not attached):
#  [1] vctrs_0.6.5      cli_3.6.3        rlang_1.1.4      stringi_1.8.4
#  [5] generics_0.1.3   glue_1.8.0       colorspace_2.1-1 hms_1.1.3
#  [9] scales_1.3.0     fansi_1.0.6      grid_4.4.2       munsell_0.5.1
# [13] tzdb_0.4.0       lifecycle_1.0.4  compiler_4.4.2   timechange_0.3.0
# [17] pkgconfig_2.0.3  R6_2.5.1         tidyselect_1.2.1 utf8_1.2.4
# [21] pillar_1.9.0     magrittr_2.0.3   tools_4.4.2      withr_3.0.2
# [25] gtable_0.3.6