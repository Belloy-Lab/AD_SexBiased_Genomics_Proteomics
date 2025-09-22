## Expanded Gene List Drug Enrichment Analysis - August 28, 2025

# Load Libraries
rm(list=ls())
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(HGNChelper) # this package is missing in fusion docker
library(epigraphdb) # this package is missing in fusion docker

# Set WD
setwd("/storage2/fs1/belloy2/Active/05_Projects/sivas/Drug_enrichment")


option_list <- list(
  make_option("--work_dir", type = "character", help = "Working directory"),
  make_option("--male_list", type = "character", help = "Path to male gene list file"),
  make_option("--female_list", type = "character", help = "Path to female gene list file"),
  make_option("--male_out", type = "character", help = "Male expanded output file (tsv)"),
  make_option("--female_out", type = "character", help = "Female expanded output file (tsv)")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$work_dir)) stop("--work_dir is required")
dir.create(opt$work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$work_dir)


##################################################################################################################
## Male
###############################################################################################################################
# Read in Male gene list (filtered 2 gene list = Genes with priority score 1 + MAPT [lit support])
gene_male <- read.table(file.path(opt$work_dir, "gene_lists", opt$male_list), header = TRUE)

# Check Gene Symbols
check_gene_male <- checkGeneSymbols(
  unique(gene_male$gene_name),
  chromosome = NULL,
  unmapped.as.na = TRUE,
  map = NULL,
  species = "human",
  expand.ambiguous = FALSE
)
# All FALSE genes had NA for SUggested Symbol -- No changes needed
colnames(gene_male) <- "hgnc_names"


# druggable tiers file (same as the package)
druggable_tiers <- readxl::read_excel("aag1166_Table S1.xlsx")
druggable_tiers <- as.data.frame(druggable_tiers) %>% 
  dplyr::select(hgnc_names,druggability_tier)

# Check Gene Symbols in Druggability tiers file
check_gene_dt <- checkGeneSymbols(
  unique(druggable_tiers$hgnc_names),
  chromosome = NULL,
  unmapped.as.na = TRUE,
  map = NULL,
  species = "human",
  expand.ambiguous = FALSE
)

check_gene_dt = check_gene_dt %>%
  dplyr::rename(hgnc_names = x)

# Merge updated symbold to druggability terms df
druggable_tiers2 <- inner_join(druggable_tiers,check_gene_dt,by='hgnc_names') %>% 
  drop_na(Suggested.Symbol) %>%
  dplyr::select(Suggested.Symbol, druggability_tier) %>% 
  dplyr::rename(hgnc_names = Suggested.Symbol)


# find druggable info for raw gene list
candidate_gene_tiers <- data.frame(hgnc_names=gene_male$hgnc_names)
candidate_gene_tiers <- left_join(candidate_gene_tiers,druggable_tiers2)

# find related genes (ppi) + druggable info
ppi_df_all <- data.frame()
for (i in 1:length(gene_male$hgnc_names)) {
  ppi_df <- query_epigraphdb(
    route = "/gene/druggability/ppi",
    params = list(gene_name = c(gene_male$hgnc_names[i])),
    mode = "table",
    method = "GET"
  )
  ppi_df_all <- rbind(ppi_df_all,ppi_df)
}

ppi_df_all$disease_related <- NA  # Create a new empty column


# Loop through each interacting gene in g2.name
for (i in seq_len(nrow(ppi_df_all))) {
  gene <- ppi_df_all$g2.name[i]
  
  res <- query_epigraphdb(
    route = "/gene/literature",
    params = list(gene_name = gene, object_name = "Alzheimer's disease"),
    mode = "table"
  )
  
  # If there is literature evidence linking to Alzheimer's, mark as Yes
  ppi_df_all$disease_related[i] <- if (nrow(res) > 0) "Yes" else "No"
}
colnames(ppi_df_all)[1] <- 'hgnc_names'

# merge 
res <- left_join(candidate_gene_tiers,ppi_df_all,by='hgnc_names') 

res_filter = res %>% 
  filter(disease_related == "Yes" & g2.druggability_tier == "Tier 1") %>% 
  dplyr::select(`g2.name`) %>% 
  dplyr::rename(hgnc_names = `g2.name`) 


# Add expanded genes to raw gene list - get unique gene_names
m_full = rbind(gene_male, res_filter) %>% 
  distinct(hgnc_names) %>% 
  dplyr::rename(gene_name = hgnc_names)

# Write out expanded gene_list
fwrite(m_full, file.path(opt$work_dir, "gene_lists", opt$male_out), col.names = TRUE, row.names = FALSE)



rm(list=ls())

##################################################################################################################
## Repeat for Female
##################################################################################################################

# Read in Female gene list (filtered 2 gene list = Genes with priority score 1 + SORL1, ALPL, HLA [lit support])
gene_female <- read.table(file.path(opt$work_dir, "gene_lists", opt$female_list), header = TRUE)

# Check Gene Symbols
check_gene_female <- checkGeneSymbols(
  unique(gene_female$gene_name),
  chromosome = NULL,
  unmapped.as.na = TRUE,
  map = NULL,
  species = "human",
  expand.ambiguous = FALSE
)

colnames(gene_female) <- "hgnc_names"


# druggable tiers file (same as the package)
druggable_tiers <- readxl::read_excel("aag1166_Table S1.xlsx")
druggable_tiers <- as.data.frame(druggable_tiers) %>% 
  dplyr::select(hgnc_names,druggability_tier)

# Check Gene Symbols in Druggability tiers file
check_gene_dt <- checkGeneSymbols(
  unique(druggable_tiers$hgnc_names),
  chromosome = NULL,
  unmapped.as.na = TRUE,
  map = NULL,
  species = "human",
  expand.ambiguous = FALSE
)

check_gene_dt = check_gene_dt %>%
  dplyr::rename(hgnc_names = x)

# Merge updated symbold to druggability terms df
druggable_tiers2 <- inner_join(druggable_tiers,check_gene_dt,by='hgnc_names') %>% 
  drop_na(Suggested.Symbol) %>%
  dplyr::select(Suggested.Symbol, druggability_tier) %>% 
  dplyr::rename(hgnc_names = Suggested.Symbol)


# find druggable info for raw gene list
candidate_gene_tiers <- data.frame(hgnc_names=gene_female$hgnc_names)
candidate_gene_tiers <- left_join(candidate_gene_tiers,druggable_tiers2)


# find related genes (ppi) + druggable info
ppi_df_all <- data.frame()
for (i in 1:length(gene_female$hgnc_names)) {
  ppi_df <- query_epigraphdb(
    route = "/gene/druggability/ppi",
    params = list(gene_name = c(gene_female$hgnc_names[i])),
    mode = "table",
    method = "GET"
  )
  ppi_df_all <- rbind(ppi_df_all,ppi_df)
}

ppi_df_all$disease_related <- NA  # Create a new empty column



# Loop through each interacting gene in g2.name
for (i in seq_len(nrow(ppi_df_all))) {
  gene <- ppi_df_all$g2.name[i]
  
  res <- query_epigraphdb(
    route = "/gene/literature",
    params = list(gene_name = gene, object_name = "Alzheimer's disease"),
    mode = "table"
  )
  
  # If there is literature evidence linking to Alzheimer's, mark as Yes
  ppi_df_all$disease_related[i] <- if (nrow(res) > 0) "Yes" else "No"
}
colnames(ppi_df_all)[1] <- 'hgnc_names'

# merge 
res <- left_join(candidate_gene_tiers,ppi_df_all,by='hgnc_names') 

res_filter = res %>% 
  filter(disease_related == "Yes" & g2.druggability_tier == "Tier 1") %>% 
  dplyr::select(`g2.name`) %>% 
  dplyr::rename(hgnc_names = `g2.name`) 


# Add expanded genes to raw gene list - get unique gene_names
f_full = rbind(gene_female, res_filter) %>% 
  distinct(hgnc_names) %>% 
  dplyr::rename(gene_name = hgnc_names)

# Write out expanded gene_list
fwrite(f_full, file.path(opt$work_dir, "gene_lists", opt$female_out), col.names = TRUE, row.names = FALSE)



rm(list=ls())


## End of code
## Next reference "2_drug_enrichment_script.R" to run enrichment analysis on the expanded gene lists

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
# [1] epigraphdb_0.2.3       HGNChelper_0.8.15      stringr_1.5.1
# [4] data.table_1.16.4      tidyr_1.3.1            dplyr_1.1.4
# [7] clusterProfiler_4.12.6
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.3               gson_0.1.0              shadowtext_0.1.4
#   [4] gridExtra_2.3           httr2_1.0.7             readxl_1.4.3
#   [7] rlang_1.1.4             magrittr_2.0.3          DOSE_3.30.5
#  [10] compiler_4.4.2          RSQLite_2.3.9           png_0.1-8
#  [13] vctrs_0.6.5             reshape2_1.4.4          pkgconfig_2.0.3
#  [16] crayon_1.5.3            fastmap_1.2.0           XVector_0.44.0
#  [19] ggraph_2.2.1            splitstackshape_1.4.8   utf8_1.2.4
#  [22] enrichplot_1.24.4       UCSC.utils_1.0.0        purrr_1.0.2
#  [25] bit_4.5.0.1             zlibbioc_1.50.0         cachem_1.1.0
#  [28] aplot_0.2.3             GenomeInfoDb_1.40.1     jsonlite_1.8.9
#  [31] blob_1.2.4              BiocParallel_1.38.0     tweenr_2.0.3
#  [34] parallel_4.4.2          R6_2.5.1                stringi_1.8.4
#  [37] RColorBrewer_1.1-3      cellranger_1.1.0        GOSemSim_2.30.2
#  [40] Rcpp_1.0.13-1           R.utils_2.12.3          IRanges_2.38.1
#  [43] Matrix_1.7-1            splines_4.4.2           igraph_2.1.1
#  [46] tidyselect_1.2.1        qvalue_2.36.0           viridis_0.6.5
#  [49] codetools_0.2-20        curl_6.0.1              lattice_0.22-6
#  [52] tibble_3.2.1            plyr_1.8.9              Biobase_2.64.0
#  [55] treeio_1.28.0           withr_3.0.2             KEGGREST_1.44.1
#  [58] gridGraphics_0.5-1      scatterpie_0.2.4        polyclip_1.10-7
#  [61] Biostrings_2.72.1       pillar_1.9.0            ggtree_3.12.0
#  [64] stats4_4.4.2            ggfun_0.1.8             generics_0.1.3
#  [67] S4Vectors_0.42.1        ggplot2_3.5.2           munsell_0.5.1
#  [70] scales_1.3.0            tidytree_0.4.6          glue_1.8.0
#  [73] lazyeval_0.2.2          tools_4.4.2             fgsea_1.30.0
#  [76] fs_1.6.5                graphlayouts_1.2.1      fastmatch_1.1-4
#  [79] tidygraph_1.3.1         cowplot_1.1.3           grid_4.4.2
#  [82] ape_5.8                 AnnotationDbi_1.66.0    colorspace_2.1-1
#  [85] nlme_3.1-166            GenomeInfoDbData_1.2.12 patchwork_1.3.0
#  [88] ggforce_0.4.2           cli_3.6.3               rappdirs_0.3.3
#  [91] fansi_1.0.6             viridisLite_0.4.2       gtable_0.3.6
#  [94] R.methodsS3_1.8.2       yulab.utils_0.1.8       digest_0.6.37
#  [97] BiocGenerics_0.50.0     ggrepel_0.9.6           ggplotify_0.1.2
# [100] farver_2.1.2            memoise_2.0.1           R.oo_1.27.0
# [103] lifecycle_1.0.4         httr_1.4.7              GO.db_3.19.1
# [106] bit64_4.5.2             MASS_7.3-61
