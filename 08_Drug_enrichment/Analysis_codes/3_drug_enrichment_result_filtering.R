## Drug Enrichment Filtering Script - August 28, 2025

#load library
library(data.table)
library(tidyverse)
library(patchwork)
library(stringr)
library(gprofiler2)
library(igraph)
library(ggtext)
library(readxl)
library(optparse)

option_list <- list(
  make_option("--work_dir", type = "character", help = "Working directory"),
  make_option("--FDA_gmt", type = "character", help = "FDA GMT file"),
  make_option("--male_list", type = "character", help = "Male drug enrichment XLSX"),
  make_option("--female_list", type = "character", help = "Female drug enrichment XLSX"),
  make_option("--female_out", type = "character", help = "Female filtered output CSV"),
  make_option("--male_out", type = "character", help = "Male filtered output CSV"),
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$work_dir)) stop("--work_dir is required")
dir.create(opt$work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$work_dir)

# Function to convert GeneRatio (e.g., "6/98") to numeric
convert_gene_ratio <- function(gene_ratio) {
  sapply(strsplit(gene_ratio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}


# Read in male and female outputs from Drug Enrichment Analysis - List 1 (No filters) expanded
f1 = read_excel(file.path(opt$work_dir, opt$female_list)) %>%
       mutate(Sex = 'Female', GeneRatio_num = convert_gene_ratio(GeneRatio))

m1 = read_excel(file.path(opt$work_dir, opt$male_list)) %>%
       mutate(Sex = 'Male',   GeneRatio_num = convert_gene_ratio(GeneRatio))

###############################################################################################################################################################
# Read in FDA approaved GMT file
gmt_df_FDA <- read.gmt(file.path(opt$work_dir, "input_files", opt$FDA_gmt))
gmt_df_FDA = as.data.frame(gmt_df_FDA) %>% 
  distinct(term)
## Filter lists to only include FDA approved drugs
# Female
f1_FDA = f1 %>% 
  filter(ID %in% gmt_df_FDA$term) 

# Male
m1_FDA = m1 %>% 
  filter(ID %in% gmt_df_FDA$term)

###############################################################################################################################################################
## Filter 2: FDA approved & p.adjust < 0.05 & GeneRatio Sex1 > 1.5 * GeneRatio Sex2

# Female
f1_FDA_filter2 <- f1_FDA %>%
  filter(p.adjust < 0.05) %>%
  left_join(m1_FDA %>% select(ID, GeneRatio_m = GeneRatio_num), by = "ID") %>%
  filter(is.na(GeneRatio_m) | GeneRatio_num >= 1.5 * GeneRatio_m) %>%
  select(-GeneRatio_m)

# Male
m1_FDA_filter2 <- m1_FDA %>%
  filter(p.adjust < 0.05) %>%
  left_join(f1_FDA %>% select(ID, GeneRatio_f = GeneRatio_num), by = "ID") %>%
  filter(is.na(GeneRatio_f) | GeneRatio_num >= 1.5 * GeneRatio_f) %>%
  select(-GeneRatio_f)


fwrite(f1_FDA_filter2, file.path(opt$work_dir, opt$female_out))
fwrite(m1_FDA_filter2, file.path(opt$work_dir, opt$male_out))

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
#  [1] readxl_1.4.3           ggtext_0.1.2           igraph_2.1.1
#  [4] gprofiler2_0.2.3       patchwork_1.3.0        lubridate_1.9.3
#  [7] forcats_1.0.0          purrr_1.0.2            readr_2.1.5
# [10] tibble_3.2.1           ggplot2_3.5.2          tidyverse_2.0.0
# [13] writexl_1.5.4          epigraphdb_0.2.3       HGNChelper_0.8.15
# [16] stringr_1.5.1          data.table_1.16.4      tidyr_1.3.1
# [19] dplyr_1.1.4            clusterProfiler_4.12.6
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      jsonlite_1.8.9          magrittr_2.0.3
#   [4] farver_2.1.2            fs_1.6.5                zlibbioc_1.50.0
#   [7] vctrs_0.6.5             memoise_2.0.1           ggtree_3.12.0
#  [10] htmltools_0.5.8.1       curl_6.0.1              cellranger_1.1.0
#  [13] gridGraphics_0.5-1      htmlwidgets_1.6.4       plyr_1.8.9
#  [16] httr2_1.0.7             plotly_4.10.4           cachem_1.1.0
#  [19] lifecycle_1.0.4         pkgconfig_2.0.3         Matrix_1.7-1
#  [22] R6_2.5.1                fastmap_1.2.0           gson_0.1.0
#  [25] GenomeInfoDbData_1.2.12 digest_0.6.37           aplot_0.2.3
#  [28] enrichplot_1.24.4       colorspace_2.1-1        AnnotationDbi_1.66.0
#  [31] S4Vectors_0.42.1        RSQLite_2.3.9           fansi_1.0.6
#  [34] timechange_0.3.0        httr_1.4.7              polyclip_1.10-7
#  [37] compiler_4.4.2          bit64_4.5.2             withr_3.0.2
#  [40] BiocParallel_1.38.0     viridis_0.6.5           DBI_1.2.3
#  [43] ggforce_0.4.2           R.utils_2.12.3          MASS_7.3-61
#  [46] rappdirs_0.3.3          tools_4.4.2             splitstackshape_1.4.8
#  [49] ape_5.8                 scatterpie_0.2.4        R.oo_1.27.0
#  [52] glue_1.8.0              nlme_3.1-166            GOSemSim_2.30.2
#  [55] gridtext_0.1.5          grid_4.4.2              shadowtext_0.1.4
#  [58] reshape2_1.4.4          fgsea_1.30.0            generics_0.1.3
#  [61] gtable_0.3.6            tzdb_0.4.0              R.methodsS3_1.8.2
#  [64] hms_1.1.3               tidygraph_1.3.1         xml2_1.3.6
#  [67] utf8_1.2.4              XVector_0.44.0          BiocGenerics_0.50.0
#  [70] ggrepel_0.9.6           pillar_1.9.0            yulab.utils_0.1.8
#  [73] splines_4.4.2           tweenr_2.0.3            treeio_1.28.0
#  [76] lattice_0.22-6          bit_4.5.0.1             tidyselect_1.2.1
#  [79] GO.db_3.19.1            Biostrings_2.72.1       gridExtra_2.3
#  [82] IRanges_2.38.1          stats4_4.4.2            graphlayouts_1.2.1
#  [85] Biobase_2.64.0          stringi_1.8.4           UCSC.utils_1.0.0
#  [88] lazyeval_0.2.2          ggfun_0.1.8             codetools_0.2-20
#  [91] ggraph_2.2.1            qvalue_2.36.0           ggplotify_0.1.2
#  [94] cli_3.6.3               munsell_0.5.1           Rcpp_1.0.13-1
#  [97] GenomeInfoDb_1.40.1     png_0.1-8               parallel_4.4.2
# [100] blob_1.2.4              DOSE_3.30.5             viridisLite_0.4.2
# [103] tidytree_0.4.6          scales_1.3.0            crayon_1.5.3
# [106] rlang_1.1.4             cowplot_1.1.3           fastmatch_1.1-4
# [109] KEGGREST_1.44.1



