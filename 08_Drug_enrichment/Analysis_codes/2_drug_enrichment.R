## Raw Gene List Drug Enrichment Analysis - August 28, 2025

# Load Libraries
rm(list=ls())
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(stringr)
library(HGNChelper)
library(writexl) # this package is missing in fusion docker

# Set WD
setwd(opt$work_dir)


# import gmt file downloaded from DSigDB
option_list <- list(
  make_option("--work_dir", type = "character", help = "Working directory"),
  make_option("--dsigdb_gmt", type = "character", help = "DSigDB GMT file"),
  make_option("--male_list", type = "character", help = "Male expanded gene list file"),
  make_option("--female_list", type = "character", help = "Female expanded gene list file"),
  make_option("--female_out", type = "character", help = "Female drug enrichment XLSX"),
  make_option("--male_out", type = "character", help = "Male drug enrichment XLSX")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$work_dir)) stop("--work_dir is required")
dir.create(opt$work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$work_dir)

gmt_df_DSigDB <- read.gmt(file.path(opt$work_dir, opt$dsigdb_gmt))
gmt_df_DSigDB$term <- sub("_.*", "", gmt_df_DSigDB$term) # Merge identical drugs from different databases by removing suffixes

# update gene symbol in gmt_df_DSigDB
updated_gene_symbol <- checkGeneSymbols(
  unique(gmt_df_DSigDB$gene),
  chromosome = NULL,
  unmapped.as.na = TRUE,
  map = NULL,
  species = "human",
  expand.ambiguous = FALSE
)
need_updated_gene_symbol <- filter(updated_gene_symbol,updated_gene_symbol$Approved==F)
need_updated_gene_symbol <- filter(need_updated_gene_symbol, !is.na(need_updated_gene_symbol$Suggested.Symbol))
colnames(need_updated_gene_symbol)[1] <- "gene"

gmt_df_DSigDB <- left_join(gmt_df_DSigDB,need_updated_gene_symbol,by='gene')
gmt_df_DSigDB$gene[!is.na(gmt_df_DSigDB$Approved) & gmt_df_DSigDB$Approved == FALSE] <- 
  gmt_df_DSigDB$Suggested.Symbol[!is.na(gmt_df_DSigDB$Approved) & gmt_df_DSigDB$Approved == FALSE]
gmt_df_DSigDB <- gmt_df_DSigDB[,1:2]


# read gene list male/female
gene_male <- read.table(file.path(opt$work_dir, "gene_lists", opt$male_list), header = TRUE)
gene_female <- read.table(file.path(opt$work_dir, "gene_lists", opt$female_list), header = TRUE)

# Check gene symbols in female and male gene lists
check_gene_male <- checkGeneSymbols(
  unique(gene_male$gene_name),
  chromosome = NULL,
  unmapped.as.na = TRUE,
  map = NULL,
  species = "human",
  expand.ambiguous = FALSE
)

check_gene_female <- checkGeneSymbols(
  unique(gene_female$gene_name),
  chromosome = NULL,
  unmapped.as.na = TRUE,
  map = NULL,
  species = "human",
  expand.ambiguous = FALSE
)

# Filter out genes with no valid symbol
check_gene_male = check_gene_male %>% 
  drop_na(Suggested.Symbol)

check_gene_female = check_gene_female %>% 
  drop_na(Suggested.Symbol)


# enrichment analysis gene female
enrich_res_female <- enricher(
  gene = check_gene_female$Suggested.Symbol,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 1,
  maxGSSize = 500,
  qvalueCutoff = 1,
  gson = NULL,
  TERM2GENE = gmt_df_DSigDB,  # gmt_df_DSigDB
  TERM2NAME = NA
)

# Write out female results
enrich_res_female <- as.data.frame(enrich_res_female)
writexl::write_xlsx(enrich_res_female, file.path(opt$work_dir, "results", opt$female_out))


# enrichment analysis gene male
enrich_res_male <- enricher(
  gene = check_gene_male$Suggested.Symbol,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 1,
  maxGSSize = 500,
  qvalueCutoff = 1,
  gson = NULL,
  TERM2GENE = gmt_df_DSigDB,  # gmt_df_DSigDB
  TERM2NAME = NA
)
enrich_res_male <- as.data.frame(enrich_res_male)
writexl::write_xlsx(enrich_res_male, file.path(opt$work_dir, "results", opt$male_out))


## End of code
## Next reference "3_drug_enrichment_result_filtering.R" to apply filters to sex-specifc drugs identifed from analyses.

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
# [1] writexl_1.5.4          epigraphdb_0.2.3       HGNChelper_0.8.15
# [4] stringr_1.5.1          data.table_1.16.4      tidyr_1.3.1
# [7] dplyr_1.1.4            clusterProfiler_4.12.6
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