## Cell-specific enrichment analysis - figure creation
# August 6, 2025

library(data.table)
library(tidyverse)
library(patchwork)
library(extrafont)
library(optparse)

option_list <- list(
  make_option("--work_dir", type = "character", help = "Working directory"),
  make_option("--results_in", type = "character", help = "Input enrichment results CSV"),
  make_option("--out_fig", type = "character", help = "Output 3-cell barplot (jpg)"),
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$work_dir)) stop("--work_dir is required")
dir.create(opt$work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$work_dir)

## Read in comprehensive results
all = fread(file.path(opt$work_dir, "results", opt$results_in), header = TRUE) %>% 
  dplyr::select(-V8:-V19)

# Filter to include results when using full gene list
all1 = all %>% 
  filter(gene_set == "All Female" | gene_set == "All Male") %>% 
  mutate(Sex = case_when(gene_set == "All Female" ~ "Female",
                         T ~ "Male"))

# Create the horizontal bar plot
p1 <- ggplot(all1, aes(x = fold_enrichment, y = cell_type, fill = Sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Female" = "#700000", "Male" = "#3780b5")) +
  geom_text(aes(label = ifelse(p_value < 0.05, "*", ""), 
                x = fold_enrichment + ifelse(fold_enrichment >= 0, 0.05, -0.05)),
            position = position_dodge(width = 0.9), 
            vjust = 0.75, hjust = 0.3, size = 8, fontface = "bold") +
  theme_minimal() +
  labs(x = "Fold Enrichment", y = "Cell Type", fill = "Sex") +
  theme(
    axis.text = element_text(size = 12, face = "bold", family = "Arial"),
    axis.title = element_text(size = 14, face = "bold", family = "Arial"),
    legend.text = element_text(size = 12, face = "bold", family = "Arial"),
    legend.title = element_text(size = 14, face = "bold", family = "Arial"),
    legend.position = "top",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(5, "pt"),
    plot.margin = margin(t = 15, r = 10, b = 15, l = 10, unit = "pt"))

# Print the plot
# print(p1)

# Filter to include results for 3 main cell-types
all2 = all1 %>% 
  filter(gene_set == "All Female" | gene_set == "All Male") %>% 
  filter(cell_type != "Endothelial Cells" & cell_type != "Oligodendrocytes") %>% 
  mutate(Sex = case_when(gene_set == "All Female" ~ "Female",
                         T ~ "Male"))

# Create the horizontal bar plot
p2 <- ggplot(all2, aes(x = fold_enrichment, y = cell_type, fill = Sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Female" = "#700000", "Male" = "#3780b5")) +
  geom_text(aes(label = ifelse(p_value < 0.05, "*", ""), 
                x = fold_enrichment + ifelse(fold_enrichment >= 0, 0.05, -0.05)),
            position = position_dodge(width = 0.9), 
            vjust = 0.75, size = 8, fontface = "bold") +
  theme_minimal() +
  labs(x = "Fold Enrichment", y = "Cell Type", fill = "Sex") +
  theme(
    axis.text = element_text(size = 12, face = "bold", family = "Arial"),
    axis.title = element_text(size = 14, face = "bold", family = "Arial"),
    legend.text = element_text(size = 12, face = "bold", family = "Arial"),
    legend.title = element_text(size = 14, face = "bold", family = "Arial"),
    legend.position = "top",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(5, "pt"),
    plot.margin = margin(t = 15, r = 10, b = 15, l = 10, unit = "pt"))

# Print the plot
print(p2)

# Save plot
file_name <- opt$out_fig
out_dir <- file.path(opt$work_dir, "results")

file_path <- file.path(out_dir, file_name)
ggsave(filename = file_path, plot = p2, device = "jpg", width = 6, height = 6, dpi = 320)

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
#  [1] extrafont_0.19    patchwork_1.3.0   lubridate_1.9.3   forcats_1.0.0
#  [5] stringr_1.5.1     dplyr_1.1.4       purrr_1.0.2       readr_2.1.5
#  [9] tidyr_1.3.1       tibble_3.2.1      ggplot2_3.5.2     tidyverse_2.0.0
# [13] data.table_1.16.4
# 
# loaded via a namespace (and not attached):
#  [1] gtable_0.3.6      compiler_4.4.2    tidyselect_1.2.1  textshaping_0.4.1
#  [5] systemfonts_1.1.0 scales_1.3.0      R6_2.5.1          labeling_0.4.3
#  [9] generics_0.1.3    munsell_0.5.1     pillar_1.9.0      tzdb_0.4.0
# [13] rlang_1.1.4       utf8_1.2.4        Rttf2pt1_1.3.12   stringi_1.8.4
# [17] timechange_0.3.0  cli_3.6.3         withr_3.0.2       magrittr_2.0.3
# [21] grid_4.4.2        hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5
# [25] glue_1.8.0        extrafontdb_1.0   farver_2.1.2      ragg_1.3.3
# [29] fansi_1.0.6       colorspace_2.1-1  tools_4.4.2       pkgconfig_2.0.3