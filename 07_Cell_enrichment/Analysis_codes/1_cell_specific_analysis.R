# Cell-type specific Enrichment Analysis Script
# August 4, 2025

# Load necessary libraries
library(data.table)
library(tidyverse)
library(HGNChelper)
library(optparse)
option_list <- list(
  make_option("--work_dir", type = "character", help = "Working directory"),
  make_option("--gene_list", type = "character", help = "Prioritized gene list CSV"),
  make_option("--background_list", type = "character", help = "Background gene list CSV"),
  make_option("--results_out", type = "character", help = "Cell-specific enrichment results CSV")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$work_dir)) stop("--work_dir is required")
dir.create(opt$work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$work_dir)


setwd(opt$work_dir)
################################################################################################################################################################
## Read in Male and Female Gene lists - full list - not expanded
gene_list = fread(file.path(opt$work_dir, "gene_lists", opt$gene_list))

# update gene symbol in prioritized gene list
updated_gene_list <- checkGeneSymbols(
  unique(gene_list$gene_name),
  chromosome = NULL,
  unmapped.as.na = TRUE,
  map = NULL,
  species = "human",
  expand.ambiguous = FALSE) 

# Filter to only approved gene symbols and rename columns
updated_gene_list =  as.data.frame(updated_gene_list) %>% 
  drop_na(Suggested.Symbol) %>% 
  dplyr::rename(gene_name = x,
                updated_gene_name = Suggested.Symbol) 

# Merge updated gene names with original gene list
final_gene_list = inner_join(gene_list, updated_gene_list, by = "gene_name") %>% 
  dplyr::select(-gene_name) %>% 
  dplyr::rename(gene_name = updated_gene_name)

# Get Female-specific gene list
f1 = final_gene_list %>% 
 # filter(discovery == "DAL_Brain" | discovery == "DAL_CSF") %>% # Customize filter to only look at sub-sets of genes
  filter(ad_sample == "Female") %>% 
  distinct(gene_name, .keep_all = TRUE)
n_f = n_distinct(f1$gene_name) # N = 112

# Get Male-specific gene list
m1 = final_gene_list %>% 
 # filter(discovery == "DAL_Brain" | discovery == "DAL_CSF") %>% # Customize filter to only look at sub-sets of genes
  filter(ad_sample == "Male") %>% 
  distinct(gene_name, .keep_all = TRUE)
n_m = n_distinct(m1$gene_name) # N = 17


################################################################################################################################################################
## Read in Background List from Dan's supplemental table S28
bg = fread(file.path(opt$work_dir, "gene_lists", opt$background_list)) %>% 
  dplyr::rename(gene_name = Gene,
                cell_type = `Max Cell Type`) %>% 
  group_by(gene_name, cell_type) %>%
  distinct(gene_name, .keep_all = TRUE)

n_distinct(bg$gene_name) # N = 19,117 genes in background list # SAT: I get 20,559 check with Noah
# update gene symbol in background gene list
updated_bg <- checkGeneSymbols(
  unique(bg$gene_name),
  chromosome = NULL,
  unmapped.as.na = TRUE,
  map = NULL,
  species = "human",
  expand.ambiguous = FALSE) 

# Filter to only approved gene symbols and rename columns
updated_bg =  as.data.frame(updated_bg) %>% 
  drop_na(Suggested.Symbol) %>% 
  dplyr::rename(gene_name = x,
                updated_gene_name = Suggested.Symbol) 

# Merge updated gene names with original gene list
final_bg = inner_join(bg, updated_bg, by = "gene_name") %>% 
  dplyr::rename(old_gene_name = gene_name,
                gene_name = updated_gene_name) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene_name, cell_type)


n_bg = n_distinct(final_bg$gene_name) # N = 19,117 genes in background list # SAT: the number id correct here.
unique(final_bg$cell_type) # 5 cell-types in data



################################################################################################################################################################
# Set empty df for results
results = data.frame()

########################################################
## Mature Astrocytes
ast = final_bg %>% 
  filter(cell_type == "Human mature astrocytes") %>% 
  distinct(gene_name)

# Make sure there are no duplicates
n_distinct(ast$gene_name) == nrow(ast)
n_ast = n_distinct(ast$gene_name) # N = 1,869


########################################################
# Merge female genes with mature astrocytes
merge_f_ast = inner_join(ast, f1, by = "gene_name") 
n_f_ast = n_distinct(merge_f_ast$gene_name) # N = 18

# Calculate p-value using hypergeometric test
f_ast_p <- phyper(n_f_ast - 1, n_ast, n_bg - n_ast, n_f, lower.tail = FALSE)

# Calculate the GeneRatio
f_ast_gr <- n_f_ast / n_ast

# Calculate the fold enrichment
f_ast_fe <- ((n_f_ast / n_f) / (n_ast / n_bg)) - 1

# Store results in temp df
cell_results = data.frame(
  gene_set = "All Female",
  cell_type = "Mature Astrocytes",
  overlap = n_f_ast,
  percent_celltype = (n_f_ast/n_f)*100,
  gene_ratio = f_ast_gr,
  fold_enrichment = f_ast_fe,
  p_value = f_ast_p)

# Append results to final df
results = rbind(results, cell_results)


########################################################
# Merge male genes with mature astrocytes
merge_m_ast = inner_join(ast, m1, by = "gene_name") 
n_m_ast = n_distinct(merge_m_ast$gene_name) # N = 2

# Calculate p-value using hypergeometric test
m_ast_p <- phyper(n_m_ast - 1, n_ast, n_bg - n_ast, n_m, lower.tail = FALSE)

# Calculate the GeneRatio
m_ast_gr <- n_m_ast / n_ast

# Calculate the fold enrichment
m_ast_fe <- ((n_m_ast / n_m) / (n_ast / n_bg)) - 1

# Store results in temp df
cell_results = data.frame(
  gene_set = "All Male",
  cell_type = "Mature Astrocytes",
  overlap = n_m_ast,
  percent_celltype = (n_m_ast/n_m)*100,
  gene_ratio = m_ast_gr,
  fold_enrichment = m_ast_fe,
  p_value = m_ast_p)

# Append results to final df
results = rbind(results, cell_results)


########################################################################################################################################################################
## Endothelial Cells
end = final_bg %>% 
  filter(cell_type == "Human Endothelial") %>% 
  distinct(gene_name)

# Make sure there are no duplicates
n_distinct(end$gene_name) == nrow(end)
n_end = n_distinct(end$gene_name) # N = 1,017


########################################################
# Merge female genes with endothelial cells
merge_f_end = inner_join(end, f1, by = "gene_name") 
n_f_end = n_distinct(merge_f_end$gene_name) # N = 5

# Calculate p-value using hypergeometric test
f_end_p <- phyper(n_f_end - 1, n_end, n_bg - n_end, n_f, lower.tail = FALSE)

# Calculate the GeneRatio
f_end_gr <- n_f_end / n_end

# Calculate the fold enrichment
f_end_fe <- ((n_f_end / n_f) / (n_end / n_bg)) - 1

# Store results in temp df
cell_results = data.frame(
  gene_set = "All Female",
  cell_type = "Endothelial Cells",
  overlap = n_f_end,
  percent_celltype = (n_f_end/n_f)*100,
  gene_ratio = f_end_gr,
  fold_enrichment = f_end_fe,
  p_value = f_end_p)

# Append results to final df
results = rbind(results, cell_results)


########################################################
# Merge male genes with endothelial cells
merge_m_end = inner_join(end, m1, by = "gene_name") 
n_m_end = n_distinct(merge_m_end$gene_name) # N = 2

# Calculate p-value using hypergeometric test
m_end_p <- phyper(n_m_end - 1, n_end, n_bg - n_end, n_m, lower.tail = FALSE)

# Calculate the GeneRatio
m_end_gr <- n_m_end / n_end

# Calculate the fold enrichment
m_end_fe <- ((n_m_end / n_m) / (n_end / n_bg)) - 1

# Store results in temp df
cell_results = data.frame(
  gene_set = "All Male",
  cell_type = "Endothelial Cells",
  overlap = n_m_end,
  percent_celltype = (n_m_end/n_m)*100,
  gene_ratio = m_end_gr,
  fold_enrichment = m_end_fe,
  p_value = m_end_p)

# Append results to final df
results = rbind(results, cell_results)


########################################################################################################################################################################
## Microglia/Macrophage
mic = final_bg %>% 
  filter(cell_type == "Human Microglia/Macrophage") %>% 
  distinct(gene_name)

# Make sure there are no duplicates
n_distinct(mic$gene_name) == nrow(mic)
n_mic = n_distinct(mic$gene_name) # N = 1,789


########################################################
# Merge female genes with Microglia/Macrophages
merge_f_mic = inner_join(mic, f1, by = "gene_name") 
n_f_mic = n_distinct(merge_f_mic$gene_name) # N = 19

# Calculate p-value using hypergeometric test
f_mic_p <- phyper(n_f_mic - 1, n_mic, n_bg - n_mic, n_f, lower.tail = FALSE)

# Calculate the GeneRatio
f_mic_gr <- n_f_mic / n_mic

# Calculate the fold enrichment
f_mic_fe <- ((n_f_mic / n_f) / (n_mic / n_bg)) - 1

# Store results in temp df
cell_results = data.frame(
  gene_set = "All Female",
  cell_type = "Microglia/Macrophages",
  overlap = n_f_mic,
  percent_celltype = (n_f_mic/n_f)*100,
  gene_ratio = f_mic_gr,
  fold_enrichment = f_mic_fe,
  p_value = f_mic_p)

# Append results to final df
results = rbind(results, cell_results)


########################################################
# Merge male genes with Microglia/Macrophages
merge_m_mic = inner_join(mic, m1, by = "gene_name") 
n_m_mic = n_distinct(merge_m_mic$gene_name) # N = 2

# Calculate p-value using hypergeometric test
m_mic_p <- phyper(n_m_mic - 1, n_mic, n_bg - n_mic, n_m, lower.tail = FALSE)

# Calculate the GeneRatio
m_mic_gr <- n_m_mic / n_mic

# Calculate the fold enrichment
m_mic_fe <- ((n_m_mic / n_m) / (n_mic / n_bg)) - 1

# Store results in temp df
cell_results = data.frame(
  gene_set = "All Male",
  cell_type = "Microglia/Macrophages",
  overlap = n_m_mic,
  percent_celltype = (n_m_mic/n_m)*100,
  gene_ratio = m_mic_gr,
  fold_enrichment = m_mic_fe,
  p_value = m_mic_p)

# Append results to final df
results = rbind(results, cell_results)



########################################################################################################################################################################
## Neurons
neu = final_bg %>% 
  filter(cell_type == "Human Neurons") %>% 
  distinct(gene_name)

# Make sure there are no duplicates
n_distinct(neu$gene_name) == nrow(neu)
n_neu = n_distinct(neu$gene_name) # N = 3,272


########################################################
# Merge female genes with Neurons
merge_f_neu = inner_join(neu, f1, by = "gene_name") 
n_f_neu = n_distinct(merge_f_neu$gene_name) # N = 16

# Calculate p-value using hypergeometric test
f_neu_p <- phyper(n_f_neu - 1, n_neu, n_bg - n_neu, n_f, lower.tail = FALSE)

# Calculate the GeneRatio
f_neu_gr <- n_f_neu / n_neu

# Calculate the fold enrichment
f_neu_fe <- ((n_f_neu / n_f) / (n_neu / n_bg)) - 1

# Store results in temp df
cell_results = data.frame(
  gene_set = "All Female",
  cell_type = "Neurons",
  overlap = n_f_neu,
  percent_celltype = (n_f_neu/n_f)*100,
  gene_ratio = f_neu_gr,
  fold_enrichment = f_neu_fe,
  p_value = f_neu_p)

# Append results to final df
results = rbind(results, cell_results)


########################################################
# Merge male genes with Neurons
merge_m_neu = inner_join(neu, m1, by = "gene_name") 
n_m_neu = n_distinct(merge_m_neu$gene_name) # N = 2

# Calculate p-value using hypergeometric test
m_neu_p <- phyper(n_m_neu - 1, n_neu, n_bg - n_neu, n_m, lower.tail = FALSE)

# Calculate the GeneRatio
m_neu_gr <- n_m_neu / n_neu

# Calculate the fold enrichment
m_neu_fe <- ((n_m_neu / n_m) / (n_neu / n_bg)) - 1

# Store results in temp df
cell_results = data.frame(
  gene_set = "All Male",
  cell_type = "Neurons",
  overlap = n_m_neu,
  percent_celltype = (n_m_neu/n_m)*100,
  gene_ratio = m_neu_gr,
  fold_enrichment = m_neu_fe,
  p_value = m_neu_p)

# Append results to final df
results = rbind(results, cell_results)



########################################################################################################################################################################
## Oligodendrocytes
oli = final_bg %>% 
  filter(cell_type == "Human Oligodendrocytes") %>% 
  distinct(gene_name)

# Make sure there are no duplicates
n_distinct(oli$gene_name) == nrow(oli)
n_oli = n_distinct(oli$gene_name) # N = 749


########################################################
# Merge female genes with Oligodendrocytes
merge_f_oli = inner_join(oli, f1, by = "gene_name") 
n_f_oli = n_distinct(merge_f_oli$gene_name) # N = 1

# Calculate p-value using hypergeometric test
f_oli_p <- phyper(n_f_oli - 1, n_oli, n_bg - n_oli, n_f, lower.tail = FALSE)

# Calculate the GeneRatio
f_oli_gr <- n_f_oli / n_oli

# Calculate the fold enrichment
f_oli_fe <- ((n_f_oli / n_f) / (n_oli / n_bg)) - 1

# Store results in temp df
cell_results = data.frame(
  gene_set = "All Female",
  cell_type = "Oligodendrocytes",
  overlap = n_f_oli,
  percent_celltype = (n_f_oli/n_f)*100,
  gene_ratio = f_oli_gr,
  fold_enrichment = f_oli_fe,
  p_value = f_oli_p)

# Append results to final df
results = rbind(results, cell_results)


########################################################
# Merge male genes with Oligodendrocytes
merge_m_oli = inner_join(oli, m1, by = "gene_name") 
n_m_oli = n_distinct(merge_m_oli$gene_name) # N = 0

# Calculate p-value using hypergeometric test
m_oli_p <- phyper(n_m_oli - 1, n_oli, n_bg - n_oli, n_m, lower.tail = FALSE)

# Calculate the GeneRatio
m_oli_gr <- n_m_oli / n_oli

# Calculate the fold enrichment
m_oli_fe <- ((n_m_oli / n_m) / (n_oli / n_bg)) - 1

# Store results in temp df
cell_results = data.frame(
  gene_set = "All Male",
  cell_type = "Oligodendrocytes",
  overlap = n_m_oli,
  percent_celltype = (n_m_oli/n_m)*100,
  gene_ratio = m_oli_gr,
  fold_enrichment = m_oli_fe,
  p_value = m_oli_p)

# Append results to final df
results = rbind(results, cell_results) %>% 
  arrange(gene_set, cell_type)


## Write out results
fwrite(results, file.path(opt$work_dir, "results", opt$results_out))


### End of code



