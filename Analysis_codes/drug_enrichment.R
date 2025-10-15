# Drug Repurposing Enrichment analyses base script

# Libraries
library(data.table)
library(tidyverse)
library(clusterProfiler)
library(epigraphdb)
library(HGNChelper)


###############################################################################################################################
## Get expanded gene list

# Read in gene list 
gene_list = fread("/path/to/gene/list")


# Read in druggable tiers file from epigraphdb
druggable_tiers <- readxl::read_excel("aag1166_Table S1.xlsx")
druggable_tiers <- as.data.frame(druggable_tiers) %>% 
  dplyr::select(hgnc_names,druggability_tier)

# Check Gene Symbols in Druggability tiers file using HGNC helper
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


# Merge updated symbols to druggability terms df
druggable_tiers2 <- inner_join(druggable_tiers ,check_gene_dt ,by='hgnc_names') %>% 
  drop_na(Suggested.Symbol) %>%
  dplyr::select(Suggested.Symbol, druggability_tier) %>% 
  dplyr::rename(hgnc_names = Suggested.Symbol)


# find druggable info for raw gene list
candidate_gene_tiers <- data.frame(hgnc_names=gene_male$hgnc_names)
candidate_gene_tiers <- left_join(candidate_gene_tiers,druggable_tiers2)


## find related genes (ppi) + druggable info
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


## Get disese information
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


## merge canidate gene tiers with ppi information
res <- left_join(candidate_gene_tiers,ppi_df_all,by='hgnc_names') 

# Filter to only inclde AD related and druggability tier 1
res_filter = res %>% 
  filter(disease_related == "Yes" & g2.druggability_tier == "Tier 1") %>% 
  dplyr::select(`g2.name`) %>% 
  dplyr::rename(hgnc_names = `g2.name`) 


# Add expanded genes to raw gene list - get unique gene_names
expanded_gene_list = rbind(gene_list, res_filter) %>% 
  distinct(hgnc_names) %>% 
  dplyr::rename(gene_name = hgnc_names)


###############################################################################################################################
## Run enrichment analyses for drug repurposing

# import gmt file downloaded from DSigDB
gmt_df_DSigDB <- read.gmt("DSigDB_All.gmt")
gmt_df_DSigDB$term <- sub("_.*", "", gmt_df_DSigDB$term) # Merge identical drugs from different databases by removing suffixes


# enrichment analysis
enrich_re <- enricher(
  gene = expanded_gene_list$gene_name,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 1,
  maxGSSize = 500,
  qvalueCutoff = 1,
  gson = NULL,
  TERM2GENE = gmt_df_DSigDB, 
  TERM2NAME = NA
)
