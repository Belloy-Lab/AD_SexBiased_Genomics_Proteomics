# Drug Repurposing Enrichment analyses base script

# Libraries
library(data.table)
library(epigraphdb)
BiocManager::install("clusterProfiler"); library(clusterProfiler)



######################################################
## Read in input files
# Gene list
gene_list = fread("/path/to/gene/list")

# import gmt file downloaded from DSigDB
gmt_df_DSigDB = read.gmt("/path/to/gmt/file.gmt")



######################################################
## Run enrichment analyses for drug repurposing

enrich_re <- enricher(
  gene = gene_list$gene_name,
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
