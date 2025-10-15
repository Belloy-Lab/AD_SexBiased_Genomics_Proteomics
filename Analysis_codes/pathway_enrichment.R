## Pathway enrichment (Gene Ontology) base script

# Libraries - Installrequired packages via BioConductor
library(data.table)
BiocManager::install("clusterProfiler"); library(clusterProfiler)
BiocManager::install("org.Hs.eg.db"); library(org.Hs.eg.db)


######################################################
## Read in gene list 
gene_list = fread("/path/to/gene/list")

# Convert gene lists into ENTREZ IDs
gene_list = bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb ="org.Hs.eg.db")
genes_entrez = gene_list$ENTREZID


######################################################
## Input raw gene list into enrichGO

# GO overrepresentation - raw list
go_fa = enrichGO(
  gene = genes_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = 'ENTREZID',
  ont = 'BP',
  pvalueCutoff = 1,
  pAdjustMethod = 'BH',
  qvalueCutoff = 1,
  readable = TRUE
)