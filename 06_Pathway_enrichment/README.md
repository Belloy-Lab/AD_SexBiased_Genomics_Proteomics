
# Pathway Enrichment Analysis

This repository contains scripts for performing and visualizing Gene Ontology (GO) pathway enrichment analysis in a sex-stratified manner. The workflow includes clustering similar GO terms and generating figures and word clouds to summarize biological signals from female and male datasets.

## Overview of Scripts

### 1. `Run_pathway_enrichment.R`

Performs GO enrichment analysis for each sex-specific gene set and saves filtered enrichment results.

### 2. `Make_clusters.R`

Clusters GO terms by semantic similarity using `rrvgo` and `simplifyEnrichment`. Generates dendrograms and heatmaps for both sexes.

### 3. `Enrichment_visualization.R`

Generates boxplots, violin plots, and word clouds for each GO cluster from female and male data. Visualization helps identify representative functional groups.

## Dependencies

Ensure the following R packages are installed:

* `optparse`, `data.table`, `clusterProfiler`, `rrvgo`, `org.Hs.eg.db`, `simplifyEnrichment`
* `ggplot2`, `ggraph`, `igraph`, `tm`, `wordcloud`, `stringr`, `RColorBrewer`

To install a package from Bioconductor (e.g., `simplifyEnrichment`):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("simplifyEnrichment")
```
---
**Citation:** If you use these scripts, please cite our PWAS paper (in preparation).  
**License:** MIT (see [main repository README](../README.md) for full text).
