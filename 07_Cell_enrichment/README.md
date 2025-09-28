**AD Sex-Biased Genomics & Proteomics**

# Cell enrichment analysis.
This document outlines the step-by-step command-line usage of the cell enrichment analysis.

## Introductions
We assessed the cellular specificity of sex-specific prioritized genes using the method of Western et al. (2024). Cell type–specific expression data for endothelial cells, oligodendrocytes, astrocytes, neurons, and microglia/macrophages were obtained from published resources (Zhang et al., 2016). For each gene, mean expression per cell type was calculated, then proportions across the five cell types were derived. A gene was classified as cell type–specific if expression in the top cell type was ≥1.5-fold higher than in the second.

Sex-specific gene sets were mapped to these cell type assignments. For each sex, we counted genes specific to each cell type and tested for enrichment using the hypergeometric test (phyper in R). Only cell types with at least five genes across both sexes were considered. This analysis identifies whether sex-prioritized Alzheimer’s disease genes are enriched in particular brain cell types.

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)