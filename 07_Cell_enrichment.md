**AD Sex-Biased Genomics & Proteomics**

## Cell-Type Enrichment
We assessed cellular specificity of sex-specific prioritized genes using the method of [Western et al. (2024)](https://www.nature.com/articles/s41588-024-01972-8). Cell type–specific expression data for endothelial cells, oligodendrocytes, astrocytes, neurons, and microglia/macrophages were obtained from [Zhang et al.(2016)](https://www.sciencedirect.com/science/article/pii/S0896627315010193?via%3Dihub). For each gene, mean expression per cell type was calculated, and a gene was considered cell type–specific if expression in the top cell type was ≥1.5-fold higher than in the second cell type.

```bash
Rscript Analysis_codes/cell_enrichment.R

```

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)