**AD Sex-Biased Genomics & Proteomics**

# Cell enrichment analysis.

This document outlines the step-by-step command-line usage of the cell enrichment analysis.

## Introductions
We assessed the cellular specificity of sex-specific prioritized genes using the method of Western et al. (2024). Cell type–specific expression data for endothelial cells, oligodendrocytes, astrocytes, neurons, and microglia/macrophages were obtained from published resources (Zhang et al., 2016). For each gene, mean expression per cell type was calculated, then proportions across the five cell types were derived. A gene was classified as cell type–specific if expression in the top cell type was ≥1.5-fold higher than in the second.

Sex-specific gene sets were mapped to these cell type assignments. For each sex, we counted genes specific to each cell type and tested for enrichment using the hypergeometric test (phyper in R). Only cell types with at least five genes across both sexes were considered. This analysis identifies whether sex-prioritized Alzheimer’s disease genes are enriched in particular brain cell types.

## Analysis
### Step 1 — Cell-specific enrichment
```bash
Rscript analysis_codes/1_cell_specific_analysis.R \
    --work_dir working/directory \
    --gene_list GWAS_PWAS_final_enrichment_list.csv \
    --background_list background_genes.csv \
    --results_out cell_specific_enrichment_results_All.csv
```

### Step 2 — Cell-specific figures
```bash
Rscript analysis_codes/2_cell_specific_figures.R \
    --work_dir working/directory \
    --results_in cell_specific_enrichment_results_All.csv \
    --out_fig cell-type_enrichment_barplot_Filtered.jpg
```

![**Figure.**:Bar plot](results/cell-type_enrichment_barplot_Filtered.jpg)

---
**Citation:** If you use these scripts, please cite our PWAS paper (in preparation).  
**License:** MIT (see [main repository README](../README.md) for full text).