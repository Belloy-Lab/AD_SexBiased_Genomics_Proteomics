**AD Sex-Biased Genomics & Proteomics**

## Locus zoom plots
Multi-panel LocusZoom plots were generated to visualize the genetic architecture around PWAS-prioritized genes. Two plot types are supported:

1. Single-tag plots showing AD GWAS and pQTL signals for one lead variant.

2. Multi-tag plots integrating AD GWAS, pQTL, and additional QTL datasets for up to three lead variants.

Plots use LD data from AD GWAS Stage 1, but can be adapted to 1000 Genomes or other LD references.


### Scripts for Supplementary Figure S14

- locus_zoom_plot.R – main driver

- functions_locus_zoom_plot.R – helper functions for windowing, merging, LD, and plotting

```bash
Rscript Analysis_codes/locus_zoom_plot.R
```

### Scripts for Figure 4 (Overlay Plots)

- overlay_plot_CD84.R – CD84 locus overlay

- overlay_plot_PSEN1.R – PSEN1 locus overlay

- functions_overlay_plot.R – shared utilities

```bash
Rscript Analysis_codes/overlay_plot_CD84.R

Rscript Analysis_codes/overlay_plot_PSEN1.R
```

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)