**AD Sex-Biased Genomics & Proteomics**

## xQTL Colocalization 
The aim of this pipeline is to utilize a wide range of QTL datasets to run QTL colocalization analyses with a target dataset. Colocalization is conducted using the [Coloc](https://cran.r-project.org/web/packages/coloc/index.html) CRAN package.

### Datasets 
The datasets included in the xQTL pipeline was described in supplementary table S23.

### Coloc analysis
```bash
# To run ABF colocalization analysis
Rscript Analysis_codes/coloc_abf.R

# To run SuSiE colocalization analysis
Rscript Analysis_codes/coloc_susie.R
```

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)
