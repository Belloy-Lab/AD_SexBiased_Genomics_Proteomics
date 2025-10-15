**AD Sex-Biased Genomics & Proteomics**

## Drug enrichment analysis.
We conducted drug enrichment analyses using sex-specific gene lists and the networks of druggable genes were expanded using EpiGraphDB protein–protein interactions. Only Tier 1 druggable genes ([Finan et al. 2017](https://www.science.org/doi/10.1126/scitranslmed.aag1166?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)) with Alzheimer’s disease associations were retained, generating final female and male gene sets.

Drug enrichment was then performed with DSigDB via clusterProfiler[clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html).

```bash
Rscript Analysis_codes/drug_enrichment.R
```

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)