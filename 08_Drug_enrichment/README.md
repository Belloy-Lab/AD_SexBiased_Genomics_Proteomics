**AD Sex-Biased Genomics & Proteomics**

# Drug enrichment analysis.
We performed drug enrichment analyses using sex-specific gene lists derived from our prioritization framework. Genes with the highest prioritization score (score = 1) were selected separately for females and males. Gene symbols were validated using the HGNChelper R package, and these lists were expanded into networks of druggable genes by integrating protein–protein interaction data from EpiGraphDB. Only genes with Tier 1 druggability (Finan et al., 2017) and confirmed Alzheimer’s disease literature associations in EpiGraphDB were retained, resulting in expanded female and male gene sets.

We then performed drug enrichment using the DSigDB database with the clusterProfiler::enricher function. Enrichment identified thousands of candidate drugs in both sexes, which were subsequently filtered for FDA approval, statistical significance, and sex-specificity. Finally, sex hormone–related drugs identified in the female-specific results were manually curated and visualized within a protein–drug interaction network to highlight hormone-driven pathways relevant to Alzheimer’s disease.

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)