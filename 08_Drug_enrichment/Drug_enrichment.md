**AD Sex-Biased Genomics & Proteomics**

# Drug enrichment analysis.
This document outlines the step-by-step command-line usage of the drug enrichment analysis.

## Introductions
To explore potential drug repurposing opportunities in Alzheimer’s disease, we performed drug enrichment analyses using sex-specific gene lists derived from our prioritization framework. Genes with the highest prioritization score (score = 1) were selected separately for females and males. Gene symbols were validated using the HGNChelper R package, and these lists were expanded into networks of druggable genes by integrating protein–protein interaction data from EpiGraphDB. Only genes with Tier 1 druggability (Finan et al., 2017) and confirmed Alzheimer’s disease literature associations in EpiGraphDB were retained, resulting in expanded female and male gene sets.

We then performed drug enrichment using the DSigDB database with the clusterProfiler::enricher function. Enrichment identified thousands of candidate drugs in both sexes, which were subsequently filtered for FDA approval, statistical significance, and sex-specificity. Finally, sex hormone–related drugs identified in the female-specific results were manually curated and visualized within a protein–drug interaction network to highlight hormone-driven pathways relevant to Alzheimer’s disease.

## Analysis

```bash
# Extend the gene list.
  Rscript analysis_codes/1_get_expanded_gene_list.R \
    --work_dir working/directory \
    --input_gene_list input_files/GWAS_PWAS_final_enrichment_list.csv \
    --male_out Male_genes_filtered_expanded.txt \
    --female_out Female_genes_filtered_expanded.txt

# Run the drug entichment analysis
Rscript analysis_codes/2_drug_enrichment.R \
    --work_dir working/directory \
    --dsigdb_gmt input_files/DSigDB_All.gmt \
    --male_list Male_genes_filtered_expanded.txt \
    --female_list Female_genes_filtered_expanded.txt \
    --male_out DsigDB_drug_enrich_res_male_expanded.xlsx \
    --female_out DsigDB_drug_enrich_res_female_expanded.xlsx

# Select onky FDA approved drugs
Rscript analysis_codes/3_drug_enrichment_result_filtering.R \
    --work_dir working/directory \
    --female_in DsigDB_drug_enrich_res_female_expanded.xlsx \
    --male_in DsigDB_drug_enrich_res_male_expanded.xlsx \
    --FDA_gmt input_files/FDA_All.gmt \
    --female_filter2 FDA_drug_enrich_res_female_filtered.csv \
    --male_filter2 FDA_drug_enrich_res_male_filtered.csv \
 
# Make gene-drug network.
Rscript analysis_codes/4_hormone_drug_network.R \
    --work_dir working/directory \
    --hormone_csv Classified_Sex_Hormone_Drugs.csv \
    --female_filter2 FDA_drug_enrich_res_female_filtered.csv \
    --male_filter2 FDA_drug_enrich_res_male_filtered.csv \
    --female_out hormone_drugs_female_filtered.csv \
    --male_out hormone_drugs_male_filtered.csv \
    --network_out ppi_network_input_for_cytoscape.csv \
    --node_out ppi_node_type_for_cytoscape.csv
```

---
**Citation:** If you use these scripts, please cite our PWAS paper (in preparation).  
**License:** MIT (see [main repository README](../README.md) for full text).