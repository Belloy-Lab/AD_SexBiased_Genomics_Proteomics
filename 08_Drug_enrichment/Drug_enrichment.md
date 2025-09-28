**AD Sex-Biased Genomics & Proteomics**

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

# Select only FDA approved drugs
Rscript analysis_codes/3_drug_enrichment_result_filtering.R \
    --work_dir working/directory \
    --female_in DsigDB_drug_enrich_res_female_expanded.xlsx \
    --male_in DsigDB_drug_enrich_res_male_expanded.xlsx \
    --FDA_gmt input_files/FDA_All.gmt \
    --female_out FDA_drug_enrich_res_female_filtered.csv \
    --male_out FDA_drug_enrich_res_male_filtered.csv
 
# Make gene-drug network.
Rscript analysis_codes/4_hormone_drug_network.R \
    --work_dir working/directory \
    --hormone_csv Classified_Sex_Hormone_Drugs.csv \
    --female_input FDA_drug_enrich_res_female_filtered.csv \
    --ppi_df_female ppi_df_female_list.csv \
    --female_out hormone_drugs_female_filtered.csv \
    --network_out ppi_network_input_for_cytoscape.csv \
    --node_out ppi_node_type_for_cytoscape.csv
```

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)