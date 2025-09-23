**Sex Strat AD xQTL**

# Drug enrichment analysis.

This document outlines the step-by-step command-line usage of the drug enrichment analysis.

## Introductions
To explore potential drug repurposing opportunities in Alzheimer’s disease, we performed drug enrichment analyses using sex-specific gene lists derived from our prioritization framework. Genes with the highest prioritization score (score = 1) were selected separately for females and males. Gene symbols were validated using the HGNChelper R package, and these lists were expanded into networks of druggable genes by integrating protein–protein interaction data from EpiGraphDB. Only genes with Tier 1 druggability (Finan et al., 2017) and confirmed Alzheimer’s disease literature associations in EpiGraphDB were retained, resulting in expanded female and male gene sets.

We then performed drug enrichment using the DSigDB database with the clusterProfiler::enricher function. Enrichment identified thousands of candidate drugs in both sexes, which were subsequently filtered for FDA approval, statistical significance, and sex-specificity. Finally, sex hormone–related drugs identified in the female-specific results were manually curated and visualized within a protein–drug interaction network to highlight hormone-driven pathways relevant to Alzheimer’s disease.


### Environment Setup (Required Before Running Any Script)

Before launching jobs or interactive sessions, set the following environment variables:

```bash
export LSF_DOCKER_VOLUMES="/storage1/fs1/belloy/Active:/storage1/fs1/belloy/Active \
/storage2/fs1/belloy2/Active:/storage2/fs1/belloy2/Active /scratch1/fs1/belloy:/scratch1/fs1/belloy $HOME:$HOME"

export CONDA_ENVS_DIRS="/storage1/fs1/belloy/Active/conda/envs/"
export CONDA_PKGS_DIRS="/storage1/fs1/belloy/Active/conda/pkgs/"
export PATH="/opt/conda/bin:$PATH"
export LSF_DOCKER_ENTRYPOINT=/bin/bash
```

### Launch interactive session:
Run the interactive section of fusion project docker if it is needed.

```bash
bsub -Is -G compute-belloy-t1 -q subscription -R 'rusage[mem=20GB]' -a 'docker(dmr07083/fusion-project:4.3.2)' /bin/bash
```


### Set Working Directory
Navigate to your project directory (choose the appropriate path based on your system). Please replace '$USER' with your actual username or preferred folder name. Also, make sure to replace '$USER' with your actual username or preferred folder name in all R and Bash scripts located in the analysis_codes directory.

```bash
cd /storage1/fs1/belloy/Active/05_Projects/$USER/
mkdir -p Drug_enrichment/gene_lists
mkdir -p Drug_enrichment/results
mkdir -p Drug_enrichment/network

# OR
cd /storage2/fs1/belloy2/Active/05_Projects/$USER/
mkdir -p Drug_enrichment/gene_lists
mkdir -p Drug_enrichment/results
mkdir -p Drug_enrichment/network
```

### Set Up Code Directory
Navigate to 04_Code directory from storage1 or storage2 and set up project-specific code folders:
```bash
cd /storage1/fs1/belloy/Active/04_Code/$USER/
mkdir Drug_enrichment
# OR
cd /storage2/fs1/belloy2/Active/04_Code/$USER/
mkdir Drug_enrichment
```

Please download the analysis_codes folder from this repository and copy its entire contents into the 04_Code/$USER/Drug_enrichment/ directory on the server.
Similarly, download the curated gene lists from the genelist folder and copy into gene_lists folder.
If any folder paths differ from the default setup, make sure to update them accordingly in all R and Bash scripts.

## Analysis

```bash
# Extend the gene list.
  Rscript /storage1/fs1/belloy/Active/04_Code/$USER/Drug_enrichment/1_get_expanded_gene_list.R \
    --work_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/Drug_enrichment \
    --male_list Male_gene_list3_filter2.txt \
    --female_list Female_gene_list3_filter2.txt \
    --male_out Male_gene_list3_filter2_expanded.txt \
    --female_out Female_gene_list3_filter2_expanded.txt

# Run the drug entichment analysis
Rscript /storage1/fs1/belloy/Active/04_Code/$USER/Drug_enrichment/2_drug_enrichment.R \
    --work_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/Drug_enrichment \
    --dsigdb_gmt DSigDB_All.gmt \
    --male_list Male_gene_list3_filter2_expanded.txt \
    --female_list Female_gene_list3_filter2_expanded.txt \
    --male_out DsigDB_drug_enrich_res_male_expanded.xlsx \
    --female_out DsigDB_drug_enrich_res_female_expanded.xlsx

# Select onkly FDA approved drugs
Rscript /storage1/fs1/belloy/Active/04_Code/$USER/Drug_enrichment/3_drug_enrichment_result_filtering.R \
    --work_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/Drug_enrichment \
    --female_in DsigDB_drug_enrich_res_female_expanded.xlsx \
    --male_in DsigDB_drug_enrich_res_male_expanded.xlsx \
    --FDA_gmt FDA_approved.gmt \
    --female_filter1 FDA_drug_enrich_res_female_filter1.csv \
    --male_filter1 FDA_drug_enrich_res_male_filter1.csv \
    --female_filter2 FDA_drug_enrich_res_female_filter2.csv \
    --male_filter2 FDA_drug_enrich_res_male_filter2.csv \
    --female_filter3 FDA_drug_enrich_res_female_filter3.csv \
    --male_filter3 FDA_drug_enrich_res_male_filter3.csv

# Make gene-drug network.
Rscript /storage1/fs1/belloy/Active/04_Code/$USER/Drug_enrichment/4_hormone_drug_network.R \
    --work_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/Drug_enrichment \
    --hormone_csv Classified_Sex_Hormone_Drugs.csv \
    --female_filter2 FDA_drug_enrich_res_female_filter2_1.5-fold.csv \
    --male_filter2 FDA_drug_enrich_res_male_filter2_1.5-fold.csv \
    --female_out hormone_drugs_female_filter2.csv \
    --male_out hormone_drugs_male_filter2.csv \
    --network_out ppi_network_input_for_cytoscape.csv \
    --node_out ppi_node_type_for_cytoscape.csv
```

---
**Citation:** If you use these scripts, please cite our PWAS paper (in preparation).  
**License:** MIT (see [main repository README](../README.md) for full text).