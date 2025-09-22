
# Enrichment Analysis Workflow

This document outlines the step-by-step command-line usage of the enrichment analysis pipeline, consisting of three modular R scripts:

---

## 📌 Environment Setup (Required Before Running Any Script)

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

```bash
bsub -Is -G compute-belloy-t1 -q subscription -R 'rusage[mem=20GB]' -a 'docker(dmr07083/fusion-project:4.3.2)' /bin/bash
```
---

## Project Directory Structure for Analysis
### Set Working Directory
Navigate to your project directory (choose the appropriate path based on your system). Please replace '$USER' with your actual username or preferred folder name. Also, make sure to replace '$USER' with your actual username or preferred folder name in all R and Bash scripts located in the analysis_codes directory.

```bash
cd /storage1/fs1/belloy/Active/05_Projects/$USER/
mkdir Pathway
# OR
cd /storage2/fs1/belloy2/Active/05_Projects/$USER/
mkdir Pathway
```

### Set Up Code Directory
Navigate to 04_Code directory from storage1 or storage2 and set up project-specific code folders:
```bash
cd /storage1/fs1/belloy/Active/04_Code/$USER/
mkdir Pathway
# OR
cd /storage2/fs1/belloy2/Active/04_Code/$USER/
mkdir Pathway
```

Please download the analysis_codes folder from this repository and copy its entire contents into the 04_Code/$USER/Pathway/ directory on the server.
If any folder paths differ from the default setup, make sure to update them accordingly in all R and Bash scripts.

---
## 1. Functional Enrichment Analysis

Run the main enrichment script to perform GO term enrichment for both female and male gene sets.

```bash
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/Pathway/Run_pathway_enrichment.R \
  --work_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/Pathway \
  --gene_list GWAS_PWAS_final_enrichment_list.csv \
  --female_out Female_GO_output.csv \
  --male_out Male_GO_output.csv \
  --female_filtered Female_GO_filtered.csv \
  --male_filtered Male_GO_filtered.csv
```

## 2. Clustering of Enriched GO Terms

Generate clusters based on semantic similarity of GO terms and output visualizations.

```bash
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/Pathway/Enrichment_clusters.R \
  --work_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/Pathway \
  --female_filtered Female_GO_filtered.csv \
  --male_filtered Male_GO_filtered.csv \
  --female_out_csv Female_clusters_final.csv \
  --male_out_csv Male_clusters_final.csv \
  --female_dendro_pdf female_dendrogram.pdf \
  --male_dendro_pdf male_dendrogram.pdf \
  --female_heatmap_pdf female_heatmap.pdf \
  --male_heatmap_pdf "male_heatmap.pdf
```

## 3. Visualization of Clustered Results

Generate violin/box plots and word clouds for each enriched cluster.

```bash
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/Pathway/Enrichment_visualization.R \
--work_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/Pathway \
--female_input Female_clusters_final.csv \
--male_input Male_clusters_final.csv
```

## Plots and heatmaps

![**Figure 1.**: Violin plot showing enriched pathways on sex-biased prioritized AD genes - female. ](figures/Female_raw_cluster_pathways_violin.jpg)


![**Figure 1.**: Violin plot showing enriched pathways on sex-biased prioritized AD genes - male. ](figures/Male_raw_cluster_pathways_violin.jpg)


![**Figure 1.**: Heatmap showing enriched pathways - female. ](figures/female_heatmap.jpg)


![**Figure 1.**: Heatmap showing enriched pathways - female. ](figures/male_heatmap.jpg)

---

This modular pipeline supports scalable and reproducible enrichment analysis for sex-stratified gene sets.