**Sex Strat AD xQTL**

# Cell enrichment analysis.

This document outlines the step-by-step command-line usage of the cell enrichment analysis.

## Introductions
We assessed the cellular specificity of sex-specific prioritized genes using the method of Western et al. (2024). Cell type–specific expression data for endothelial cells, oligodendrocytes, astrocytes, neurons, and microglia/macrophages were obtained from published resources (
    Zhang et al., 2016). For each gene, mean expression per cell type was calculated, then proportions across the five cell types were derived. A gene was classified as cell type–specific if expression in the top cell type was ≥1.5-fold higher than in the second.

Sex-specific gene sets were mapped to these cell type assignments. For each sex, we counted genes specific to each cell type and tested for enrichment using the hypergeometric test (phyper in R). Only cell types with at least five genes across both sexes were considered. This analysis identifies whether sex-prioritized Alzheimer’s disease genes are enriched in particular brain cell types.

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
mkdir -p Cell_enrichment/gene_lists
mkdir -p Cell_enrichment/results
mkdir -p Cell_enrichment/logs

# OR
cd /storage2/fs1/belloy2/Active/05_Projects/$USER/
mkdir -p Cell_enrichment/gene_lists
mkdir -p Cell_enrichment/results
mkdir -p Cell_enrichment/logs
```

### Set Up Code Directory
Navigate to 04_Code directory from storage1 or storage2 and set up project-specific code folders:
```bash
cd /storage1/fs1/belloy/Active/04_Code/$USER/
mkdir Cell_enrichment
# OR
cd /storage2/fs1/belloy2/Active/04_Code/$USER/
mkdir Cell_enrichment
```

Please download the analysis_codes folder from this repository and copy its entire contents into the 04_Code/$USER/Cell_enrichment/ directory on the server.
Similarly, download the curated gene lists from the genelist folder and copy into gene_lists folder.
If any folder paths differ from the default setup, make sure to update them accordingly in all R and Bash scripts.

## Analysis

### Step 1 — Cell-specific enrichment
```bash
bsub -J cell_step1 -G compute-belloy-t1 -q general -n 1 -M 20000 -R "rusage[mem=20000]" \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/Cell_enrichment/logs/cell_step1.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/Cell_enrichment/logs/cell_step1.%J.err \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  Rscript /storage1/fs1/belloy/Active/04_Code/$USER/Cell_enrichment/1_cell_specific_analysis.R \
    --work_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/Cell_enrichment \
    --gene_list GWAS_PWAS_final_enrichment_list.csv \
    --background_list background_genes.csv \
    --results_out cell_specific_enrichment_results_All.csv
```

### Step 2 — Cell-specific figures
```bash
bsub -J cell_step2 -G compute-belloy-t1 -q general -n 1 -M 16000 -R "rusage[mem=16000]" \
  -o /storage2/fs1/belloy2/Active/05_Projects/$USER/Cell_enrichment/logs/cell_step2.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/$USER/Cell_enrichment/logs/cell_step2.%J.err \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  Rscript /storage1/fs1/belloy/Active/04_Code/$USER/Cell_enrichment/2_cell_specific_figures.R \
    --work_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/Cell_enrichment \
    --results_in cell_specific_enrichment_results_All.csv \
    --plot3_out cell-type_enrichment_barplot_Filtered.jpg

# Option to plot cell plot without any filtering step (add this flag to the same command if you want it):
    --plot5_out cell-type_enrichment_barplot_NoFilter.jpg
```

![**Figure.**:Bar plot](results/cell-type_enrichment_barplot_Filtered.jpg)

---
**Citation:** If you use these scripts, please cite our PWAS paper (in preparation).  
**License:** MIT (see [main repository README](../README.md) for full text).