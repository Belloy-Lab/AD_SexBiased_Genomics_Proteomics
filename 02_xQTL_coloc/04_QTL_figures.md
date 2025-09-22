**Sex Strat AD xQTL**

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

## xQTL Figures
### Step 1 — Preprocess inputs for xQTL figures
```bash
bsub -J xqtl_prep -G compute-belloy-t1 -q general -n 1 -M 40000 -R "rusage[mem=40000]" \
  -o /storage2/fs1/belloy2/Active/05_Projects/sivas/xQTL/logs/xqtl_prep.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/sivas/xQTL/logs/xqtl_prep.%J.err \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  Rscript /storage2/fs1/belloy2/Active/04_Codes/sivas/xQTL/figures/01_xQTL_incGWASPWAS_figures_preprocessing.optparse_exact.R \
    --work_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/xQTL/ \
    --abf_dir ABF \
    --adj_abf_dir Adj_ABF \
    --susie_dir SUSIE \
    --out_GP_qtl PWAS_GWAS_qtl_results_annotated.csv \
    --out_pwas PWAS_qtl_results_annotated.csv \
    --out_gwas GWAS_qtl_results_annotated.csv \
    --out_mqtl ref_PWAS_GWAS_top_mQTL_results.csv \
    --out_haqtl ref_PWAS_GWAS_top_haQTL_results.csv \
    --out_caqtl ref_PWAS_GWAS_top_caQTL_results.csv\
```
Tip: --abf_dir and --susie_dir are directories relative to --work_dir. Place your input CSVs there or adjust the names to match your layout.

### Step 2 — Generate xQTL figures
```bash

bsub -J xqtl_fig -G compute-belloy-t1 -q general -n 1 -M 20000 -R "rusage[mem=20000]" \
  -o /storage2/fs1/belloy2/Active/05_Projects/sivas/xQTL/logs/xqtl_fig.%J.out \
  -e /storage2/fs1/belloy2/Active/05_Projects/sivas/xQTL/logs/xqtl_fig.%J.err \
  -a 'docker(dmr07083/fusion-project:4.3.2)' \
  Rscript /storage2/fs1/belloy2/Active/04_Codes/sivas/xQTL/figures/02_xQTL_incGWASPWAS_figures.optparse_exact.R \
    --work_dir /storage2/fs1/belloy2/Active/05_Projects/sivas/xQTL \
    --out_GP_qtl annotated_HP_WDFY1.csv \
    --mqtl_in ref_PWAS_GWAS_top_mQTL_results.csv \
    --haqtl_in ref_PWAS_GWAS_top_haQTL_results.csv \
    --caqtl_in ref_PWAS_GWAS_top_caQTL_results.csv \
    --gnames_in updated_gene_names.csv \
    --index_in index_genes.csv \
    --plot_out Fig3.png

```
![**Figure 1.**: Gene prioritization at novel sex-biased loci through multi-tissue, multi-QTL colocalization and protein differential abundance analyses. ](results/Fig3.png)

## 
