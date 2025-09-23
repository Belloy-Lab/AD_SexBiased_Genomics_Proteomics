**AD sex-biased Genomics & Proteomics**

## Project Directory Structure for Analysis

This project involves GWAS summary statistics from multiple datasets, including **ADSP**, **ADGC**, **UK Biobank**, and **FinnGen**. To streamline analysis and keep outputs organized, we define a structured folder hierarchy for each analysis scenario.

The directory structure includes:
- EU_all-analysis (all datasets combined)
- EU_noUKB-analysis excluding UK Biobank from the EU_all
- Analysis focused on African-admixed (AFR) datasets
- Sex-stratified and non-stratified results
- Dedicated folders for scripts used across different datasets

### Set Working Directory
Navigate to your project directory (choose the appropriate path based on your system).

```bash
mkdir -p PWAS
mkdir -p PWAS/EU_all/logs
mkdir -p PWAS/EU_noUKB/logs
mkdir -p PWAS/AFR/logs
```

### Folder for PWAS figures and tables
Please create a separate folder to create figures and tables. 
```bash
mkdir -p PWAS/Figure_Tables/logs
```

### Set Up Code Directory
Please make sure to set the correct path in all R and Bash scripts located in the analysis_codes directory.

```bash
mkdir -p PWAS/analysis_codes
```

Please copy all scripts from the analysis_codes directory into their corresponding target folders.

### You're All Set
After running the above commands, you’ll have a clean, organized folder layout ready for:  
- Input data  
- Intermediate outputs  
- Final results  
- Project-specific scripts  
This structure ensures reproducibility, clarity, and ease of collaboration across datasets and analysis stages.

**Citation:** If you use these scripts, please cite our PWAS paper (in preparation).  
**License:** MIT (see [main repository README](../README.md) for full text).
