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
Navigate to your project directory (choose the appropriate path based on your system). Please replace '$USER' with your actual username or preferred folder name.

```bash
cd /storage1/fs1/belloy/Active/05_Projects/$USER/
mkdir -p PWAS
mkdir -p PWAS/EU_all/logs
mkdir -p PWAS/EU_noUKB/logs
mkdir -p PWAS/AFR/logs

# OR

cd /storage2/fs1/belloy2/Active/05_Projects/$USER/
mkdir -p PWAS
mkdir -p PWAS/EU_all/logs
mkdir -p PWAS/EU_noUKB/logs
mkdir -p PWAS/AFR/logs
```

### Folder for PWAS figures and tables
Please create a separate folder to create supporting and final figures and tables. 
```bash
mkdir -p Figure_Tables/logs
```

### Set Up Code Directory
Navigate to 04_Code directory and create folder for PWAS project scripts. Please replace '$USER' with your actual username or preferred folder name. Also, make sure to set the correct path in all R and Bash scripts located in the analysis_codes directory.

```bash
cd /storage1/fs1/belloy/Active/04_Code/$USER/
mkdir PWAS
cd PWAS

# OR

cd /storage2/fs1/belloy2/Active/04_Code/$USER/
mkdir PWAS
cd PWAS
```

Please copy all scripts from the analysis_codes directory into their corresponding target folders. After copying, submit the scripts for execution. Please ensure that '$USER' is replaced with your preferred directory name in all R and Bash scripts to ensure proper execution.

### You're All Set
After running the above commands, you’ll have a clean, organized folder layout ready for:  
- Input data  
- Intermediate outputs  
- Final results  
- Project-specific scripts  
This structure ensures reproducibility, clarity, and ease of collaboration across datasets and analysis stages.

## License (MIT)

Copyright (c) 2025 Sathesh K. Sivasankaran
