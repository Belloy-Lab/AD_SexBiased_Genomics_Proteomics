**AD Sex-Biased Genomics & Proteomics**

## xQTL Colocalization Analysis

---

### 1. Job Submission Scripts
The analysis is initiated by submitting jobs via the Bash scripts located in the folder:

```bash
analysis_codes/Job_scripts/
```

This folder contains four Bash scripts:

* xQTL_abf.sh – Submits jobs for colocalization using the ABF method for 159 xQTL datasets.

* xQTL_susie.sh – Submits jobs for colocalization using the SUSiE method for the same 159 xQTL datasets.

* other_xQTL_abf.sh – Submits jobs for colocalization using ABF for 10 additional xQTL datasets that require a different processing workflow.

* other_xQTL_susie.sh – Same as above, but using the SUSiE method.

⚠️ IMPORTANT: Before running these scripts, please update the following:

-   xQTL data file names and paths

-   Log and error output paths (i.e., *.out, *.err)

Replace all instances of $USER with your actual username or preferred directory name


### 2. Main R Scripts for Colocalization
The following R scripts, located in 04_code/$USER/xQTL/, are called by the Bash scripts and handle the main colocalization logic:

* xQTL_master_abf.R – Performs ABF-based coloc analysis for 159 xQTL datasets

* xQTL_master_susie.R – Performs SUSiE-based coloc analysis for 159 xQTL datasets

* other_xQTL_master_abf.R – ABF-based analysis for the remaining 10 xQTL datasets

* other_xQTL_master_susie.R – SUSiE-based analysis for the remaining 10 xQTL datasets

Each R script will need the following paths correctly specified:

-   GWAS summary statistics file(s)

-   File listing the genes or SNPs to be used as input for Coloc

Output folder for results and LocusCompare (LC) plots

🔁 Again, please replace all instances of $USER with your actual username or desired folder name across all R and Bash scripts in the repository.

### 3. Create Analysis Folder Structure
This directory structure organizes the results of various xQTL analysis.

Run the following commands to set up the full directory layout:
```bash
mkdir -p logs
mkdir -p ABF
mkdir -p ABF_LCplots
mkdir -p SUSIE
mkdir -p SUSIE_LCplots
```
Please update the directory paths in the job submittion bash scripts.

### 3. To Run the Analysis
Submit the jobs using:
```bash
bash xQTL_abf.sh
bash xQTL_susie.sh
bash other_xQTL_abf.sh
bash other_xQTL_susie.sh
```
These Bash scripts will call the corresponding R scripts, which in turn source all necessary helper scripts to execute the Coloc analysis for the specified combinations of GWAS and xQTL datasets.

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)
