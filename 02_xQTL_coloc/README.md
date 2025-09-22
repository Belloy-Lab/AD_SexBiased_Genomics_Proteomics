# xQTL Colocalization Pipeline
The aim of this is to utilize a wide range of QTL datasets to run QTL colocalization analyses with a target dataset. In the example codes provided, the target dataset is a non-stratified GWAS of Alzheimer's disease.

## Requirements
* Bash
* R version 4.4 or later
* R Packages:
   * data.table
   * tidyverse
   * locuscomparer
   * coloc
   * stringr
   * foreach
 * Target dataset (GWAS summary stats for specific phenotype)
 * QTL datasets - see QTL_datasets.md for list of possible QTL datasets and tissues

 
## Workflow
* Create input reference file for the bash scipt with the following variables:
  1. study - Study from which QTL dataset comes from
  2. tissue - specific tissue/cell/brain region the QTL dataset investigates
  3. filepath - path to QTL dataset on local or remote server
  4. qtl_type - what QTL type you are investigating (eQTL, pQTL, sQTL, mQTL, haQTL, caQTL)
     
* Create gene reference file to loop through in R script
  1. chrom - chromosome number the locus of interest falls on
  2. bp_38_start - starting basepair position for larger window approach (GRch38)
  3. bp_38_end - ending basepair position for larger window approach (GRch38)
  4. bp_38_med - middle base pair position for smaller window approach (GRch38)
  5. locus_index - locus of interest for QTL colocalization analysis
  6. stratum - indicate which stratum oof the target dataset that will be used
  7. discovery - (optional) used to label output data and locus compare plots
 

