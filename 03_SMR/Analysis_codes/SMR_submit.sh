#!/bin/bash

# Define variables
SMR_BIN="/storage2/fs1/belloy2/Active/04_Code/$USER/SMR/smr"
BFILE="/storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_reference_genotype_data/EU_TOPMed/TOPMed_merged"
PQTL_BESD="/storage2/fs1/belloy2/Active/05_Projects/$USER/SMR/CSF_Brain_pQTL_besd"
OUT_DIR="/storage2/fs1/belloy2/Active/05_Projects/$USER/SMR"
THREADS=10

# Make SMR binary executable
chmod +x $SMR_BIN

# Run SMR Female
$SMR_BIN --bfile $BFILE \
  --gwas-summary $OUT_DIR/AD_female.ma \
  --beqtl-summary $PQTL_BESD \
  --peqtl-smr 1e-3 \
  --out $OUT_DIR/smr_res_female \
  --thread-num $THREADS

# Run SMR Male
$SMR_BIN --bfile $BFILE \
  --gwas-summary $OUT_DIR/AD_male.ma \
  --beqtl-summary $PQTL_BESD \
  --peqtl-smr 1e-3 \
  --out $OUT_DIR/smr_res_male \
  --thread-num $THREADS