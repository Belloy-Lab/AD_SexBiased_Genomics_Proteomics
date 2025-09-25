#!/bin/bash

# Define variables
SMR_BIN="smr"
BFILE="/EUR_1000G_ref_panel_inPLINK_format_all_chr_merged/including/full/path"
PQTL_BESD="/Path/where/pQTL/output/saved/CSF_Brain_merged_pQTL_besd"
OUT_DIR="/Path/to/save/output/files"
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