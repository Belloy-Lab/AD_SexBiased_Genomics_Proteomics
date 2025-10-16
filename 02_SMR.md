**AD Sex-Biased Genomics & Proteomics**

## Summary-based Mendelian Randomization (SMR) Analysis Pipeline
SMR analyses were conducted across GWAS and both Brain and CSF pQTL datasets to validate genes identified in AD PWAS.

Tool Reference: [SMR & HEIDI Software](https://yanglab.westlake.edu.cn/software/smr/#SMR&HEIDIanalysis)

```bash

# Run SMR Female
smr --bfile 1KGP_ref_panel_in_plink_format/including/full/path \
  --gwas-summary Female_GWAS_summstat_in_MA_format/including/full/path \
  --beqtl-summary pQTL_in_BESD_format/with/full/file/path \
  --peqtl-smr 1e-3 \
  --out /Path/to/save/output/files/smr_res_female \

# Run SMR Male
smr --bfile 1KGP_ref_panel_in_plink_format/including/full/path \
  --gwas-summary Male_GWAS_summstat_in_MA_format/including/full/path \
  --beqtl-summary pQTL_in_BESD_format/with/full/file/path \
  --peqtl-smr 1e-3 \
  --out /Path/to/save/output/files/smr_res_male \
  ```


---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)