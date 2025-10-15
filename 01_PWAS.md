**AD Sex-Biased Genomics & Proteomics**

## Protein Wide Association Study (PWAS)
PWAS were conducted by combining sex-stratified AD GWAS with sex-matched and non-sex-stratified protein-specific variant weights, respectively, in both brain and CSF samples.

### Pre-clean-up
GWAS summary statistics were cleaned up using the “mungestats.py” utility from the LD Score (LDSC) software [https://github.com/bulik/ldsc](https://github.com/bulik/ldsc).

```bash
munge_sumstats.py --sumstats GWAS_summstat_data \
  --keep-maf \
  --maf-min 0.001 \
  --a1 ALLELE1 \
  --a2 ALLELE0 \
  --snp SNP \
  --p P \
  --frq A1FREQ \
  --N sample_size \
  --out GWAS_summstat_data_cleaned
  ```

### PWAS association
PWAS were conducted via FUSION.assoc_test.R scripts from [FUSION](https://github.com/gusevlab/fusion_twas/tree/master) package.

```bash
Rscript FUSION.assoc_test.R \
        --sumstats GWAS_summstat_data_cleaned \
        --weights Protein_weights_data \
        --weights_dir Path/to/protein/weights/ \
        --ref_ld_chr Path/to/LD/reference/panel/file \
        --out /Full/path/to/store/output/file
```

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)