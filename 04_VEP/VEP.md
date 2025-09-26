**AD Sex-Biased Genomics & Proteomics**

## Variant Effect Predictor

Variant Effect Predictor (VEP) analysis was performed on both GWAS amd PWAS loci. Variants must be preprocessed into standard VCF format before being submitted to the Ensembl VEP web tool. Both GWAS and PWAS variants are given in the Variants_list folder in both VCF and CSV format for VEP analysis.

# split here into readme.me and analysis.md

**AD Sex-Biased Genomics & Proteomics**

## Variant Effect Predictor Analysis

- Upload the VCF files directly to the Ensembl VEP web (https://useast.ensembl.org/Tools/VEP) tool.
- Use default parameters, with the following modifications:
    - Species: Homo sapiens
    - Under "Additional configurations":
        - Navigate to "Additional annotations" and enable LOEUF.
        - Under "Predictions", enable CADD and REVEL.
- Submit the job and export the annotated results in .txt format in the "variants_list" folder.

Both resulting .txt files from the VEP analysis (GWAS and PWAS) will be processed using the R script provided below.

### Post-processing

Matrixes show VEP-annotated variant consequences for top GWAS and PWAS loci. Columns represent loci (gene names), and rows show functional consequence types. Red blocks indicate at least one SNP at that locus with the corresponding annotation. The plot highlights the range of functional impacts, including both coding and non-coding variants.

```bash
Rscript analysis_codes/VEP_figures.R \
	--VEP_Dir variants_list/ \
	--Results_dir Results/ \
    --GWAS_VEP  GWAS_VEP_output.txt \
    --GWAS_Var gwas_final_variants.csv \
    --PWAS_VEP PWAS_VEP_output.txt \
    --PWAS_var pwas_final_variants.csv
```

![**Figure 1.**:SNP Consequence Matrix by Locus](Results/PWAS_VEP_matrix_1a.png)

![**Figure 2.**:Plot](Results/PWAS_VEP_matrix_2a.png)

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)
