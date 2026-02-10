**AD Sex-Biased Genomics & Proteomics**

## Variant Effect Predictor
Variant Effect Predictor (VEP) analysis was performed on both GWAS amd PWAS loci. Variants were submitted to the Ensembl [VEP](https://useast.ensembl.org/Tools/VEP) web tool.

- Upload the VCF files directly to the Ensembl VEP web tool.
- Use default parameters, with the following modifications:
    - Species: Homo sapiens
    - Under "Additional configurations":
        - Navigate to "Additional annotations" and enable LOEUF.
        - Under "Predictions", enable CADD and REVEL.
- Submit the job and export the annotated results in .txt format.

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)
