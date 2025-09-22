
**Sex Strat AD PWAS**
---
Title: Integrated multi-omics analyses unravel sex-biased causal genes and drug targets in Alzheimer’s disease.
---

## Background
To elucidate sex differences in Alzheimer’s disease (AD), we conducted the largest-to-date sex-stratified genome-wide association study (GWAS) of AD. To further increase power and identify sex-specific, potentially druggable AD causal proteins, we performed proteome-wide association studies (PWAS) integrating GWAS with proteogenomic (i.e., protein quantitative trait locus [pQTL]) brain and cerebrospinal fluid (CSF) datasets.

## Data availability
You can download the GWAS summary statistics from the following links:  
| Dataset Name              | Download Link                             |
|--------------------------|-------------------------------------------|
| ADSP_ADGC_UKB_FinnGen    | [Download](https://your-download-link.com) |
| ADSP_ADGC_FinnGen        | [Download](https://your-download-link.com) |
| AFRad                    | [Download](https://your-download-link.com) |

Brain and CSF PWAS weights were derived from following references
 1. Wingo, A. P. et al. Sex differences in brain protein expression and disease. Nat Med 29, 2224–2232 (2023).
 2. Western, D. et al. Proteogenomic analysis of human cerebrospinal fluid identifies neurologically relevant regulation and implicates causal proteins for Alzheimer’s disease. Nat Genet 56, 2672–2684 (2024).  

You can download the PWAS weights from the following links:  
| Tissue | Download Link                             |
|--------|-------------------------------------------|
| Brain  | [Download](https://your-download-link.com) |
| CSF    | [Download](https://your-download-link.com) |

## Docker Images
This project utilizes the following Docker images for various data processing tasks:

1. dmr07083/fusion-project:4.3  
All data processing steps—including GWAS summary statistics liftover, LD reference panel intersection, GWAS/PWAS association using the FUSION script, and post-analysis processing—were performed using the fusion-project Docker image developed by our lab.

You can download and run the image with the following commands:
```bash
docker pull dmr07083/fusion-project:4.3
docker run -i -t dmr07083/fusion-project:4.3 /bin/bash
```
For more details about this image and its usage, please refer to the Fusion_project.md file.

2. continuumio/anaconda:latest  
Some additional processing steps, such as mungestat cleaning using a Conda environment, were performed using the publicly available Anaconda Docker image.

To use this image, run:
```bash
docker pull continuumio/anaconda:latest
docker run -i -t continuumio/anaconda:latest /bin/bash
```

## 

