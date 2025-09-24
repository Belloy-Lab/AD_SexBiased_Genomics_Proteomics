Integrative Genetic, Proteogenomic, and Multi-omics Analyses Reveal Sex-Biased Causal Genes and Drug Targets in Alzheimer’s Disease.
---

## Background
To elucidate sex differences in Alzheimer’s disease (AD), we conducted the largest-to-date sex-stratified genome-wide association study (GWAS) and Proteome-wide association study (PWAS) of AD in brain and cerebrospinal fluid (CSF) proteogenomic datasets. This repository provides codes and files related to these analyses and a series of downstream analyses.

## Data availability
GWAS summary statistics are available from the following Zenodo repository:  
[Download](https://your-download-link.com) 

Brain and CSF PWAS weights were derived from following references:
 1. Wingo, A. P. et al. Sex differences in brain protein expression and disease. Nat Med 29, 2224–2232 (2023).
 2. Western, D. et al. Proteogenomic analysis of human cerebrospinal fluid identifies neurologically relevant regulation and implicates causal proteins for Alzheimer’s disease. Nat Genet 56, 2672–2684 (2024).  

You can download the proteogenomic data from the following links:  
| Tissue | Download Link                                          |
|--------|--------------------------------------------------------|
| Brain  | https://doi.org/10.7303/syn51150434                    |
| CSF    | https://neurogenomics.wustl.edu/open-science/raw-data/ |

## Docker Images
This project utilizes the following Docker images for various data processing tasks:

1. dmr07083/fusion-project:4.3  
Presented analyses were performed using this fusion-project Docker image developed by our lab.

You can download and run the image with the following commands:
```bash
docker pull dmr07083/fusion-project:4.3
docker run -i -t dmr07083/fusion-project:4.3 /bin/bash
```
For more details about this image and its usage, please refer to the Fusion_project.md file.

2. continuumio/anaconda:latest  
Some processing steps, such as mungestat cleaning, were performed using the publicly available Anaconda Docker image.

To use this image, run:
```bash
docker pull continuumio/anaconda:latest
docker run -i -t continuumio/anaconda:latest /bin/bash
```

**Citation:** If you use these scripts, please cite our paper:
[Add preprint citation once available]

## License (MIT)

Copyright (c) 2025 Michael Belloy

Permission is hereby granted, free of charge, to any person obtaining a copy
of these softwares, scripts, pipelines, and associated documentation files (the "Code"),
to deal in the Code without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Code, and to permit persons to whom the Code is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Code.

THE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE CODE OR THE USE OR OTHER DEALINGS IN
THE CODE.

