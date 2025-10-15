**AD Sex-Biased Genomics & Proteomics**

# 🧬 Multi-omics toolkit Docker Image (v0.1)

This Docker image provides a pre-configured environment for running our PWAS, xQTL and other downstream multi-omics analyses provided in this repository.

It includes:
- R 4.5 with key CRAN + Bioconductor packages
- Python 3 (with `ldetect` in virtualenv)
- GWAS tools: PLINK, GEMMA, MR-MEGA
- Fully reproducible setup

---

## 🐋 Build the Image Locally

To build this image yourself:

```bash
git clone https://github.com/Belloy-Lab/AD_SexBiased_Genomics_Proteomics
docker build -t multiomics-toolkit:0.1 .
```
## Pull Pre-Built Image
Alternatively, you can download and run the pre-built Docker image from Docker Hub:

```bash
docker pull satheshsiva27/multiomics-toolkit:0.1

```

## Run the Container
Launch an interactive session:

```bash
docker run -it multiomics-toolkit:0.1 bash
#or 
docker run -it satheshsiva27/multiomics-toolkit:0.1 /bin/bash
```
Or mount a working directory:

```bash
docker run -it -v /your/data:/data multiomics-toolkit:0.1 bash
#or 
docker run -it -v /your/local/path:/data satheshsiva27/multiomics-toolkit:0.1 bash
```

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)
