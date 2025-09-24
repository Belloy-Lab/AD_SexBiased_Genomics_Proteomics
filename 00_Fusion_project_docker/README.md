# 🧬 FUSION Project Docker Image (v4.3.2)

This Docker image provides a pre-configured environment for running the [FUSION](http://gusevlab.org/projects/fusion/) pipeline and other downstream analyses provided in this repository.

It includes:
- R 4.4 with key CRAN + Bioconductor packages
- Python 3 (with `ldetect` in virtualenv)
- GWAS tools: PLINK, GEMMA, MR-MEGA, SHAPEIT2, RFMix, IMPUTE2
- Fully reproducible setup

---

## 🐋 Build the Image Locally

To build this image yourself:

```bash
git clone https://github.com/Belloy-Lab/Sex_Strat_AD_GWAS_PWAS/fusion_project_docker
cd Sex_Strat_AD_GWAS_PWAS/fusion_project_docker
docker build -t fusion-project:4.3.2 .
```
## Pull Pre-Built Image
Alternatively, you can download and run the pre-built Docker image from Docker Hub:

```bash
docker pull belloylab/fusion-project:4.3.2

```

## Run the Container
Launch an interactive session:

```bash
docker run -it fusion-project:4.3.2 bash
#or 
docker run -it belloylab/fusion-project:4.3.2 /bin/bash
```
Or mount a working directory:

```bash
docker run -it -v /your/data:/data fusion-project:4.3.2 bash
#or 
docker run -it -v /your/local/path:/data belloylab/fusion-project:4.3.2 bash
```

**Citation:** If you use these scripts, please cite our PWAS paper (in preparation).  
**License:** MIT (see [main repository README](../README.md) for full text).
