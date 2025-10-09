# MultiOmics toolkit Docker - Version 0.1
# Base: R 4.5 exact (compatible w/ Bioconductor 3.21)

# Use the latest Debian base image
FROM --platform=linux/amd64 debian:bookworm

# Set environment variables
ENV R_HOME="/usr/lib/R" \
    PATH="${R_HOME}/bin:/usr/bin:/usr/local/bin:$PATH" \
    LANG=en_US.utf8 \
    R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

# Install system dependencies in a single layer
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libudunits2-dev libgdal-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libgit2-dev libglpk-dev libmagick++-dev \
    python3 python3-pip git plink2 gemma less plink1.9 cmake gnupg2 dirmngr \
    ca-certificates wget curl make gcc g++ gfortran locales python3-pandas \
    python3-numpy bcftools \
 && rm -rf /var/lib/apt/lists/*

# editors & text utilities
RUN apt-get update && apt-get install -y --no-install-recommends \
    vim dos2unix \
 && rm -rf /var/lib/apt/lists/*

# Configure R repository
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends gnupg2 dirmngr ca-certificates curl; \
    install -d -m0755 /etc/apt/keyrings; \
    gpg --batch --keyserver keyserver.ubuntu.com \
        --recv-keys 95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7; \
    gpg --batch --export 95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7 \
        | gpg --dearmor -o /etc/apt/keyrings/cran.gpg; \
    echo "deb [signed-by=/etc/apt/keyrings/cran.gpg] http://cloud.r-project.org/bin/linux/debian bookworm-cran40/" \
        > /etc/apt/sources.list.d/cran.list; \
    apt-get update; \
    apt-get install -y --no-install-recommends r-base=4.5*; \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libgsl-dev

#add R to system PATH
ENV R_HOME="/usr/lib/R/bin"
ENV PATH="${R_HOME}:${PATH}"

#set locale for use in RIS
# make the "en_US.UTF-8" locale so postgres will be utf-8 enabled by default
RUN apt-get update && apt-get install -y locales && rm -rf /var/lib/apt/lists/* \
        && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8
ENV LANG en_US.utf8

# remotes first (used below)
RUN R -e 'install.packages("remotes", repos="https://cloud.r-project.org")'

# Bioconductor 3.20 + your Bioc packages
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.21', ask = FALSE)"
RUN R -e 'BiocManager::install(pkgs=c("clusterProfiler","AnnotationDbi","DOSE","enrichplot","GO.db","GOSemSim","qvalue","ReactomePA","graphite","reactome.db","HDO.db","BiocParallel","fgsea","AnnotationHub","org.Hs.eg.db","ensembldb","EnsDb.Hsapiens.v75","sva","genefilter","BiocFileCache","biomaRt","rrvgo"), ask=FALSE, update=FALSE)'
RUN R -e 'BiocManager::install(c("EnsDb.Hsapiens.v86","ComplexHeatmap","simplifyEnrichment"), dependencies=TRUE, update=TRUE, force=TRUE)'

# CRAN set
RUN R -e 'install.packages(c("optparse","RColorBrewer","glmnet","usethis","pkgdown","devtools","readr","ggtext","haven","ggrepel","plyr","openxlsx","igraph","downloader","tidyr","yulab.utils","enrichR","gson","ggraph","reshape2","qqman","tidyverse","CMplot","UpSetR","plinkFile","naniar","ggnewscale","hrbrthemes","susieR","mixsqp","coloc","microbenchmark","Rfast","mvtnorm","pheatmap","wordcloud","treemap","heatmaply"), dependencies=TRUE, repos="https://cloud.r-project.org")'

# env flag (fixed quoting) + plink2R from GitHub
RUN R -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="TRUE"); \
          remotes::install_github("gabraham/plink2R", subdir="plink2R", ref="master")'

# pin ggplot2 3.4.4
RUN R -e 'remotes::install_version("ggplot2", version="3.4.4", repos="https://cloud.r-project.org")'

# more GitHub packages (as you listed)
RUN R -e 'remotes::install_github(c("boxiangliu/locuscomparer","myles-lewis/locuszoomr","MRCIEU/epigraphdb-r"))'
RUN R -e 'devtools::install_github("kassambara/ggpubr")'

# remaining CRAN pkgs from your list
RUN R -e 'install.packages(c("ggtree","writexl","readxl","ggpattern","vroom","data.table"), dependencies=TRUE, repos="https://cloud.r-project.org")'
RUN R -e 'install.packages(c("RcppGSL","RcppZiggurat"), dependencies=TRUE, repos="https://cloud.r-project.org")'
RUN R -e 'install.packages(c("Rfast","HGNChelper"), dependencies=TRUE, repos="https://cloud.r-project.org")'

# final Bioc package from your list
RUN R -e 'BiocManager::install("rtracklayer", dependencies=TRUE, update=FALSE)'

RUN mkdir MR-MEGA \
 && cd MR-MEGA \
 && wget https://tools.gi.ut.ee/tools/MR-MEGA_v0.2.zip \
 && unzip MR-MEGA_v0.2.zip \
 && make \
 && wget https://tools.gi.ut.ee/tools/fixP.r \
 && wget https://tools.gi.ut.ee/tools/manh.r \
 && wget https://tools.gi.ut.ee/tools/qq.r \
 && sed -i 's/read.table/fread/g' fixP.r \
 && sed -i 's/write.table/fwrite/g' fixP.r \
 && sed -i 's/data<-/require(data.table)\\ndata<-/g' fixP.r \
 && cd .. \
 && chmod a+x /MR-MEGA/

#set environment variables for tool paths to the system PATH
ENV PLINK2_HOME="usr/bin/plink2"
ENV GEMMA_HOME="/usr/bin/gemma"
ENV PYTHON_HOME="/usr/bin/python3"

#add tool binary paths to the system PATH
ENV PATH="${PLINK2_HOME}:${GEMMA_HOME}:${PYTHON_HOME}:${PATH}"

#install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda init

# install impute2
RUN mkdir /impute2 && cd /impute2 \
    && wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz \
    && tar -zxvf impute_v2.3.2_x86_64_static.tgz \
    && chmod -R a+x /impute2/impute_v2.3.2_x86_64_static/impute2
ENV IMPUTE2_HOME="/impute2/impute_v2.3.2_x86_64_static"
ENV PATH="${IMPUTE2_HOME}:${PATH}"

#Cleanup
RUN apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
