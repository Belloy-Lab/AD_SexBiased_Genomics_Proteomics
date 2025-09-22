# FUSION Project Docker - Version 4.3.2
# Includes: R 4.4, Python (with ldetect), PLINK2, GEMMA, MR-MEGA, SHAPEIT2, RFMix, IMPUTE2, and Tractor Conda env

FROM --platform=linux/amd64 debian:bookworm

# Install system dependencies, R build tools, and core bioinformatics tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libudunits2-dev libgdal-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libgit2-dev libglpk-dev libmagick++-dev \
    libssl-dev openssl apt-utils libgsl-dev cmake python3 python3-pip \
    python3-full git wget plink2 gemma software-properties-common less \
    plink1.9 gwama tabix samtools libhts3 libncurses6 cwltool libhtscodecs2 \
    bcftools vcftools libgcc-s1 picard-tools default-jre libpicard-java \
    bwa bolt-lmm bedtools libboost-iostreams1.74.0 libboost-program-options1.74.0 \
    libnlopt0 libnode-dev autoconf automake libtool build-essential \
    locales gnupg2 dirmngr ca-certificates apt-transport-https && \
    apt-get clean

# Set locale for UTF-8 support
RUN apt-get install -y locales && rm -rf /var/lib/apt/lists/* && \
    localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8
ENV LANG=en_US.utf8

# Python virtual environment setup and ldetect install
RUN python3 -m venv /opt/venv && \
    /opt/venv/bin/pip install ldetect
ENV PYTHON_VENV="/opt/venv/bin"
ENV PATH="${PYTHON_VENV}:${PATH}"

# Add CRAN and install R 4.4
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7 && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 0xB8F25A8A73EACF41 && \
    echo 'deb http://cloud.r-project.org/bin/linux/debian bookworm-cran40/' > /etc/apt/sources.list.d/cran.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends r-base=4.4* make gcc g++ gfortran

ENV R_HOME="/usr/lib/R/bin"
ENV PATH="${R_HOME}:${PATH}"

# Install R packages
RUN R -e "install.packages('remotes', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c( \
    'optparse','RColorBrewer','glmnet','usethis','pkgdown','devtools','readr','ggtext','haven', \
    'ggrepel','plyr','openxlsx','igraph','downloader','tidyr','yulab.utils','enrichR','gson', \
    'ggraph','reshape2','qqman','tidyverse','CMplot','UpSetR','plinkFile','naniar','ggnewscale', \
    'hrbrthemes','susieR','mixsqp','coloc','microbenchmark','Rfast','mvtnorm','pheatmap', \
    'wordcloud','treemap','heatmaply'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "remotes::install_github('gabraham/plink2R', subdir='plink2R', ref='master')" && \
    R -e "remotes::install_version('ggplot2', version='3.4.4', repos='http://cran.rstudio.com/')"

# Install Bioconductor packages
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')" && \
    R -e "BiocManager::install(version = '3.19', ask = FALSE)" && \
    R -e "BiocManager::install(pkgs=c( \
        'clusterProfiler','AnnotationDbi','DOSE','enrichplot','GO.db','GOSemSim','qvalue', \
        'ReactomePA','graphite','reactome.db','HDO.db','BiocParallel','fgsea','AnnotationHub', \
        'org.Hs.eg.db','ensembldb','EnsDb.Hsapiens.v75','sva','genefilter','BiocFileCache', \
        'biomaRt','rrvgo'), ask = FALSE)"

# Additional GitHub R packages
RUN R -e "install.packages(c('locuszoomr', 'ggtree'), dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    R -e "devtools::install_github('boxiangliu/locuscomparer')"

# MR-MEGA install
RUN mkdir MR-MEGA && cd MR-MEGA && \
    wget https://tools.gi.ut.ee/tools/MR-MEGA_v0.2.zip && \
    unzip MR-MEGA_v0.2.zip && \
    make && \
    wget https://tools.gi.ut.ee/tools/fixP.r && \
    wget https://tools.gi.ut.ee/tools/manh.r && \
    wget https://tools.gi.ut.ee/tools/qq.r && \
    sed -i 's/read.table/fread/g' fixP.r && \
    sed -i 's/write.table/fwrite/g' fixP.r && \
    sed -i 's/data<-/require(data.table)\ndata<-/g' fixP.r && \
    chmod a+x /MR-MEGA/

# SHAPEIT2 install
RUN mkdir /shapeit2 && cd /shapeit2 && \
    wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.17.linux.tar.gz && \
    tar -zxvf shapeit.v2.r904.glibcv2.17.linux.tar.gz && \
    chmod a+x /shapeit2/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit
ENV SHAPEIT_HOME="/shapeit2/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin"

# RFMix install
RUN mkdir /rfmix && cd /rfmix && \
    wget https://github.com/slowkoni/rfmix/archive/refs/heads/master.zip && \
    unzip master.zip && cd rfmix-master && \
    autoreconf --force --install && \
    ./configure && make && \
    chmod a+x rfmix
ENV RFMIX_HOME="/rfmix/rfmix-master"

# Install Miniconda and Tractor environment
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda init && \
    wget https://raw.githubusercontent.com/Atkinson-Lab/Tractor/master/conda_py3_tractor.yml -O /opt/tractor_env.yml && \
    /opt/conda/bin/conda env create -f /opt/tractor_env.yml
ENV CONDA_HOME="/opt/conda/envs/py3_tractor"

# IMPUTE2 install
RUN mkdir /impute2 && cd /impute2 && \
    wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz && \
    tar -zxvf impute_v2.3.2_x86_64_static.tgz && \
    chmod -R a+x /impute2/impute_v2.3.2_x86_64_static
ENV IMPUTE2_HOME="/impute2/impute_v2.3.2_x86_64_static"

# Final PATH setup
ENV PATH="${IMPUTE2_HOME}:${CONDA_HOME}/bin:${RFMIX_HOME}:${SHAPEIT_HOME}:${R_HOME}:${PYTHON_VENV}:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

# Set default working directory
WORKDIR /app

# Clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

CMD ["bash"]