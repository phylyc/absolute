FROM ubuntu:20.04
WORKDIR /build
ENV DEBIAN_FRONTEND=noninteractive

# Update package list and install R + build dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libcurl4-gnutls-dev \
        libssl-dev \
        libxml2-dev \
        zlib1g-dev \
        r-base=3.6.3-* && \
    Rscript -e ' \
        install.packages(c("optparse", "data.table", "matrixStats", "reshape2", "doMC", "BiocManager"), \
                         repos="https://cran.r-project.org"); \
        BiocManager::install("GenomicRanges"); \
        for (pkg in c("optparse", "data.table", "matrixStats", "reshape2", "doMC", "GenomicRanges")) { \
            if (!library(pkg, character.only = TRUE, logical.return = TRUE)) quit(status = 10) \
        }' && \
    apt-get purge -y build-essential && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir -p /opt/absolute/library /work
WORKDIR /work
COPY library /opt/absolute/library
CMD ["/bin/bash"]