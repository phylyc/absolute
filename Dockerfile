FROM rocker/r-ver:4.3.1

WORKDIR /build
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        zlib1g-dev && \
    Rscript -e ' \
        install.packages(c("optparse", "data.table", "matrixStats", "reshape2", "doMC", "BiocManager")); \
        BiocManager::install("GenomicRanges", update = FALSE, ask = FALSE); \
        for (pkg in c("optparse", "data.table", "matrixStats", "reshape2", "doMC", "GenomicRanges")) { \
            if (!library(pkg, character.only = TRUE, logical.return = TRUE)) quit(status = 10) \
        }' && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* /tmp/* /var/tmp/*

RUN mkdir -p /opt/absolute/library /work
WORKDIR /work
COPY library /opt/absolute/library
CMD ["/bin/bash"]