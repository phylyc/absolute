# Specify architecture
#FROM --platform=$BUILDPLATFORM ubuntu:20.04\
FROM ubuntu:20.04
WORKDIR /build
ENV DEBIAN_FRONTEND=noninteractive

# Update package list and install R dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libcurl4-gnutls-dev \
        libssl-dev \
        libxml2-dev \
        zlib1g-dev \
        r-base=3.6.3-*
# R packages
RUN Rscript -e 'install.packages("optparse", repos="https://cran.r-project.org"); \
                if (!library(optparse, logical.return = TRUE)) quit(status = 10)'
RUN Rscript -e 'install.packages("data.table", repos="https://cran.r-project.org"); \
                if (!library(data.table, logical.return = TRUE)) quit(status = 10)'
RUN Rscript -e 'install.packages("matrixStats", repos="https://cran.r-project.org"); \
                if (!library(matrixStats, logical.return = TRUE)) quit(status = 10)'
RUN Rscript -e 'install.packages("reshape2", repos="https://cran.r-project.org"); \
                if (!library(reshape2, logical.return = TRUE)) quit(status = 10)'
RUN Rscript -e 'install.packages("doMC", repos="https://cran.r-project.org"); \
                if (!library(doMC, logical.return = TRUE)) quit(status = 10)'
RUN Rscript -e 'install.packages("BiocManager", repos="https://cran.r-project.org"); \
                if (!library(BiocManager, logical.return = TRUE)) quit(status = 10)'
RUN Rscript -e 'BiocManager::install("GenomicRanges"); \
                if (!library(GenomicRanges, logical.return = TRUE)) quit(status = 10)'
# RUN Rscript -e 'install.packages("devtools", repos="https://cran.r-project.org"); if (!library(devtools, logical.return = TRUE)) quit(status = 10)'
# Reduce size of image
RUN apt-get purge -y build-essential && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /app
RUN mkdir -p /library
COPY library /library
CMD ["/bin/bash"]