#!/bin/bash

# This is for installation on the DFCI DS kraken cluster.

#conda create -n absolute -y
#conda activate absolute
#source /data/rheinbay_lab/phahnel/<miniforge>/bin/activate /carterlab/phahnel/<miniforge>/envs/absolute

echo $PATH

echo "Installing R ... "
conda install r-base=3.6.3 -y
which R

Rscript -e 'install.packages("optparse", repos="https://cran.r-project.org"); \
  if (!library(optparse, logical.return = TRUE)) quit(status = 10)'
Rscript -e 'install.packages("data.table", repos="https://cran.r-project.org"); \
  if (!library(data.table, logical.return = TRUE)) quit(status = 10)'
Rscript -e 'install.packages("matrixStats", repos="https://cran.r-project.org"); \
  if (!library(matrixStats, logical.return = TRUE)) quit(status = 10)'
Rscript -e 'install.packages("reshape2", repos="https://cran.r-project.org"); \
  if (!library(reshape2, logical.return = TRUE)) quit(status = 10)'
Rscript -e 'install.packages("doMC", repos="https://cran.r-project.org"); \
  if (!library(doMC, logical.return = TRUE)) quit(status = 10)'
Rscript -e 'install.packages("BiocManager", repos="https://cran.r-project.org"); \
  if (!library(BiocManager, logical.return = TRUE)) quit(status = 10)'
Rscript -e 'BiocManager::install("GenomicRanges"); \
  if (!library(GenomicRanges, logical.return = TRUE)) quit(status = 10)'

echo "SUCCESS"
