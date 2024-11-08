#!/bin/bash

source /carterlab/phahnel/miniforge3/bin/activate /carterlab/phahnel/miniforge3/envs/somix   # LOCAL / KRAKEN
#source /cga/scarter/phahnel/miniconda3/bin/activate /cga/scarter/phahnel/miniconda3/envs/somix   # BROAD CLUSTER

#echo $PATH

# Using anaconda, we may not have access to the /tmp directory
mkdir -p tmp
export TMPDIR=`pwd`/tmp


##### LOCAL & KRAKEN
echo "R at location $(which R)"
R --version
Rscript /carterlab/phahnel/tools/absolute/library/scripts/run_absolute.R "$@" && echo "SUCCESS"


##### BROAD CLUSTER
#. /broad/software/scripts/useuse
#use R-4.1
#use .julia-1.3.1
#export R_LD_LIBRARY_PATH=/broad/software/free/Linux/redhat_7_x86_64/pkgs/julia/julia_1.3.1/lib/julia
#Rscript /carterlab/phahnel/tools/absolute/library/scripts/run_absolute.R "$@" && echo "SUCCESS"


##### DFCI DS KRAKEN CLUSTER
#module add gcc/9.2.0
#module add julia/1.9.0
#module add python/3.9.1
#module add R/4.1
#Rscript /carterlab/phahnel/tools/absolute/library/scripts/run_absolute.R "$@" && echo "SUCCESS"


conda deactivate
