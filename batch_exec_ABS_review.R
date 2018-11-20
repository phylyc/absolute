#!/usr/bin/Rscript
options(error=recover)
.libPaths(c('/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.1.1-bioconductor-3.0/bin/R'))
# reuse .r-3.1.1-bioconductor-3.0


# identify yourself with your initials
analyst.id = "RHK"

# change these to match the project, pair set cohort, etc. 
#obj.name = "icgc_tcga_pdac"
# obj.name = "pdac_targeted_panel_6_26_17"
obj.name = "SDS_batch1"
#obj.name = "full_pdac_8_23_17"
# obj.name = "full_pdac_plus_germline"
#obj.name = "pdac_targeted_panel_4_3_17"
#obj.name = "pdac_targeted_panel_02_27_17"
#obj.name = "pdac_targeted_panel_02_27_17_failed"
# obj.name = "pdac_targeted_panel_01_25_17"
#obj.name = "pdac_targeted_panel_11_16_16"
#obj.name = "full_pdac_11_16_16_no_organoid"
#obj.name = "tcga_pdac_v2"
#obj.name = "full_pdac_10_15_16_no_organoid"
#obj.name = "full_pdac_10_7_16"
#obj.name = "PR_CCPM_PANCREAS_WOLPIN_Capture_All_Pairs_07_25"
#obj.name = "CanSeq_CRSP_Panc_Pairs_07_25"

# software location 
soft.dir = "/cga/scarter/rklein/Workflows/"
#soft.dir = "/xchip/scarter/Software/"


# maf_col = "mutect_maf_file"
# maf_col = "maf_file_oxoG3_capture"
#maf_col = "union_maf_file_forcecalled"
# maf_col = "maf_file_capture_realign"
maf_col = "combined_snvs_indels"


# TRUE, then use a calls file that's been created already and solutions hand picked by user. otherwise, use FALSE
# PP_CALLS_FN = "/xchip/scarter/ncamarda/projects/paad_ccpm/full_pdac_7_13_17.NDC.ABSOLUTE.table.forcecall.txt"
# PP_CALLS_FN = "/cga/scarter/ncamarda/test_snps_msacs/ABSOLUTE_results/test_snps_msacs/man_review.txt"
#PP_CALLS_FN = "/xchip/scarter/ncamarda/projects/paad_ccpm/full_pdac_8_23_17.NDC.ABSOLUTE.table.forcecall.txt"
# PP_CALLS_FN= "/xchip/scarter/ncamarda/projects/paad_ccpm/full_pdac_4_10_17.NDC.ABSOLUTE.table.forcecalled_2.txt"
# PP_CALLS_FN = "/xchip/scarter/ncamarda/projects/paad_ccpm/ABSOLUTE_results/full_pdac_4_10_17/reviewed/full_pdac_4_10_17.NDC.ABSOLUTE.table.txt"  
#PP_CALLS_FN = "/xchip/scarter/ncamarda/projects/paad_ccpm/tcga_pp_calls_2.txt" # ABSOLUTE_results/tcga_pdac_v2/reviewed/tcga_pdac_v2.unmatched.NDC.ABSOLUTE.table.txt"
#PP_CALLS_FN = "/xchip/scarter/ncamarda/projects/paad_ccpm/ABSOLUTE_results/full_pdac_11_16_16_no_organoid/full_pdac_11_16_16_man_review.txt"
#PP_CALLS_FN = "/xchip/scarter/ncamarda/projects/paad_ccpm/ABSOLUTE_results/pdac_targeted_panel_01_25_17/pdac_targeted_panel_01_25_17.PP-calls_tab.txt"
# PP_CALLS_FN = "/xchip/scarter/ncamarda/projects/paad_ccpm/pdac_targeted_panel_02_27_17.passed.table.txt"
# PP_CALLS_FN="/xchip/scarter/ncamarda/projects/paad_ccpm/ABSOLUTE_results/pdac_targeted_panel_02_27_17_failed/force_call.txt"
PP_CALLS_FN = "/cga/scarter/rklein/Projects/SDS/workflows/abs_phy/ABSOLUTE_results/SDS_batch1/man_review.txt"
print("PP CALLS FN:")
print(PP_CALLS_FN)

summarize = T
exclude_called=F

# version control 
ACS.vers = "AllelicCapseg"
#ACS.vers = "ACS_v1.1"
# ACS.vers = "ACS_v1.1.1"
#ABS.vers = "ABSOLUTEv1.3"
ABS.vers = "absolutev1.4"

################################

## only need doMC if you want to use multicore execution
library(doMC)
library(devtools)
threads <- detectCores()
if (threads > 2){ 
  t = floor(threads * (5/8))
} else { t = threads }
registerDoMC(cores = t)
print(paste0("Number of threads: ", t))

# To run this script, you must have ACS_v1.1/library/ and ABSOLUTEv1.4/library copied to your working directory

init.dir = paste(getwd(),"/",sep="")
setwd(init.dir)
ACS.dir = paste(soft.dir,ACS.vers,sep="")
ABS.dir = paste(soft.dir,ABS.vers,sep="")

# load_all(ACS.dir, export_all=FALSE )
load_all(paste(ABS.dir, "/library/",sep=""),export_all=TRUE)


R_STUB_FN = paste(ABS.dir, "/ABSOLUTE_stub.R", sep="")
### setup arguments to RunAbsolute
ABSOLUTE_argv = list()
ABSOLUTE_argv$primary.disease = NA
ABSOLUTE_argv$platform = "Illumina_WES"
ABSOLUTE_argv$copy_num_type = "allelic"
ABSOLUTE_argv$genome_build = "hg19"  ## Note: no other genomes are supported in this release
ABSOLUTE_argv$N_threads = threads
ABSOLUTE_argv$sigma.h = 0.0

# Default values - uncomment to change these
#   ABSOLUTE_argv$min.ploidy = 1.1
#   ABSOLUTE_argv$max.ploidy = 6.0

# not passed to RunAbsolute - only used in stub, and only in dev mode
#ABSOLUTE_argv$CGA_DIR = CGA_DIR
###
SIF_FN = "/cga/scarter/rklein/Projects/SDS/workflows/abs_phy/SIF.tsv"
#SIF_FN = paste(init.dir,sprintf("exfiss/ex_%s.txt",obj.name),sep="")
# SIF_FN = paste(init.dir,"exfiss/ex_full_pdac_7_17_17.txt",sep="")
# SIF_FN = paste(init.dir,"exfiss/ex_pdac_targeted_panel_02_27_17.txt",sep="")
print("SIF FN:")
print(SIF_FN)

#SIF_FN = paste(sprintf("/xchip/scarter/ncamarda/projects/paad_ccpm/wolpin_targeted_panel/final_data_11_16_16/sif.txt",obj.name),sep="")
SIF = read.delim(SIF_FN, row.names=1, check.names=FALSE, stringsAsFactors=FALSE )
results.dir = paste(init.dir,"ABSOLUTE_results/",obj.name,"/",sep="")
sink( file=paste(results.dir, obj.name, "_", Sys.Date(), ".farmer.R.out", sep=""), split=TRUE, append=TRUE )

if (is.na(PP_CALLS_FN)){
  force_call=F
  PP_CALLS_FN = paste(results.dir,"man_review.txt",sep="")
} else { 
  force_call=T 
}

## Step 1: Process the SIF and format args for ABSOLUTE on each sample
##         'calls_FN' should be NA the 1st time this is run
## 		FORCE_CALL = FALSE on first run
var_bsub_argv  = firehose_CAPSEG_SIF( SIF, PP_CALLS_FN=PP_CALLS_FN, FORCE_CALL=force_call, EXCLUDE_CALLED=exclude_called, EXCLUDE_PASSED=FALSE, MAF_COL=maf_col )



## Step 4: Use your 'man_review.txt' to extract the preferred solutions for each sample.
## This will cause the directory:  ./ABSOLUTE_resutls/obj.name/reviewed to be created. This directory contains the ABSOLUTE SCNA and SSVN results for the selected solution of each sample.  Human readable tsv files are in the SEG_MAF/ subdir.   The samples/ subdir contains binary result files that can be input into Phylogic.   In addition, the file: reviewed/obj.name.analyst.id.ABSOLUTE.table.txt will contain a .tsv table of the ABSOLUTE selected solution for each sample.
apply_review_and_extract( pp.review.fn=PP_CALLS_FN, obj.name=obj.name, analyst.id=analyst.id, copy_num_type="allelic"  )

