#!/usr/bin/Rscript
options(error=recover)
.libPaths(c('/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.1.1-bioconductor-3.0/bin/R'))
# unuse R-3.3
# reuse .r-3.1.1-bioconductor-3.0
#' @NOTE:  if force calling, you need to first run review then clear results/ cache!!

#' @concepts 
# 1st - use call_fiss.py or fiss directly to create SIF fiss.txt file in fiss/ directory.
# 2nd - use copy_FH_annotation.R to create the ex_fiss.txt file in the exfiss/ directory.
# 3rd a) - use THIS file, with dry_run = FALSE to submit jobs to uger cluster
# 3rd b1)- use THIS file, batch_exec_ABS.R, with dry_run = TRUE, to create R files
# 3rd b2) - use run_local_ABS.R to run each R file on the local machine.
# 4th - man_review.txt file created using make_override_column.R -> edit the manual_review.txt with selected solutions
# 5th - use batch_exec_ABS_review.R to generate final output

## change these to match the project, pair set cohort, etc.
# obj.name ="ex_extra_samples_pdac_7_15_17"
obj.name = "SDS_batch1"

dry_run = FALSE# if FALSE, submitting jobs to UGER cluster
queue = 'broad'

# software location
#acs.soft.dir = "/xchip/scarter/ncamarda/software/"
soft.dir = "/cga/scarter/rklein/Workflows/"

# version control
ACS.vers = "AllelicCapseg"
#ACS.vers = "ACS_v1.1"
# ACS.vers = "ACS_v1.1.1"
#ABS.vers = "ABSOLUTEv1.3"
ABS.vers = "absolutev1.4"

PP_CALLS_FN = NA
#PP_CALLS_FN = "/xchip/scarter/ncamarda/projects/paad_ccpm/full_pdac_8_23_17.NDC.ABSOLUTE.table.forcecall.txt"
# PP_CALLS_FN = "/xchip/scarter/ncamarda/projects/paad_ccpm/ABSOLUTE_results/pdac_targeted_panel_6_26_17/pdac_targeted_panel_6_26_17.NDC.ABSOLUTE.table.forcecalled.txt"
print("PP Calls FN:")
print(PP_CALLS_FN)

# df <- tibble(alleliccapseg_tsv = dir("/cga/scarter/ncamarda/test_snps_msacs/data/tcr_seg_only/output/MA01-098-8/MSACS/samples", full.names = T, pattern = ".abs.v1.4.seg")) %>%
#   mutate(sid = gsub(".abs.v1.4.seg","", basename(alleliccapseg_tsv)),
#          allelic_capseg_rds = NA,
#          AllelicCapseg_skew = NA)
# new_df <- inner_join(df, read_tsv("/cga/scarter/ncamarda/sds/merge_snvs_and_indels_workflow/run_04.30.18/sample_mafs/cmb_muts_sif.txt") %>% dplyr::select(sid = Tumor_Sample_Barcode,combined_snvs_indels))
# rnames <- new_df$sid
# 
# nn_df <- new_df %>% dplyr::select(-sid) %>% as.data.frame()
# rownames(nn_df) <- rnames
# write.table(nn_df, file = "/cga/scarter/ncamarda/test_snps_msacs/MA01_abs_sif.txt", quote = F, sep = "\t", row.names = T, col.names = T)

# rnames allelic_capseg_rds (NA)      AllelicCapseg_skew (NA)    alleliccapseg_tsv (abs.v1.4.filt.seg)   combined_snvs_indels (maf)
summarize = F
exclude_called = F

#maf_col = NA
# maf_col = "mutect_maf_file"
#maf_col ="maf_file_oxoG3_capture"
maf_col = "combined_snvs_indels"
############################

library(doMC)
library(devtools)
# greater than 1 thread is unstable with qsub
threads = 1
source("/cga/scarter/rklein/Projects/SDS/workflows/abs_phy/make_override_column.R")

# To run this script, you must have ACS_v1.1/library/ and ABSOLUTEv1.4/library copied to software dir
init.dir = paste0(getwd(), "/")
ACS.dir = paste(soft.dir,ACS.vers,sep="")
ABS.dir = paste(soft.dir,ABS.vers,sep="")

# load_all(ACS.dir, export_all=FALSE )
load_all(paste(ABS.dir,"/library/",sep=""),export_all=TRUE)

SIF_FN = paste(init.dir,"SIF.tsv",sep="")
print("SIF FN:")
print(SIF_FN)

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
#SIF_FN = paste(init.dir,"exfiss/ex_", obj.name, ".txt",sep="") 

#SIF_FN = paste(init.dir,"exfiss/ex_","tcga_pdac_v2", ".txt",sep="")
SIF = read.delim(SIF_FN, row.names=1, check.names=FALSE, stringsAsFactors=FALSE )
#results.dir = paste(init.dir,"ABSOLUTE_results/",obj.name,"/",sep="")
# sink( file=paste(init.dir, obj.name, "_", Sys.Date(), ".farmer.R.out", sep=""), split=TRUE, append=TRUE )

## Step 1: Process the SIF and format args for ABSOLUTE on each sample
##         'calls_FN' should be NA the 1st time this is run
##              FORCE_CALL = FALSE if you don't have a calls file. if you do and want to use your solutions, then use TRUE
if (is.na(PP_CALLS_FN)){
  force_call = F
} else { force_call = T }
var_bsub_argv  = firehose_CAPSEG_SIF( SIF, PP_CALLS_FN=PP_CALLS_FN, FORCE_CALL=force_call, EXCLUDE_CALLED=exclude_called, EXCLUDE_PASSED=FALSE, MAF_COL=maf_col )

# write a simple if-else statement that runs batch_exec_ABSOLUTE if the man_review.txt file doesn't exist for the current project and runs_apply_review_and_extract if it does exist

## Step 2: run ABSOLUTE on each sample - create a review table / plot when finished
## resutls will be in ./ABSOLUTE_resutls/obj.name/results/
## The default behavior of this function is to dispatch jobs for each sample and then to wait for them all to finish, at which point the results are gathered and some summaries are output.  You can change this using different arguemnts

# first time through will submit jobs to cluster and try to aggregate - it will inevitably fail because it needs the review process, which can't happen unless you run the dry_run= TRUE case. so run these in order!!

batch_exec_ABSOLUTE( ABSOLUTE_argv, obj.name, var_bsub_argv, R_STUB_FN, queue=queue, dry_run=dry_run, EXE_ENGINE="UGER" )
#batch_exec_ABSOLUTE( ABSOLUTE_argv, obj.name, var_bsub_argv, R_STUB_FN, queue="short", dry_run=(!dry_run), EXE_ENGINE="UGER" )


## Step 3: review .pdf files in the results output and manually select solutions, if neccesary.
##      a. Create a 'man_review.txt' file by copying: ./ABSOLUTE_resutls/obj.name/results/obj.name.PP-calls_tab.txt
##         and inserting a new 1st column named 'override'.  For each sample (row), enter the # of the solution you want to pick or else leave blank to accept the automatic call.  You can also enter special labels such as 'low purity', or 'FAILED' to prevent calls from being made on these samples
## See the tutorial http:// for tips on manually reviewing solutions

## write a script that automatically creates this file
make_override_column(obj.name)


## Step 4: Use your 'man_review.txt' to extract the preferred solutions for each sample.
## This will cause the directory:  ./ABSOLUTE_resutls/obj.name/reviewed to be created. This directory contains the ABSOLUTE SCNA and SSVN results for the selected solution of each sample.  Human readable tsv files are in the SEG_MAF/ subdir.   The samples/ subdir contains binary result files that can be input into Phylogic.   In addition, the file: reviewed/obj.name.analyst.id.ABSOLUTE.table.txt will contain a .tsv table of the ABSOLUTE selected solution for each sample.
# apply_review_and_extract( pp.review.fn=PP_CALLS_FN, obj.name=obj.name, analyst.id="NDC", copy_num_type="allelic"  )
