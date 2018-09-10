### CHANGE THESE TO YOUR OWN 
R_STUB_FN = "/fullpathto/ABSOLUTE_stub.R"       ## Change to point to your own install directory
SIF_FN = "/yourproject/SIF.txt"                 ## your own Sample info file
obj.name = "initial_ABSv1.3_run"  ## Speficies subdir for ABS results.  Rename this to test different runs on the same project 

## NOTE:It is now possible to specify runtimes and runmem for absolute qsub.
## To do so, alter those values in the batch_exec_ABSOLUTE function below
## Where time is in the format HH:MM:SS and memory is a numeric amount followed
## by a unit -- G is recommended


## Example SIF columns:
#SID	AllelicCapseg_skew	alleliccapseg_tsv	maf_file_SSNV	maf_file_indel
# SID: unique sample id
# AllelicCapseg_skew: real number, output from Allelic Capseg
# alleliccapseg_tsv: tsv file with allelic segmention summary from AllelicCapseg
# maf_file_SSNV: tsv file with one line per SSNV 
# maf_file_indel: tsv file with one line per indel

###

# SLC test:
#SIF_FN = "~/Projects/PAAD_CCPM/exFH_SIF_10.14.2015.txt"
#obj.name = "All_pairs_test_ABSv1.3"
#R_STUB_FN = "~/Projects/ABSv1.3/ABSOLUTE_stub.R"    

###











options(error=recover)

require(doMC)
registerDoMC(16)

library(ABSOLUTE)

#library(devtools)
#pkg_loc = "~/ABS_lib/"
#load_all(pkg_loc, export_all=TRUE )

### setup arguments to RunAbsolute
ABSOLUTE_argv = list()
ABSOLUTE_argv$primary.disease = NA
ABSOLUTE_argv$platform = "Illumina_WES"
ABSOLUTE_argv$copy_num_type = "allelic"
ABSOLUTE_argv$genome_build = "hg19"  ## Note: no other genomes are supported in this release

# Default values - uncomment to change these
#   ABSOLUTE_argv$min.ploidy = 1.1
#   ABSOLUTE_argv$max.ploidy = 6.0

# not passed to RunAbsolute - only used in stub, and only in dev mode
#ABSOLUTE_argv$CGA_DIR = CGA_DIR
###

# Extract data from Firehose:
#fiss annot_get Analysis_CCPM_PANCREAS_WOLPIN pair pset=PR_CCPM_PANCREAS_WOLPIN_Capture_All_Pairs maf_file_capture_master_filter_removed strelka_passed_somatic_indel_maf_file_capture_pair alleliccapseg_tsv AllelicCapseg_skew > FH_SIF_10.14.2015.txt


SIF = read.delim(SIF_FN, row.names=1, check.names=FALSE, stringsAsFactors=FALSE )
sink( file=paste(obj.name, "_", Sys.Date(), ".farmer.R.out", sep=""), split=TRUE, append=TRUE )


#PP_CALLS_FN = "All_pairs_10.14.2015.SLC.ABSOLUTE.table.txt"
PP_CALLS_FN = NA

## Step 1: Process the SIF and format args for ABSOLUTE on each sample
##         'calls_FN' should be NA the 1st time this is run
var_bsub_argv  = firehose_CAPSEG_SIF( SIF, PP_CALLS_FN=PP_CALLS_FN, FORCE_CALL=FALSE, EXCLUDE_CALLED=FALSE, EXCLUDE_PASSED=FALSE )

## Step 2: run ABSOLUTE on each sample - create a review table / plot when finished
## resutls will be in ./ABSOLUTE_resutls/obj.name/results/
## The default behavior of this function is to dispatch jobs for each sample and then to wait for them all to finish, at which point the results are gathered and some summaries are output.  You can change this using different arguemnts
batch_exec_ABSOLUTE( ABSOLUTE_argv, obj.name, var_bsub_argv, R_STUB_FN, queue="short", dry_run=FALSE, EXE_ENGINE="UGER", run_time="6:00:00", run_mem="10G" )

## Step 3: review .pdf files in the results output and manually select solutions, if neccesary.
## 	a. Create a 'man_review.txt' file by copying: ./ABSOLUTE_resutls/obj.name/results/obj.name.PP-calls_tab.txt
##	   and inserting a new 1st column named 'override'.  For each sample (row), enter the # of the solution you want to pick or else leave blank to accept the automatic call.  You can also enter special labels such as 'low purity', or 'FAILED' to prevent calls from being made on these samples 
## See the tutorial http:// for tips on manually reviewing solutions



## Step 4: Use your 'man_review.txt' to extract the preferred solutions for each sample.
## This will cause the directory:  ./ABSOLUTE_resutls/obj.name/reviewed to be created. This directory contains the ABSOLUTE SCNA and SSVN results for the selected solution of each sample.  Human readable tsv files are in the SEG_MAF/ subdir.   The samples/ subdir contains binary result files that can be input into Phylogic.   In addition, the file: reviewed/obj.name.analyst.id.ABSOLUTE.table.txt will contain a .tsv table of the ABSOLUTE selected solution for each sample.
apply_review_and_extract( pp.review.fn=PP_CALLS_FN, obj.name=obj.name, analyst.id="SLC", copy_num_type="allelic"  )



# Optional: you can also use a reviewed calls table to 'FORCE' absolute to call samples and extract data at only pre-specified purity/ploidy solutions.   This is occasionaly useful when samples are very impure or you want to evalute purity/ploidy calls obtained by other means.
# 
#calls_FN = paste( "ABSOLUTE_results/",  obj.name,  "reviewed/",   obj.name, ".", analyst.id, sep="" )
## Edit the calls_FN by modifying the "purity" and "tau"  columns - note the 'tau' column is where you should set the desired ploidy.
#var_bsub_argv  = firehose_CAPSEG_SIF( SIF, calls_FN=calls_FN, FORCE_CALL=TRUE, EXCLUDE_CALLED=FALSE, EXCLUDE_PASSED=FALSE )
#batch_exec_ABSOLUTE( ABSOLUTE_argv, obj.name, var_bsub_argv, R_STUB_FN, queue="hour", dry_run=FALSE, EXE_ENGINE="LSF" )
#apply_review_and_extract( pp.review.fn=calls_FN, obj.name=obj.name, analyst.id="SLC", copy_num_type="allelic"  )


## Alternately, if you already have a calls_FN file, from another program or from a previous version of ABSOLUTE, you can use it to 'liftover' solutions to the current version of ABSOLUTE.   In this mode, ABSOLUTE will search for candidate solutions that very close to the 'purity', 'tau' columns of the calls_FN file.   If it finds a match, it will accept the solution, otherwise, the sample will be uncalled.  After the matching step, data for matched samples will be extracted.
## To run in this mode, use a calls_FN file that matches the format produced in the ./ABSOLUTE_resutls/obj.name/results/obj.name.PP-calls_tab.txt file.   Note that you do NOT have to insert the 'override' column as when you are manually selecting solutions to pick.
#apply_review_and_extract( pp.review.fn=calls_FN, obj.name=obj.name, analyst.id="SLC", copy_num_type="allelic"  )

## Note - you can set up the batch args to only run samples that have not yet been called (i.e. samples that do not appear in the calls_FN file.   This is useful if you want to generate results for review only on a subset of samples, e.g. those that did not match during liftover.
#var_bsub_argv  = firehose_CAPSEG_SIF( SIF, calls_FN=calls_FN, FORCE_CALL=FALSE, EXCLUDE_CALLED=TRUE, EXCLUDE_PASSED=FALSE )
## now run the batch_exec_ABSOLUTE function and apply_review_and_extract functions as above




## step x: gene-level SCNA calling
## Threshold parameters for calling SCNAs
   SCNA_thresholds = get_SCNA_thresholds()

## Experiment with trackData.  needed for gr.match()
# export GIT_HOME = "/cga/meyerson/home/marcin/DB/git"...
   Sys.setenv( 'GIT_HOME'=  "/cga/meyerson/home/marcin/DB/git" )
   ISVA.HOME = paste(Sys.getenv('GIT_HOME'), 'isva', sep = "/");
   source(paste(ISVA.HOME, 'grUtils.R', sep = "/"))
   source(paste(ISVA.HOME, 'trackData.R', sep = "/"))

## Obtain gUtils library from: https://github.com/mskilab/gUtils


   out.dir.base = file.path( "ABSOLUTE_results", obj.name )
   transcript_GRs = get_GENCODE_transcript_GRs()
   gene_SCNA_calls = genotype_transcript_SCNAs_in_called_ABS_files( file.path(out.dir.base, "reviewed", "samples"), transcript_GRs, SCNA_thresholds, analyst_id = "SLC" )
   saveRDS( gene_SCNA_calls, file.path(out.dir.base, "reviewed", paste(obj.name, "_SCNA_genotypes.rds", sep="")) )


