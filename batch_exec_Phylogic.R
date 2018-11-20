rm(list=ls())
options(error=recover) 
.libPaths(c("/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.1.1-bioconductor-3.0/lib64/R/library"))
# run with 1 thread if using UGER, known issues with multi-threading
#threads = 1
#registerDoMC(threads)
library(doMC)
library(devtools)
threads = 1
# who are you? 
analyst_id = "RHK"

# name of the pair_set in the project, required always. be sure to change this for each  
obj.name = "SDS_batch1"
# obj.name = "full_pdac_8_23_17"
# obj.name = "full_pdac_11_16_16_no_organoid"
#obj.name = "tcga_pdac_v2"
#obj.name = "full_pdac_10_15_16_no_organoid"
#obj.name = "pdac_tcga_mafs"
#obj.name = "PR_CCPM_PANCREAS_WOLPIN_Capture_All_Pairs_07_25"
#obj.name = "CanSeq_CRSP_Panc_Pairs_07_25"

# these patient and sample files need to be created up front
# NOTE: working on a way to have these be automatically created and updated. 
# Patient_SIF_FN: tab-delimited file with header = patient_id, Histology group, sample_type, Whole genome doubling
# store these files in folder work.dir/phylogic_sample_info/obj.name/
#patient_SIF = "full_pdac_10_7_16_patient_SIF.txt" 
patient_SIF = sprintf("%s_patient_SIF.txt",obj.name)

# Sample_SIF_FN: tab-delimited file with header = sample, pair_id, patient_id, sample_type
#sample_SIF = "full_pdac_10_7_16_Phylogic_SIF.txt" 
sample_SIF = sprintf("%s_Phylogic_SIF.txt",obj.name)

# what kind of samples are these? 
sample_types = c("Primary", "Met")
sample_type_colors = c("dodgerblue", "darkviolet") 
#sample_types = c("Primary", "Organoid")
#sample_type_colors = c("dodgerblue", "darkviolet") 

# select TARGET mutation file
#curated_mutation_FN = "PDAC_TARGET_v9_061716.txt"
# curated_mutation_FN="/xchip/scarter/ncamarda/projects/paad_ccpm/PDAC_TARGET_v11_83117.txt"
# print("TARGET DB:")
# print(curated_mutation_FN)
# software location
soft.dir = "/cga/scarter/rklein/Workflows/"
#soft.dir = "/xchip/scarter/Software/"


# version control
ACS.vers = "AllelicCapseg"
# ACS.vers = "ACS_v1.1"
#ACS.vers = "ACS_v1.1.1"
#ABS.vers = "ABSOLUTEv1.3"
ABS.vers = "absolutev1.4"
phylogic.vers = "phylogic"


dry_run = FALSE
# ============================================================================================


# make sure you have Phylogic library copied into your software directory
init.dir = paste(getwd(),"/",sep="")
setwd(init.dir)
phylogic.path = paste(soft.dir,phylogic.vers,sep="")

# link to exfiss file, created during ABSOLUTE run
# exfiss.dir = paste(init.dir,"exfiss/",sep="")
# FH_SIF_FN = paste(exfiss.dir,"ex_",obj.name,".txt",sep="")
FH_SIF_FN = paste0(init.dir, "SIF.tsv")

patient_SIF_FN = patient_SIF
sample_SIF_FN = sample_SIF

# assume that Absolute and ACS libraries are already copied into the current directory...
phylogic.lib = paste(phylogic.path,"/library/",sep="")
#ACS.lib = paste(soft.dir,ACS.vers,sep="")
ABS.lib = paste(soft.dir,ABS.vers,"/library/",sep="")

#load_all( ACS.lib,export_all=FALSE )
load_all( ABS.lib, export_all = TRUE)
load_all( phylogic.lib, export_all=TRUE )
# load_all("~/../scarter/Phylogic_lib", export_all = TRUE)
load_all( paste(soft.dir,"mixr/library/",sep=""), export_all=TRUE ) 

### setup arguments to Phylogic
Phylogic_argv = list()
#Phylogic_argv$filter_sites_FN = "filter_Brastianos_2015_PoN_sites.tsv"
  ## other args will be set to defaults

# find the ABSOLUTE results dir with phylogic input samples for this project
# id.called.RData file -- absolute seg file 
ABSOLUTE_DIR = sprintf("ABSOLUTE_results/%s/reviewed/samples/",obj.name)
#ABSOLUTE_DIR = "ABSOLUTE_results/full_pdac_10_4_16/reviewed/samples/"

# generate Phylogic SIF FN using obj.name and analyst id
Phylogic_SIF = make_Phylogic_SIF( patient_SIF_FN, sample_SIF_FN, FH_SIF_FN, ABSOLUTE_DIR, analyst_id )
Phylogic_SIF_FN = "test_Phylogic_SIF.tsv"
write.table( Phylogic_SIF, file=Phylogic_SIF_FN, quote=FALSE, sep="\t", row.names=FALSE )

names(sample_type_colors) = sample_types

Phylogic_argv$sample_type_colors = sample_type_colors
Phylogic_argv$model_type="sibling_bottleneck"
#Phylogic_argv$model_type="free_ND" 

sink( file=paste(obj.name, "_", Sys.Date(), ".Phylogic.farmer.R.out", sep=""), split=TRUE, append=TRUE )

# EXE_ENGINE defaults to MULTICORE - untested

Phylogic_mode="ND"
Phylogic_argv[["N_iter"]] = 500
Phylogic_argv[["N_threads"]] = 1
Phylogic_argv[["N_chains"]] = 1
#Phylogic_mode="2D_histogram"


results.main  =paste(init.dir, "Phylogic_", Phylogic_mode, "_results", sep="") 
out.dir= file.path( results.main, obj.name )
results.dir = file.path( out.dir, 'results' )

# Read target.txt file 
# TARGET = read.delim(curated_mutation_FN, check.names=FALSE, stringsAsFactors=FALSE )
data("VanAllen2014_TARGET", package="ABSOLUTE")
TARGET <- TARGET[,-5]
# colnames(TARGET)[colnames(TARGET)=="PKB Therapy Group"] = "Therapy group"

# TCGA_ICGC_ix = TARGET[,"PDAC_TCGA_ICGC"]=="YES"
# TCGA_ICGC = TARGET[TCGA_ICGC_ix,"Gene"]
# 
# #PDAC_Germline_Risk_ix = TARGET["PDAC_Germline_Risk"] == "YES"
# #PDAC_Germline_Risk = TARGET[PDAC_Germline_Risk_ix, "Gene"]
# 
# Clinically_Actionable_ix = TARGET[,"Clinically_Actionable"]=="YES"
# Clinically_Actionable = TARGET[Clinically_Actionable_ix,"Gene"]
# 
# DDR_ix = TARGET[,"DDR"]=="YES"
# DDR = TARGET[DDR_ix,"Gene"]
# 
genelists = list( Gene = TARGET$Gene )

summarize_argv = list(Phylogic_SIF_FN, TARGET, patient_SIF_FN, out.dir, ABSOLUTE_DIR, analyst_id, results.dir, genelists)
names(summarize_argv) = c("Phylogic_SIF_FN", "TARGET", "patient_SIF_FN", "out.dir", "ABSOLUTE_DIR", "analyst_id", "results.dir", "genelists")
EXOME_SIZE=32928527 
summarize_argv$covered_bases = EXOME_SIZE 
summarize_argv$sample_type_colors = sample_type_colors
summarize_argv$amp_gene_list = get_TARGET_amp_genes(TARGET)
summarize_argv$del_gene_list = get_TARGET_del_genes(TARGET)

# set summarize to False first...

batch_exec_Phylogic( Phylogic_argv, Phylogic_SIF_FN, ABSOLUTE_DIR, analyst_id, obj.name, Phylogic_mode=Phylogic_mode, queue="broad", overwrite=FALSE, dry_run=dry_run, wait=TRUE, summarize=FALSE, EXE_ENGINE="UGER", groupname="" )


# after .RData files have been created, summarizes, for all sample ids
# run batch_exec_Phylogic_review.R
# summarize_ND_Phylogic( summarize_argv )
#save(summarize_argv, file=paste0(results.dir,"/summarize_argv.RData") )

