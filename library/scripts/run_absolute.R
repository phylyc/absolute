#!/usr/bin/env Rscript

options(warn=1)

suppressPackageStartupMessages({
  #library(ABSOLUTE)
  library(optparse)
  library(data.table)
  library(reshape2)
  library(matrixStats)
  library(doMC)
  library(GenomicRanges)
  # library(devtools)
})

option_list <- list(
  make_option("--results_dir", type = "character", default = NULL, help="results local directory [required]", metavar = "string"),
  make_option("--sample", type = "character", default = NULL, help = "sample name [required]", metavar = "string"),
  make_option("--rds", type = "character", default = NA, help = "Path to rds file of segmentation data", metavar = "string"),
  make_option("--seg_dat_fn", type = "character", default = NA, help = "Path to tab delimited segmentation data", metavar = "string"),
  make_option("--maf", type = "character", default = NA, help = "somatic SNV MAF file [optional]", metavar = "string"),
  make_option("--indel_maf", type = "character", default = NA, help = "somatic indel MAF file [optional]", metavar = "string"),
  make_option("--gender", type = "character", default = NA, help = "gender [Male or Female] or biological sex of individual [XY or XX]", metavar = "string"),
  make_option("--alpha", type = "double", default = NA, help = "alpha (purity)", metavar = "number"),
  make_option("--tau", type = "double", default = NA, help = "tau (ploidy)", metavar = "number"),
  make_option("--ssnv_skew", type = "double", default = 0.9883274, help = "skew", metavar = "number"),
  make_option("--min_ploidy", type = "double", default = 1.1, help = "minimum ploidy [default= %default]", metavar = "number"),
  make_option("--max_ploidy", type = "double", default = 6, help = "maximum ploidy [default= %default]", metavar = "number"),
  make_option("--b_res", type = "double", default = 0.025, help="resolution of b grid during provisional mode sweep", metavar = "number"),
  make_option("--delta_res", type = "double", default = 0.01, help="resolution of delta grid during provisional mode sweep", metavar = "number"),
  make_option("--copy_num_type", type = "character", default = "allelic", help = "type: allelic or total [default= %default]", metavar = "string"),
  make_option("--primary_disease", type = "character", default = NA, help = "Disease type of the primary tumor [default= %default]", metavar = "string"),
  make_option("--apply_karyotype_model", type = "logical", action = "store_true", default = FALSE, help = "Apply chromosome arm-level SCNA priors based on disease type [default= %default]", metavar = "string"),
  make_option("--genome_build", type = "character", default = "hg19", help = "build of the genome: hg18, hg19, hg38, mm9 [default= %default]", metavar = "string"),
  make_option("--pkg_dir", type = "character", default = ".", help = "package directory", metavar = "string")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

cat("Input options:\n")
for (name in names(opt)) {
  cat(sprintf("  %s: %s\n", name, toString(opt[[name]])))
}

if (is.null(opt$results_dir)) {
  print_help(opt_parser)
  stop("results_dir must be provided", call. = FALSE)
} else if (is.null(opt$sample)) {
  print_help(opt_parser)
  stop("sample must be provided", call. = FALSE)
} else if (is.null(opt$rds) || is.null(opt$seg_dat_fn)) {
  print_help(opt_parser)
  stop("rds or seg_dat_fn must be provided", call. = FALSE)
}

primary.disease <- opt$primary_disease
apply_karyotype_model <- opt$apply_karyotype_model
platform <- "Illumina_WES"
copy_num_type <- opt$copy_num_type
genome_build <- opt$genome_build
N_threads <- 1
results.dir <- opt$results_dir
min.ploidy <- opt$min_ploidy
max.ploidy <- opt$max_ploidy
max.as.seg.count <- 5000
max.non.clonal <- 0.99
max.neg.genome <- 0.05
maf.fn <- opt$maf
indel.maf.fn <- opt$indel_maf
min.mut.af <- 0
output.fn.base <- opt$sample
min_probes <- 1
max_sd <- 100
sigma.h <- 0.01
SSNV_skew <- opt$ssnv_skew
b.res <- opt$b_res
d.res <- opt$delta_res
filter_segs <- TRUE
force.alpha <- opt$alpha
force.tau <- opt$tau
allelic_capseg_rds <- opt$rds
verbose <- TRUE
sample.name <- opt$sample
seg.dat.fn <- opt$seg_dat_fn
gender <- opt$gender
pkg_dir <- file.path(opt$pkg_dir, "library")

# load_all(file.path(pkg_dir, "library"), export_all = FALSE)

print(paste0("sourcing files in ", file.path(pkg_dir, "R")))
rr = dir(file.path(pkg_dir, "R") , full.names=TRUE, pattern = "*.R" )
for( i in 1:length(rr) ) {
  # print(paste0("sourcing ", rr[i]))
  source(rr[i])
}


RunAbsolute(
  seg.dat.fn, primary.disease, platform, sample.name, results.dir, copy_num_type, genome_build, gender,
  min.ploidy, max.ploidy, max.as.seg.count, max.non.clonal, max.neg.genome, maf.fn, indel.maf.fn, min.mut.af,
  output.fn.base, min_probes, max_sd, sigma.h, SSNV_skew, b.res, d.res, filter_segs, force.alpha, force.tau,
  allelic_capseg_rds, apply_karyotype_model, N_threads, verbose
)
