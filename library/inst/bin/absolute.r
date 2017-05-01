#!/usr/bin/env Rscript

library(ABSOLUTE, quietly = TRUE)
library(optparse, quietly = TRUE)

option_list <- list(make_option(c("--results_dir"), type = "character", default = NULL, 
                               help="results local directory [required]", metavar = "string"),
                   make_option(c("--sample"), type = "character", default = NULL, 
                               help = "sample name [required]", metavar = "string"),
                   make_option(c("--rds"), type = "character", default = NA, 
                               help = "RDS file [required]", metavar = "string"),
                   make_option(c("--maf"), type = "character", default = NA, 
                               help = "somatic SNV MAF file [optional]", metavar = "string"),
                   make_option(c("--indel_maf"), type = "character", default = NA, 
                               help = "somatic indel MAF file [optional]", metavar = "string"),
                   make_option(c("--gender"), type = "character", default = NA, 
                               help = "gender [Male or Female] or biological sex of individual [XY or XX]", metavar = "string"),
                   make_option(c("--alpha"), type = "double", default = NA, 
                               help = "alpha", metavar = "number"),
                   make_option(c("--tau"), type = "double", default = NA, 
                               help = "tau", metavar = "number"),
                   make_option(c("--min_ploidy"), type = "double", default = 1.1, 
                               help = "minimum ploidy [default= %default]", metavar = "number"),
                   make_option(c("--max_ploidy"), type = "double", default = 6, 
                               help = "maximum ploidy [default= %default]", metavar = "number"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$results_dir)) {
  print_help(opt_parser)
  stop("results_dir must be provided", call. = FALSE)
} else if (is.null(opt$sample)) {
  print_help(opt_parser)
  stop("sample must be provided", call. = FALSE)
} else if (is.null(opt$rds)) {
  print_help(opt_parser)
  stop("rds must be provided", call. = FALSE)
}

primary.disease <- NA
platform <- "Illumina_WES"
copy_num_type <- "allelic"
genome_build <- "hg19"
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
SSNV_skew <- 0.95
filter_segs <- TRUE
force.alpha <- opt$alpha
force.tau <- opt$tau
allelic_capseg_rds <- opt$rds
verbose <- TRUE
sample.name <- opt$sample
seg.dat.fn <- NA
gender <- opt$gender

RunAbsolute(seg.dat.fn, primary.disease, platform, sample.name, results.dir, copy_num_type, genome_build, gender,
            min.ploidy, max.ploidy, max.as.seg.count, max.non.clonal, max.neg.genome, maf.fn, indel.maf.fn, min.mut.af,
            output.fn.base, min_probes, max_sd, sigma.h, SSNV_skew, filter_segs, force.alpha, force.tau,
            allelic_capseg_rds, N_threads, verbose)
