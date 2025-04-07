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
  make_option("--analyst_id", type = "character", default = NULL, help="acronym of analyst who called the solution [required]", metavar = "string"),
  make_option("--solution_num", type = "integer", default = NULL, help = "ordinal number of solution [required]", metavar = "number"),
  make_option("--sample", type = "character", default = NULL, help = "sample name [required]", metavar = "string"),
  make_option("--rdata", type = "character", default = NULL, help = "Path to rdata file of Absolute output [required]", metavar = "string"),
  make_option("--copy_num_type", type = "character", default = "allelic", help = "type: allelic or total [default= %default]", metavar = "string"),
  make_option("--genome_build", type = "character", default = "hg19", help = "build of the genome: hg18, hg19, mm9 [default= %default]", metavar = "string"),
  make_option("--pkg_dir", type = "character", default = ".", help = "package directory", metavar = "string")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$results_dir)) {
  print_help(opt_parser)
  stop("results_dir must be provided", call. = FALSE)
} else if (is.null(opt$sample)) {
  print_help(opt_parser)
  stop("sample must be provided", call. = FALSE)
} else if (is.null(opt$rdata)) {
  print_help(opt_parser)
  stop("rdata must be provided", call. = FALSE)
}

analyst.id <- opt$analyst_id
copy_num_type <- opt$copy_num_type
genome_build <- opt$genome_build
N_threads <- 1
results.dir <- opt$results_dir
output.fn.base <- opt$sample
rdata <- opt$rdata
verbose <- TRUE
sample.name <- opt$sample
solution_num <- opt$solution_num
pkg_dir <- file.path(opt$pkg_dir, "library")

# load_all(file.path(pkg_dir, "library"), export_all = FALSE)

print(paste0("sourcing files in ", file.path(pkg_dir, "R")))
rr = dir(file.path(pkg_dir, "R") , full.names=TRUE, pattern = "*.R" )
for( i in 1:length(rr) ) {
  # print(paste0("sourcing ", rr[i]))
  source(rr[i])
}
registerDoMC(N_threads)


apply_review_and_extract(
  pp.solution.num=solution_num, results.dir=results.dir, rdata=rdata, obj.name=sample.name, analyst.id=analyst.id, copy_num_type=copy_num_type, genome_build=genome_build
)
