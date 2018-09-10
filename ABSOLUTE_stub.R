library(devtools)
# Instead of command installing the absolute library, Nick loads it in using devtools, I prefer this approach
# though it should be better documented


soft.dir = "/cga/scarter/rklein/Workflows/"
# load_all( paste(soft.dir,"AllelicCapseg/",sep=""), export_all=FALSE )

pkg_loc = paste(soft.dir,"absolutev1.4/library/",sep="")
load_all(pkg_loc, export_all=TRUE )


#library(ABSOLUTE)

RunAbsolute( seg.dat.fn, primary.disease, platform, sample.name, results.dir, copy_num_type,
             genome_build, gender, min.ploidy, max.ploidy,
             max.as.seg.count, max.non.clonal, max.neg.genome,
             maf.fn, indel.maf.fn, min.mut.af,
             output.fn.base, min_probes, max_sd, sigma.h, SSNV_skew,
             filter_segs, force.alpha, force.tau, allelic_capseg_rds, N_threads, verbose )
