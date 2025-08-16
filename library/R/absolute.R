## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


.onLoad = function(...) 
{
#   packageStartupMessage("ABSOLUTE v1.2 [BETA] succesfully loaded.")
#    packageStartupMessage("ABSOLUTE v1.4 [ALPHA] succesfully loaded.")
   packageStartupMessage("ABSOLUTE v1.6 [ALPHA] succesfully loaded.")
   options(error = recover)
}

 compute_cached = function( results.dir, file.base, file.ext, fn, verbose, ... )
  {
     result_FN = file.path(results.dir, paste(file.base, "_", file.ext, ".Rds", sep = ""))
     if( file.exists(result_FN))
     {
        if( verbose ) { print( paste( "loading cached ", file.ext, " result", sep="")) }
        obj = try(readRDS( result_FN))

        if( class(obj) == "try-error" ) { print( "load failed!!"); cached=FALSE }
        else { cached = TRUE }
     }
     else{ cached = FALSE }
     if( !cached )
     {
        if( verbose ) { print( paste( "Computing ", file.ext, " result", sep="")) }
        obj = fn(..., verbose=verbose)
        saveRDS( obj, file=result_FN )
     }
     return( obj )
  }


RunAbsolute = function(seg.dat.fn, primary.disease, platform, sample.name, results.dir, copy_num_type, genome_build, gender=NA, min.ploidy=1, max.ploidy=8, max.as.seg.count=1500, max.non.clonal=0.8, max.neg.genome=0.005, maf.fn = NULL, indel.maf.fn = NULL, min.mut.af = NULL, output.fn.base=NULL, min_probes=10, max_sd=100, sigma.h=0.01, SSNV_skew=1, b.res=0.1, d.res=0.01, filter_segs=TRUE, force.alpha=NA, force.tau=NA, allelic_capseg_rds=NA, apply_karyotype_model=FALSE, N_threads=1, verbose = FALSE)
{  
  print( paste("Registering ", N_threads, " threads.", sep=""))
  registerDoMC(N_threads)
  
  genome_build = match.arg(genome_build, c("mm9", "hg18", "hg19", "hg38") )
  
 ##	   3) genome_HSCR_seg_plot.R is currently fixed to hg18 data (in genome.R)

  # if( genome_build == "hg18") { data(hg18_ChrArmsDat, package="ABSOLUTE") }
  # if( genome_build == "hg19") { data(hg19_ChrArmsDat, package="ABSOLUTE") }
  # if( genome_build == "hg38") { data(hg38_ChrArmsDat, package="ABSOLUTE") }
  # if( genome_build == "mm9" ) { data(mm9_ChrArmsDat, package="ABSOLUTE") }

  if ( genome_build == "hg18") { chr.arms.dat.file = file.path(pkg_dir, "data", "hg18_ChrArmsDat.RData") }
  else if ( genome_build == "hg19") { chr.arms.dat.file = file.path(pkg_dir, "data", "hg19_ChrArmsDat.RData") }
  else if ( genome_build == "hg38") { chr.arms.dat.file = file.path(pkg_dir, "data", "hg38_ChrArmsDat.RData") }
  else if ( genome_build == "mm9" ) { chr.arms.dat.file = file.path(pkg_dir, "data", "mm9_ChrArmsDat.RData") }
  else {}

  print(paste("Loading", chr.arms.dat.file))
  load(chr.arms.dat.file)

  platform = match.arg(platform, c("SNP_6.0", "Illumina_WES"))
  if (platform == "SNP_6.0") {
    filter_segs = TRUE
  } else if (platform == "Illumina_WES") {
    filter_segs = TRUE
  } else {
    stop("Unsupported platform: ", platform)
  }
  
  if (copy_num_type == "total") {
    set_total_funcs()
  } else if (copy_num_type == "allelic") {
    set_allelic_funcs()
  } else {
    stop("Unsupported copy number type: ", copy_num_type)
  }

  if (copy_num_type == "total") {
      MakeSegObj <<- total_make_seg_obj
    } else if (copy_num_type == "allelic") {
      MakeSegObj <<- AllelicMakeSegObj
  }
  

## Note - soon we will switch to ASCII HAPSEG output for 6.0, then the HAPSEG functions below will be deprecated and replaced by the Allelic versions   The code block below can then be replaced by platform_funcs.R
  if( platform == "SNP_6.0" ) 
  {
## extract segtab form .RData binary
    segtab = extract_HAPSEG_segtab(seg.dat.fn, verbose=verbose) 
    MakeSegObj <<- AllelicMakeSegObj

#    MakeSegObj <<- HAPSEGMakeSegObj
  }
  else
  {
    if( is.na(SSNV_skew) && !is.na(allelic_capseg_rds) && file.exists(allelic_capseg_rds) )
    {
       print("Overriding seg.dat.fn with allelic_capseg_rds")
       
#       ACS_result = extract_ACS_result(allelic_capseg_rds)
       ACS_result = readRDS(allelic_capseg_rds)
       print("Overriding SSNV_skew with ACS value")
       SSNV_skew = ACS_result[["capture.em.fit"]][["Theta"]][["f_skew"]]

       seg.dat = import_AllelicCapseg_data( ACS_result, gender, min_probes=min_probes, max_sd=max_sd, filter_segs=filter_segs, verbose=verbose)
    }
    else
    {
      if( !is.na(seg.dat.fn) && !file.exists(seg.dat.fn)) { stop("seg.dat.fn does not exist") }

      segtab = read.delim( seg.dat.fn, row.names=NULL, stringsAsFactors=FALSE, check.names=FALSE)

      if (copy_num_type == "allelic") {
        nix = is.na(segtab[,"f"])
        print( paste( "Removing ", sum(nix), " of ", length(nix), " segs with NA f", sep="") )
        if (all(nix)) {
          print("No segments left. Aborting.")
          return(TRUE)
        }
        segtab = segtab[!nix,]
      }

      seg.dat = MakeSegObj(segtab, gender, min_probes=min_probes, max_sd=max_sd,
                           filter_segs=filter_segs, verbose=verbose)
    }
  }

  ##  set up SCNA and SSNV model parameters
  SCNA.argv = list( copy_num_type, min.ploidy, max.ploidy, sigma.h )
  names(SCNA.argv) = c( "copy_num_type", "min.ploidy", "max.ploidy", "sigma.h" )
  SCNA_model = SCNA_model_setup( SCNA.argv, verbose )
  
  tmp.dir = file.path(results.dir, "tmp")
  dir.create(tmp.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results.dir, recursive=TRUE, showWarnings=FALSE)
  file.base = paste(output.fn.base, copy_num_type, "ABSOLUTE", sep = ".")

  seg.dat[["primary.disease"]] = primary.disease
  seg.dat[["group"]] = DetermineGroup(primary.disease)
  seg.dat[["platform"]] = platform
  seg.dat[["copy_num_type"]] = copy_num_type
  seg.dat[["sample.name"]] = as.character(sample.name)
  if (is.null(seg.dat$array.name)) {
    seg.dat$array.name = seg.dat$sample.name
  }
  seg.dat[["maf.fn"]] = maf.fn
  seg.dat[["indel.maf.fn"]] = indel.maf.fn
  
  ## either allelic or total CR, to be modeled.
  seg.dat[["obs.scna"]] = ExtractSampleObs(seg.dat)
  seg.dat[["obs.total.scna"]] = ExtractTotalObs(seg.dat)

  SCNA_model[["N_probes"]] = seg.dat[["obs.scna"]][["N_probes"]]

  ## check for QC failure modes
  if (verbose) {
    print(paste("Expected copy-ratio = ", round( seg.dat[["obs.scna"]][["e.cr"]], 5), sep=""))
  }

  mode.res = list(mode.flag = NA)
  
  if ( length(seg.dat[["obs.scna"]][["W"]] ) > max.as.seg.count) {
    mode.res[["mode.flag"]] = "OVERSEG"
  }
  if ((seg.dat[["obs.scna"]][["e.cr"]] < 0.5) || (seg.dat[["obs.scna"]][["e.cr"]] > 1.5)) {
    mode.res[["mode.flag"]] = "E_CR_SCALE"
  }

  # ## check for QC failure modes
  # if (copy_num_type == "allelic") {
  #   if ((seg.dat[["obs.scna"]][["e.cr"]] < 0.5) || (seg.dat[["obs.scna"]][["e.cr"]] > 1.5)) {
  #     mode.res[["mode.flag"]] = "E_CR_SCALE"
  #   }
  # }
  # if (copy_num_type == "total") {
  #   if ((seg.dat[["obs.scna"]][["e.cr"]] < 1.5) || (seg.dat[["obs.scna"]][["e.cr"]] > 2.5)) {
  #     mode.res[["mode.flag"]] = "E_CR_SCALE"
  #   }
  # }
  
  if (is.na(mode.res[["mode.flag"]])) {
    ## check for MAF describing somatic mutations
    maf = NULL
    if ((!is.na(maf.fn)) && (file.exists(maf.fn))) {
      maf = read.delim(maf.fn, row.names = NULL, stringsAsFactors = FALSE, 
                        check.names = FALSE, na.strings = c("NA", "---"),
                        blank.lines.skip=TRUE, comment.char="#")
    } else {
       print(paste("MAF file: ", maf.fn, " not found.", sep = ""))
    }
    
    indel.maf = NULL
    if (!is.na(indel.maf.fn) && file.exists(indel.maf.fn)) {
      indel.maf = read.delim(indel.maf.fn, row.names = NULL, stringsAsFactors = FALSE, 
                        check.names = FALSE, na.strings = c("NA", "---"),
                        blank.lines.skip=TRUE, comment.char="#")
    } else {
       print(paste("Indel MAF file: ", indel.maf.fn, " not found.", sep = ""))
    }

## find initial purity/ploidy solutions for debugger

#    data(ChrArmsDat, package = "ABSOLUTE")
    if ((!is.null(maf)) && (nrow(maf) > 0)) 
    {
      if (is.na(SSNV_skew)){
	      message("Applying default f_skew value of 0.95")
        SSNV_skew = 0.95
      }
      SSNV_model = init_SSNV_model( SCNA_model[["kQ"]], SSNV_skew, nrow(maf) )
#      mut.cn.dat = classic_CreateMutCnDat(maf, indel.maf, seg.dat, min.mut.af, verbose=verbose)
      mut.cn.dat = classic_CreateMutCnDat(maf, indel.maf, seg.dat, verbose=verbose)
    } else {
      mut.cn.dat = NA
      SSNV_model = NA
    }

  ## Caching for mode.tab
    mode.tab = compute_cached( results.dir, file.base, "mode.tab", ProvisionalModeSweep, verbose, 
               seg.dat, SCNA_model, mut.cn.dat, SSNV_model, force.alpha, force.tau, b.res, d.res, chr.arms.dat )

# For debugging - only process 1st mode:
#   n.modes=1  ## for debugging
#   mode.tab = mode.tab[c(1:n.modes),,drop=FALSE]


## Caching for mode.res
    mode.res = compute_cached( results.dir, file.base, "mode.res", fit_modes_SCNA_models, verbose, 
                               seg.dat, mode.tab, SCNA_model, mut.cn.dat, chr.arms.dat )
  }

  if (is.na(mode.res[["mode.flag"]])) 
  {
    bad.ix = GenomeHetFilter(seg.dat[["obs.scna"]], mode.res, max.non.clonal,
                              max.neg.genome, SCNA_model[["kQ"]], verbose=verbose)
    if (sum(bad.ix) == nrow(mode.res[["mode.tab"]])) {
      mode.res = list(mode.flag="ALPHA_TAU_DOM")
    } else {
      mode.res = ReorderModeRes(mode.res, !bad.ix)
    }
  }

  if (is.na(mode.res[["mode.flag"]]))
  {
    ## 1 - apply karyotype model
    ## Kar model only defined for human cancers
    if (genome_build %in% c("hg18", "hg19", "hg38"))
    {
       # data(ChrArmPriorDb, package="ABSOLUTE")
      load(file.path(pkg_dir, "data", "ChrArmPriorDb.RData"))

      if (seg.dat[["group"]] %in% names(train.obj)) {
        model.id = seg.dat[["group"]]
      } else {
        model.id = "Primary"
        print(paste("Disease type", seg.dat[["group"]], "not in ChrArmPriorDb.RData, set to default model:", model.id))
      }
      mode.res = ApplyKaryotypeModel(mode.res, model.id, train.obj, apply_karyotype_model=apply_karyotype_model)
    } else {
      seg.dat[["group"]] = ""
    }

    ## 2 - apply mutation model
    if ((!is.null(maf)) && (nrow(maf) > 0) && (copy_num_type == "allelic"))
    {
      seg.dat[["mut.cn.dat"]] = mut.cn.dat

## Caching for updated mode.res with SSNV results
      mode.res = compute_cached( results.dir, file.base, "SSNV.mode.res", ApplySSNVModel, verbose, 
                                 mode.res, mut.cn.dat, SSNV_model )

      # bad.ix = ClonalSSNVFilter(mode.res, mut.cn.dat)
      # if (sum(bad.ix) == nrow(mode.res[["mode.tab"]])) {
      #   mode.res = list(mode.flag="ALPHA_TAU_DOM")
      # } else {
      #   mode.res = ReorderModeRes(mode.res, !bad.ix)
      # }
    }

    # bad.ix = RealAlphaFilter(mode.res)
    # if (sum(bad.ix) == nrow(mode.res[["mode.tab"]])) {
    #   mode.res = list(mode.flag="ALPHA_TAU_DOM")
    # } else {
    #   mode.res = ReorderModeRes(mode.res, !bad.ix)
    # }

    mode.res = WeighSampleModes(mode.res)
    mode.res[["call.status"]] = GetCallStatus(mode.res, seg.dat[["obs.scna"]][["W"]])
  }
  
  seg.dat[["mode.res"]] = mode.res
  rm(mode.res); gc()   ## try to save some mem

#  if (is.null(output.fn.base)) {
#    output.fn.base = ifelse(is.null(seg.dat$array.name), sample.name, seg.dat$array.name)
#  }
    
  save(seg.dat, file = file.path(results.dir, paste(file.base, "RData", sep = ".")))
 
## plot result
  if (is.na(seg.dat[["mode.res"]][["mode.flag"]])) {
    sample.pdf.fn = file.path(results.dir, paste(file.base, "plot.pdf", sep = "_"))
    if( verbose ) { print("Making result plot") }

    AbsoluteResultPlot(sample.pdf.fn, seg.dat, chr.arms.dat)
  } else {
    if (verbose) {
      print("Mode flag is NA, not generating plots. Sample has failed ABSOLUTE")
    }    
  }
  
  return(TRUE)
}


