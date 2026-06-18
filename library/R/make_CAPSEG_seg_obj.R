## functions for running ABSOLUTE with total CR

total_make_seg_obj = function(segs_tab, gender, filter_segs=FALSE, min_probes=NA, max_sd=NA, verbose=FALSE) {
  # segs_tab = read.delim(dat_fn, row.names=NULL, check.names=FALSE, stringsAsFactors=FALSE)
  seg_dat = list()
  
  # print("Gender not supported yet!!")
  # gender = NA

  X.ix = segs_tab[,"Chromosome"]==23
  segs_tab[X.ix, "Chromosome"] = "X"
  Y.ix = segs_tab[,"Chromosome"]==24
  segs_tab[Y.ix, "Chromosome"] = "Y"


## normalize gender col
  is.M = gender %in% c("M", "male", "Male", "XY")
  is.F = gender %in% c("F", "female", "Female", "XX")

  if(is.M) { gender = "Male" }
  if(is.F) { gender = "Female" }

  if( verbose ) { print(paste("Detected ", gender, " gender.", sep="")) }

  segs_tab = filter_sex_chromosomes( segs_tab, gender, verbose=verbose )

  # if( !is.na(gender) && gender %in% c( "Male", "Female") )
  # {
  #    nix = segs_tab[,"Chromosome"] %in% c("Y", "M", "chrY", "chrM")
  # }
  # else
  # {
  #    if( verbose ) {
  #      print("No or invalid gender specified - dropping X chromosome segs")
  #    }
  #    nix = segs_tab[,"Chromosome"] %in% c("X", "Y", "M", "chrX", "chrY", "chrM")
  # }
  seg_dat$gender = gender
  # segs_tab = segs_tab[!nix, ]



 ## NEED this for mm9 
  # segs_tab[,"Chromosome"] = as.integer(gsub( "chr", "", segs_tab[,"Chromosome"]))
  # if( any(is.na(segs_tab[,"Chromosome"]))) { stop("converted to non-integer chromosome") }



  # segtab = segs_tab[, c("Chromosome", "Start", "End", "Num_Probes")]
  # colnames(segtab) = c("Chromosome", "Start.bp", "End.bp", "n_probes")
  segtab =  segs_tab[, c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", "tau", "sigma.tau")]
  colnames(segtab) = c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", "tau", "seg_sigma")
  
  length = segtab[, "End.bp"] - segtab[, "Start.bp"]
  ## Convert from base 2 log
  # copy_num = 2^(segs_tab[, "tau"] )
  copy_num = segs_tab[, "tau"] / 2
  
  ix = copy_num > 25.0
  if (verbose) {
    print( paste( "Capping ", sum(ix), " segs at tCR = 5.0", sep=""))
  }
  copy_num[ix] = 25.0
  
  seg_sigma_num = 0.1  ## TODO - get rid of this - not used in downstream model - but crashes filtering code if missing
  seg_sigma =  seg_sigma_num / sqrt(as.numeric(segs_tab[,"n_probes"]))
#  seg_sigma = rep(NA, nrow(segs_tab))  ## calculate later in SCNA_model
  seg.ix <- cbind(c(1:nrow(segtab)))
  colnames(seg.ix) = c("seg.ix")

  segtab = cbind(segtab, length, copy_num, seg_sigma, seg.ix)

  if (filter_segs) {
    seg_dat$segtab = FilterSegs(segtab, min_probes=min_probes, max_sd=max_sd)$seg.info
  }

  W = as.numeric(seg_dat$segtab[,"length"])
  W = W / sum(W)
  seg_dat$segtab = cbind(seg_dat$segtab, W)
  colnames(seg_dat$segtab)[ncol(seg_dat$segtab)] = "W"

  seg_dat$error_model = list()

## create an object for total CN analysis
  total.seg.dat =  segtab[, c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", "tau", "seg_sigma", "seg.ix")]
  colnames(total.seg.dat) = c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", "copy_num", "seg_sigma", "seg.ix")

  W <- as.numeric(total.seg.dat[, "length"])
  W <- W / sum(W)
  total.seg.dat <- cbind(total.seg.dat, W)
  seg_dat$total.seg.dat = total.seg.dat

  return(seg_dat)
}

total_extract_sample_obs = function(seg.obj)
{
  seg.tab = seg.obj[["segtab"]]
  d = seg.tab[, "copy_num"]
  stderr <- seg.tab[, "seg_sigma"]
  W <- seg.tab[, "W"]

  if( "bi.allelic" %in% colnames(seg.tab) ) {
     bi.allelic = seg.tab[, "bi.allelic"]
  } else {
     bi.allelic = rep( FALSE, nrow(seg.tab))
  }

  ## expected copy-number, should be 1.0
  e.cr = sum(W * d )
  gender=seg.obj$gender

  X.ix = (seg.tab[,"Chromosome"]=="X")
  Y.ix = (seg.tab[,"Chromosome"]=="Y")
  normal_allele_count = rep(2, nrow(seg.tab))

  if( !is.na(gender) && gender == "Male" ) 
  {
     normal_allele_count[X.ix] = 1 
     normal_allele_count[Y.ix] = 1
  }
  if( !is.na(gender) && gender == "Female" ) 
  {
     normal_allele_count[X.ix] = 2
     normal_allele_count[Y.ix] = 0
  }

  ## FIXME: "error.model" was originally named "HSCN_params" - double check this
  obs = list(
    d=d, d.tx=d, d.stderr=stderr, W=W, seg.ix=seq_along(d), bi.allelic=bi.allelic,
    n_probes=seg.tab[,"n_probes"],
    error.model=seg.obj$error_model, e.cr=e.cr, data.type="TOTAL",
    platform=seg.obj[["platform"]], "normal_allele_count"=normal_allele_count
  )
  
  return(obs)
}

CAPSEG_get_seg_sigma = function(SCNA_model, obs)
{
  seg_sigma =  exp(SCNA_model[["Theta"]]["sigma.A"]) / (sqrt(obs[["n_probes"]] ))   
  sigma.h = SCNA_model[["sigma.h"]]

  seg_sigma = sqrt(seg_sigma^2 + sigma.h^2)
}
