filter_sex_chromosomes = function( seg.dat, gender, verbose )
{
  if( !is.na(gender) && gender %in% c( "Male", "Female") )
  {
     nix = seg.dat[,"Chromosome"] %in% c("Y", "M", "chrY", "chrM")
#     nix = !(seg.dat[,"Chromosome"] %in% c(c(1:22),"X") )
  }
  else
  {
     if( verbose ) {
       print("No or invalid gender specified - dropping X chromosome segs")
     }
     nix = seg.dat[,"Chromosome"] %in% c("X", "Y", "M", "chrX", "chrY", "chrM")
  }
  seg.dat = seg.dat[!nix, ]

  return(seg.dat)
}


AllelicMakeSegObj <- function(seg.dat, gender, filter_segs=FALSE, min_probes=NA, max_sd=NA, verbose=FALSE) 
{   
#  print("Gender not supported yet!!")
#  gender = NA

  X.ix = seg.dat[,"Chromosome"]==23
  seg.dat[X.ix, "Chromosome"] = "X"

## normalize gender col
  is.M = gender %in% c("M", "male", "Male", "XY")
  is.F = gender %in% c("F", "female", "Female", "XX")

  if(is.M) { gender = "Male" }
  if(is.F) { gender = "Female" }

  if( verbose ) { print(paste("Detected ", gender, " gender.", sep="")) }
  


## Filter out sex-appropriate chromosomes
  seg.dat = filter_sex_chromosomes( seg.dat, gender, verbose=verbose )

  nix= is.na(seg.dat[,"tau"]) | is.na(seg.dat[,"sigma.tau"])
  print( paste( "Removing ", sum(nix), " of ", length(nix), " segs with NA tCR or tCR sem", sep="") )
  seg.dat = seg.dat[!nix,]

## TODO: configure optional filtering on tCR seg-size and SD
#  if (filter_segs) {
#     seg.dat = FilterSegs(seg.dat, min_probes=min_probes, max_sd=max_sd,
#                                  verbose=verbose)$seg.info
#  }


## create an object for total CN analysis
  total.seg.dat =  seg.dat[, c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", "tau", "sigma.tau")]
  colnames(total.seg.dat) = c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", "copy_num", "seg_sigma")
    
  W <- as.numeric(total.seg.dat[, "length"])
  W <- W / sum(W)
  total.seg.dat <- cbind(total.seg.dat, W)
#  colnames(total.seg.dat)[ncol(total.seg.dat)] <- "W"
##

## Before active seg.dat filter triggered by NA allelic seg.dat fields
  hom.seg.pair.tab = matrix( NA, nrow=nrow(seg.dat), ncol=2 )


  if( !is.na(gender) & gender == "Male" )
  {
     print( paste( "Setting ", sum(X.ix), " of ", length(X.ix), " allelic seg data to NA in male", sep=""))
     X.ix = seg.dat[,"Chromosome"]=="X"
     seg.dat[X.ix,c("mu.minor","mu.major")] = NA
  }

## Remove segs with NA allelic CN
  as.nix= is.na(seg.dat[,"mu.minor"]) | is.na(seg.dat[,"mu.major"])
  print( paste( "Removing ", sum(as.nix), " of ", length(as.nix), " segs with NA minor or major allelic CN", sep="") )
  seg.dat = seg.dat[!as.nix,]


# don't remember what this does
#  if (filter_segs) {
#    as.seg.dat = FilterSegs(as.seg.dat, min_probes=min_probes, max_sd=max_sd,
#                                verbose=verbose)$seg.info
#  }


# TODO: get rid of allele.segs
 ## Construct allele.segs object for genome HSCR plotting
#  allele.segs = seg.dat[, c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", "mu.minor", "mu.major", "bi.allelic")]
#  colnames(allele.segs) = c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", "A1.Seg.CN", "A2.Seg.CN", "bi.allelic")
 

 ## reformat into concatenated seg-tab for each allele...
  if( any(colnames(seg.dat)=="Interest") ) {
     colnames(seg.dat)[colnames(seg.dat)=="Interest"] = "bi.allelic"
     seg.dat[,"bi.allelic"] = seg.dat[,"bi.allelic"] == 0 
  } else {
     seg.dat = cbind( seg.dat, "bi.allelic"=FALSE )
  }

  as.1 <- seg.dat[, c("Chromosome", "Start.bp", "End.bp",
                                       "n_probes", "length", "mu.minor", "sigma.minor", "bi.allelic")]

  as.2 <- seg.dat[, c("Chromosome", "Start.bp", "End.bp",
                                       "n_probes", "length", "mu.major", "sigma.major", "bi.allelic")]
  
  colnames(as.1)[6] <- "copy_num"
  colnames(as.2)[6] <- "copy_num"
  
  colnames(as.1)[7] <- "seg_sigma"
  colnames(as.2)[7] <- "seg_sigma"

  ## remember each segment's mate
  seg.ix <- c(1:nrow(as.1))
  as.1 <- cbind(as.1, seg.ix)
  as.2 <- cbind(as.2, seg.ix)

  t.seg.ix = cbind( "seg.ix.1"=1:nrow(as.1), "seg.ix.2"=1:nrow(as.1)+nrow(as.1) )
  hom.seg.pair.tab[ !as.nix, ] = t.seg.ix  ## matches allelic segs and total CR segs

  as.seg.dat <- rbind(as.1, as.2)
  colnames(as.seg.dat)[ncol(as.seg.dat)] <- "seg.ix"
 
  
  W <- as.numeric(as.seg.dat[, "length"])
  W <- W / sum(W)
  as.seg.dat <- cbind(as.seg.dat, W)
  colnames(as.seg.dat)[ncol(as.seg.dat)] <- "W"


  seg.obj = list()
  seg.obj$gender = gender
  seg.obj$as.seg.dat = as.seg.dat
  seg.obj$segtab = as.seg.dat
  seg.obj$total.seg.dat = total.seg.dat
#  seg.obj$allele.segs = allele.segs
  seg.obj$hom.seg.pair.tab = hom.seg.pair.tab

  return(seg.obj)
}


import_AllelicCapseg_data = function( iams.res , gender, min_probes=min_probes, max_sd=max_sd, 
                          filter_segs=filter_segs, verbose=verbose)
{
   segtab = AllelicCapseg::AbsolutePostProcess(iams.res, iams.res[["merged.segtab"]])
   seg.obj = AllelicMakeSegObj( segtab, gender, filter_segs=filter_segs, min_probes=min_probes, max_sd=max_sd, verbose=verbose ) 

   if( "H1.f.post_mat" %in% names(iams.res) )
   {
      seg.obj[["H1.f.post_mat"]] = iams.res[["H1.f.post_mat"]]
      seg.obj[["seg.f.fit"]] = iams.res[["seg.f.fit"]]
   }
   else
   {
      print("Computing ACS wes.f posteriors on the fly - run a newer version of ACS!") 

      Theta = iams.res[["capture.em.fit"]][["Theta"]]
      h.seg.dat = iams.res[["as.res"]][["h.seg.dat"]]
      nseg = length(h.seg.dat[["gh.wes.allele.d"]])
#   ev.H0 = ev.H1 = rep(NA, nseg)
      dx = 1e-2
      x = seq(0, 0.5, by=dx)
      H1.f.post_mat = matrix(NA, nrow=nseg, ncol=length(x) )
#   
      cat("Computing grid posts on f from AllelicCapseg data")
      for( i in 1:nseg )
      {
         alt <- h.seg.dat[["gh.wes.allele.d"]][[i]]["alt", ]
         ref <- h.seg.dat[["gh.wes.allele.d"]][[i]]["ref", ]
         H1.f.post_mat[i,] = AllelicCapseg:::CalcPosteriorOnF(alt, ref, Theta, dx=dx)[2,] 

#      ev.H1[i] = ACS:::GridEstimateOfFNormalizingConstant(alt, ref, Theta, dx=1e-2) 
#      ev.H0[i] = ACS:::CalcCaptureSegAllelicLogLik(alt, ref, f=0.5, Theta )
#
         cat(".")
      }
      cat("done\n")

      H1.f.post_mat = H1.f.post_mat / rowSums(H1.f.post_mat)
#      grid.f.hat = x[apply(H1.f.post_mat, 1, which.max)]

      init.wes.f = AllelicCapseg:::InitCaptureF( h.seg.dat, Theta )
      res = AllelicCapseg:::OptimizeCaptureF( h.seg.dat, init.wes.f, Theta, tol=1e-5, verbose=verbose ) 
      xx = order(abs(res[,"f.hat"] - res[,"H1.f.hat"]), decreasing=TRUE )
      print(res[xx[1:10],])

      seg.obj[["H1.f.post_mat"]] = H1.f.post_mat
      seg.obj[["seg.f.fit"]] = res
   }
 
   return(seg.obj)
}




## Get map between hom seg pairs
get_hom_seg_pair_map = function( seg.ix )
{
   uniq.seg.ix <- unique(seg.ix)
   n.seg <- length(uniq.seg.ix) 
   seg.pair.tab <- array(NA, dim=c( n.seg, 2 ) )
   colnames(seg.pair.tab) = c("seg.ix.1", "seg.ix.2")

   for (i in 1:n.seg) 
   {
     if (sum(seg.ix == uniq.seg.ix[i]) != 2 ) {
       next
     }
     
     pair.ix <- which(seg.ix == uniq.seg.ix[i]) 
     seg.pair.tab[i, ] = pair.ix
   }
   
   return( seg.pair.tab )
}



get_hom_pairs_segtab = function( seg.dat )
{
   obs = AllelicExtractSampleObs(seg.dat)
   segtab = seg.dat[["as.seg.dat"]]

   return( get_hom_pairs_segtab_internal( obs, segtab )  )
}

get_hom_pairs_segtab_internal = function( obs, as.seg.dat )
{
   AS.seg.ix = obs[["AS.seg.ix"]]
   segtab = as.seg.dat

   allele.segs = cbind( segtab[AS.seg.ix[,1], c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", "W") ],
                        "A1.Seg.CN"=obs[["d.tx"]][AS.seg.ix[,1]],
                        "A2.Seg.CN"=obs[["d.tx"]][AS.seg.ix[,2]],
                        "A1.Seg.sem"=obs[["d.stderr"]][AS.seg.ix[,1]],   
                        "A2.Seg.sem"=obs[["d.stderr"]][AS.seg.ix[,2]],
                        AS.seg.ix )

   return(allele.segs)
}


AllelicExtractSampleObs <- function(seg.obj) 
{
  seg.dat <- seg.obj[["as.seg.dat"]]
  d <- seg.dat[, "copy_num"]
  stderr <- seg.dat[, "seg_sigma"]
  W <- seg.dat[, "W"]

  if( "bi.allelic" %in% colnames(seg.dat) ) {
     bi.allelic = seg.dat[, "bi.allelic"]
  } else { 
     bi.allelic = rep( FALSE, nrow(seg.dat))
  }

  seg.ix <- seg.dat[,"seg.ix"]
  AS.seg.ix = get_hom_seg_pair_map( seg.ix )

  e.cr = sum(W*d )
  gender=seg.obj$gender
  
  X.ix = (seg.dat[,"Chromosome"]=="X")
  Y.ix = (seg.dat[,"Chromosome"]=="Y")
  normal_allele_count = rep(2, nrow(seg.dat))

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

## TODO: purge error.model from all obs.   Get it from the parent seg.dat instead
  obs <- list(d=d, d.tx=d, d.stderr=stderr, W=W, seg.ix=seg.ix, AS.seg.ix=AS.seg.ix, bi.allelic=bi.allelic,
              error.model=list(), n_probes=seg.dat[,"n_probes"], e.cr=e.cr,
	      platform=seg.obj[["platform"]], "normal_allele_count"=normal_allele_count )
  
  return(obs)
}



extract_total_copy_ratios_from_allelic_CAPSEG <- function(seg.obj) 
{
  seg.dat <- seg.obj[["total.seg.dat"]]
  d <- seg.dat[, "copy_num"]/2
  stderr <- seg.dat[, "seg_sigma"]
  W <- seg.dat[, "W"]

  e.cr = sum(W*d )
  gender=seg.obj$gender
  
  X.ix = (seg.dat[,"Chromosome"]=="X")
  Y.ix = (seg.dat[,"Chromosome"]=="Y")
  normal_allele_count = rep(2, nrow(seg.dat))

  if( !is.na(gender) && gender == "Male" ) 
  {
     normal_allele_count[X.ix] = 1 
     normal_allele_count[Y.ix] = 1
     d[X.ix] = d[X.ix]/2  ## correct for the way X-chr CRs are calculated in Males
  }
  if( !is.na(gender) && gender == "Female" ) 
  {
     normal_allele_count[X.ix] = 2
     normal_allele_count[Y.ix] = 0
  }



## TODO: purge error.model from all obs.   Get it from the parent seg.dat instead
  obs <- list(d=d, d.tx=d, d.stderr=stderr, W=W, hom.seg.pair.tab=seg.obj[["hom.seg.pair.tab"]],
              error.model=list(), n_probes=seg.dat[,"n_probes"], e.cr=e.cr,
	      platform=seg.obj[["platform"]], "normal_allele_count"=normal_allele_count )
  
  return(obs)
}



## In fact - this is the right function whenever segs have input std errors
Allelic_get_seg_sigma = function(SCNA_model, obs)
{
  seg_sigma =  exp(SCNA_model[["Theta"]]["sigma.A"]) * obs[["d.stderr"]] 
  sigma.h = SCNA_model[["sigma.h"]]

  seg_sigma = sqrt(seg_sigma^2 + sigma.h^2)
}





