## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GenomeHetFilter <- function(obs, mode.res, max.non.clonal, max.neg.genome,
                            Q, verbose=FALSE) {
  ## calculate provisional seg_Z_tab and filter out modes that imply > 50% genome het.
  ## and filter out modes with > 2.5% het genome < 0
  mode.tab <- mode.res[["mode.tab"]]
  ## init both to all zeros
  frac.het <- frac.neg.het <- rep(0, nrow(mode.tab))

  for (i in seq_len(nrow(mode.tab))) {
    delta <- mode.tab[i, "delta"]
    b <- mode.tab[i, "b"]
    
    obs[["error.model"]][["fit.at"]] <- mode.tab[1, "AT"]
    comb <-  GetCopyRatioComb(Q, delta, b, obs[["error.model"]])
#    seg.z <- mode.res[["seg.qz.tab"]][i, , Q+1]
    seg.z <- mode.res[["mode_SCNA_models"]][[i]][["seg.qz.tab"]][,Q+1]

    frac.het[i] <- sum(seg.z * obs[["W"]])
#    frac.neg.het[i] <- sum((obs[["W"]] * seg.z)[obs[["d.tx"]] < comb[1]])
    frac.neg.het[i] <- sum((obs[["W"]])[ seg.z > 0.9 & obs[["d.tx"]] < comb[1]])
  }

  if (max.non.clonal > 0) {
    nc.ix <- (frac.het > max.non.clonal)
  
    if (verbose) {
      print(paste("removing ", sum(nc.ix), " / ", length(nc.ix),
                  " modes with >", max.non.clonal*100, "% genome non-clonal.", sep=""))
    }
  } else {
    nc.ix <- rep(FALSE, length(frac.het))
  }

if( FALSE )
{
  if (max.neg.genome > 0) {
    neg.mode.ix <- (frac.neg.het > max.neg.genome) & (!nc.ix)

  if( is.na(neg.mode.ix)){ stop() }

    if (verbose) {
      print(paste("removing ", sum(neg.mode.ix), " / ", length(neg.mode.ix),
                  " modes with >", (max.neg.genome * 100) ,
                  "% genome non-clonal < 0 copies.", sep=""))
    }
  } else {
    neg.mode.ix <- rep(FALSE, length(frac.het))
  }
}

  if( any(is.na(nc.ix)) ) { stop() }

  ## return the 'bad' indices
  return(nc.ix)

}





NegGenomeFilter = function( obs, mode.tab, max.neg.genome, Q, verbose=FALSE) 
{
  frac.neg.het <- rep(0, nrow(mode.tab))
  eps = 0.1 ## this much below comb(0)

  for (i in seq_len(nrow(mode.tab))) 
  {
    delta <- mode.tab[i, "delta"]
    b <- mode.tab[i, "b"]
    
#    obs[["error.model"]][["fit.at"]] <- mode.tab[1, "AT"]
    comb_A <-  GetCopyRatioComb(Q, delta, b, obs[["error.model"]])
    comb_X <-  get_male_sex_chr_comb(Q, delta, b, obs)
#    frac.neg.het[i] <- sum((obs[["W"]] * seg.z)[obs[["d.tx"]] < comb[1]])

    one_allele_ix = obs[["normal_allele_count"]] == 1
    frac.neg.het[i] = sum( obs[["W"]][ !one_allele_ix & obs[["d.tx"]] < comb_A[1] - eps ] ) + 
                      sum( obs[["W"]][  one_allele_ix & obs[["d.tx"]] < comb_X[1] - eps ] )
  }

  if (max.neg.genome > 0) {
    neg.mode.ix <- (frac.neg.het > max.neg.genome)

  if( any(is.na(neg.mode.ix))){ stop() }

    if (verbose) {
      print(paste("removing ", sum(neg.mode.ix), " / ", length(neg.mode.ix),
                  " modes with >", (max.neg.genome * 100) ,
                  "% genome non-clonal < 0 copies.", sep=""))
    }
  } else {
    neg.mode.ix <- rep(FALSE, length(frac.neg.het))
  }

  return( neg.mode.ix )
}


ClonalSSNVFilter = function(mode.res, mut.cn.dat, max.Z=0.5, verbose=TRUE)
{
  mode.tab = mode.res[["mode.tab"]]
  n.modes = nrow(mode.tab)

  z <- rep(1, n.modes)
  for (j in 1:n.modes)
  {
    modeled <- mode.res[["modeled.muts"]][[j]]
    mut.dat <- cbind(mut.cn.dat, modeled)

    pr.clonal <- mut.dat[, "Pr_somatic_clonal"]
    n.grid = 100
    beta.grid = GetMutBetaDensities(mut.dat, n.grid=n.grid)
    grid.vals = seq_len(n.grid) / (n.grid + 1)
    clonal.grid <- matrix(NA, nrow=nrow(beta.grid), ncol=ncol(beta.grid))
    for (i in seq_len(nrow(beta.grid))) {
      clonal.grid[i, ] <- beta.grid[i, ] * pr.clonal[i]
    }
    obs.vaf = colSums(clonal.grid) * grid.vals
    obs.vaf.mean = mean(obs.vaf)
    obs.vaf.sd = sd(obs.vaf)
    cv = obs.vaf.sd / obs.vaf.mean

    alpha = mode.tab[j, "alpha"]
    SSNV_skew <- mut.dat[1, "SSNV_skew"]
    modeled.vaf = alpha * SSNV_skew / 2

    z[j] = abs(obs.vaf.mean / modeled.vaf - 1) / cv

    # print(paste(modeled.vaf, obs.vaf.mean, z[j], sep="  "))
  }

  Z = max(1.1 * min(z), max.Z)

  neg.mode.ix <- (z > Z)

  if (verbose) {
    print(paste("removing ", sum(neg.mode.ix), " / ", length(neg.mode.ix),
                " modes with modeled purity disagreeing with mean clonal SSNV VAF distribution (Z > ", round(Z, 3), ").", sep=""))
  }

  return(neg.mode.ix)
}


RealAlphaFilter = function(mode.res, verbose=TRUE)
{
  mode.tab = mode.res[["mode.tab"]]
  neg.mode.ix <- (mode.tab[, "alpha"] == 1.)

  if (all(neg.mode.ix)) {
    neg.mode.ix <- rep(FALSE, length(neg.mode.ix))
  } else {
    if (verbose) {
      print(paste("removing ", sum(neg.mode.ix), " / ", length(neg.mode.ix),
                  " modes with modeled purity == 1.", sep=""))
    }
  }

  return(neg.mode.ix)
}
