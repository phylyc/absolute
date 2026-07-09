## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


ProvisionalModeSweep <- function(seg.obj, SCNA_model, mut.cn.dat, SSNV_model, force.alpha, force.tau, b.res=0.1, d.res=0.01, chr.arms.dat, refit.alpha.tau=FALSE, refit.window.alpha=0.05, refit.window.tau=NA, verbose=FALSE)
{
  Q = SCNA_model[["kQ"]]
  obs <- seg.obj[["obs.scna"]]
  res <- FindLocationModes(obs, force.alpha, force.tau, SCNA_model, mut.cn.dat, SSNV_model, b.res=b.res, d.res=d.res,
                           refit.alpha.tau=refit.alpha.tau, refit.window.alpha=refit.window.alpha, refit.window.tau=refit.window.tau, verbose=verbose)
  
  if (!is.null(res[["mode.flag"]])) { 
      return(list("mode.flag"=res[["mode.flag"]]))
   } else {
     mode.tab <- res[["mode.tab"]]
   }

   if (identical(mode.tab, NA))  {
      return(list(mode.flag="ERROR"))
   }


## Filter modes with lots of genome mass below comb(0) 
   if( !(!is.na(force.alpha) & !is.na(force.tau) & nrow(mode.tab)==1) ) ## Skip this filter in forced-calling mode
   {
      bad.ix = NegGenomeFilter( obs, mode.tab, max.neg.genome, Q, verbose=verbose) 
      mode.tab = mode.tab[!bad.ix,, drop=FALSE]
   }
##

   return( mode.tab )
}



FindLocationModes <- function(obs, force.alpha, force.tau, SCNA_model, mut.cn.dat, SSNV_model, b.res=0.1, d.res=0.01, refit.alpha.tau=FALSE, refit.window.alpha=0.05, refit.window.tau=NA, verbose=FALSE) {

  kAlphaDom <- c(0, 1)

#  mode.tab <- MargModeFinder(obs, mut.cn.dat, SSNV_model, SCNA_model, verbose=verbose)

  if( !is.na(force.alpha) & !is.na(force.tau) )
  {
     if( isTRUE(refit.alpha.tau) )
     {
        ## Don't take the seeds at face value: locally optimize near them.
        mode.row = refit_forced_mode(obs, force.alpha, force.tau, SCNA_model,
                                     window.alpha=refit.window.alpha,
                                     window.tau=refit.window.tau, verbose=verbose)
        mode.tab = matrix( mode.row, nrow=1)
     }
     else
     {
        res = get_b_and_delta( force.alpha, force.tau)
#        mode.tab = rbind( c(res$b, log(res$delta), NA), mode.tab )
        mode.tab = matrix( c(res$b, log(res$delta)), nrow=1)
     }

     if(verbose)
     {
       msg = paste("Added force-call mode (seed) alpha=", force.alpha, ", tau=", force.tau,
                   if (isTRUE(refit.alpha.tau)) " [refit enabled]" else "", sep="")
       print(msg)
     }
  }
  else
  {
     mode.tab <- MargModeFinder(obs, mut.cn.dat, SSNV_model, SCNA_model, b.res=b.res, d.res=d.res, verbose=verbose)
  }
  
  if (nrow(mode.tab) == 0) {
    return(list(mode.flag = "DELTA_B_DOM"))
  }
  
  mode.tab <- cbind(mode.tab, NA, NA, NA)
  for (i in 1:nrow(mode.tab)) {
    par <- mode.tab[i, c(1, 2)]
    cur.par <- par
    mode.params <- list(b = par[1], delta = par[2], AT=NA)
    
    LL <- CombLL(cur.par, obs, SCNA_model )


#    mode.hess <- hessian(CombLL, x = cur.par, method = "Richardson", obs = obs, SCNA_model=SCNA_model )
    
#    if (!is.na(mode.hess[1])) {
#      mode.curv <- CalcModeLogCurv(par, mode.hess, verbose=verbose)
#    } else {
#      LL <- NA
#      mode.curv <- NA
#    }
   mode.curv = NA
    
    mode.tab[i, ] <- c(mode.params[["b"]], mode.params[["delta"]],
                       mode.params[["AT"]], LL, mode.curv)
  }
  
  b <- mode.tab[, 1]
  delta <- exp(mode.tab[, 2])
  at <- mode.tab[, 3]
  res <- GetAlphaAndTau(b, delta)
  alpha <- res[["alpha"]]
  tau <- res[["tau"]]
  
  LL <- mode.tab[, 4]
  mode.curv <- mode.tab[, 5]
  mode.tab <- cbind(alpha, tau, at, b, delta, LL, mode.curv)
  colnames(mode.tab) <- c("alpha", "tau", "AT", "b", "delta", "LL", "mode_curv")

  if (verbose) {
    print(mode.tab)
  }
  
  mode.ix <- (mode.tab[, "alpha"] >= kAlphaDom[1] &
              mode.tab[, "alpha"] <= kAlphaDom[2] & 
              mode.tab[, "tau"] >= SCNA_model[["kTauDom"]][1] &
              mode.tab[, "tau"] <= SCNA_model[["kTauDom"]][2])

  if (verbose) {
    print(paste("removing ", sum(!mode.ix), " / ",
                length(mode.ix), " modes outside of alpha/tau range.", 
                sep = ""))
  }
  mode.tab <- mode.tab[mode.ix, , drop = FALSE]
  
  if (nrow(mode.tab) == 0) {
    return(list(mode.flag = "ALPHA_TAU_DOM"))
  }

  # mode.ix <- (mode.tab[, "alpha"] == 1)
  # if (sum(!mode.ix) > 1) {
  #   if (verbose) {
  #     print(paste("removing ", sum(mode.ix), " / ",
  #                 length(mode.ix), " modes with alpha == 1.",
  #                 sep = ""))
  #   }
  #   mode.tab <- mode.tab[!mode.ix, , drop = FALSE]
  # }
  
  return(list(mode.tab = mode.tab))
}

## Local refit of a user-supplied (alpha, tau) seed. Rather than converting the
## seed straight to (b, delta), run a short bounded optimization of the provisional
## SCNA log-likelihood over a small box around the seed and return the best-fitting
## c(b, log(delta)) pair. Invoked from FindLocationModes when --refit-alpha-tau is set.
refit_forced_mode <- function(obs, force.alpha, force.tau, SCNA_model,
                              window.alpha=0.05, window.tau=NA, verbose=FALSE)
{
  ## When the tau window is not supplied (NA), derive it from the alpha window so
  ## both grant the copy-ratio comb the same wiggle. The comb spacing delta obeys
  ## 1/delta = tau + 2(1-alpha)/alpha, so d(1/delta)/dtau = 1 and
  ## d(1/delta)/dalpha = -2/alpha^2; matching the change in 1/delta gives
  ## window.tau = (2 / alpha^2) * window.alpha. This grows as purity drops, so cap
  ## it at 1.0 to keep low-purity samples from ranging over implausible ploidy.
  tau.window.auto <- is.na(window.tau)
  if (tau.window.auto) {
    window.tau <- min( (2 / force.alpha^2) * window.alpha, 1.0 )
    if (!is.finite(window.tau)) { window.tau <- 1.0 }
  }

  ## Negative provisional SCNA LL as a function of (alpha, tau); convert to
  ## (b, delta) for the LL. provisional_SCNA_LL() stop()s on non-finite values,
  ## so trap and penalize.
  neg_ll <- function(alpha, tau) {
    bd <- get_b_and_delta(alpha, tau)
    ll <- tryCatch(provisional_SCNA_LL(obs, bd[["b"]], bd[["delta"]], SCNA_model),
                   error = function(e) NA_real_)
    if (!is.finite(ll)) { return(1e12) }
    return(-ll)
  }

  ## Box = seed +/- window, clamped to the valid alpha and tau domains.
  alpha.lo <- max(1e-4,     force.alpha - window.alpha)
  alpha.hi <- min(1 - 1e-4, force.alpha + window.alpha)
  tau.lo   <- max(SCNA_model[["kTauDom"]][1], force.tau - window.tau)
  tau.hi   <- min(SCNA_model[["kTauDom"]][2], force.tau + window.tau)

  ## Guard against a degenerate box when the seed sits outside its domain.
  if (alpha.hi < alpha.lo) { alpha.lo <- alpha.hi <- min(max(force.alpha, 1e-4), 1 - 1e-4) }
  if (tau.hi   < tau.lo)   { tau.lo   <- tau.hi   <- min(max(force.tau, SCNA_model[["kTauDom"]][1]), SCNA_model[["kTauDom"]][2]) }

  ## A zero (or negative, or clamped-away) window pins that axis. L-BFGS-B
  ## degenerates when a coordinate has lower==upper (it makes no progress on the
  ## free axis either), so only optimize the free axes: 2-D via L-BFGS-B, 1-D via
  ## optimize(), or 0-D (return the seed) when both are pinned.
  free.alpha <- alpha.hi > alpha.lo
  free.tau   <- tau.hi   > tau.lo

  best.alpha <- min(max(force.alpha, alpha.lo), alpha.hi)
  best.tau   <- min(max(force.tau,   tau.lo),   tau.hi)

  if (free.alpha && free.tau) {
    opt <- tryCatch(
      optim(par = c(best.alpha, best.tau),
            fn = function(par) neg_ll(par[1], par[2]), method = "L-BFGS-B",
            lower = c(alpha.lo, tau.lo), upper = c(alpha.hi, tau.hi),
            control = list(maxit = 200)),
      error = function(e) {
        if (verbose) { print(paste("Refit optim failed; falling back to seed. Error:", conditionMessage(e))) }
        NULL })
    if (!is.null(opt)) { best.alpha <- opt[["par"]][1]; best.tau <- opt[["par"]][2] }
  } else if (free.tau) {
    best.tau   <- optimize(function(tau)   neg_ll(best.alpha, tau), lower = tau.lo,   upper = tau.hi)[["minimum"]]
  } else if (free.alpha) {
    best.alpha <- optimize(function(alpha) neg_ll(alpha, best.tau), lower = alpha.lo, upper = alpha.hi)[["minimum"]]
  }
  ## else: both axes pinned -> keep the (clamped) seed

  bd <- get_b_and_delta(best.alpha, best.tau)

  if (verbose) {
    print(sprintf("Refit alpha %.4f -> %.4f, tau %.4f -> %.4f (LL=%.4f, window +/-%.3f alpha / +/-%.3f tau%s)",
                  force.alpha, best.alpha, force.tau, best.tau, -neg_ll(best.alpha, best.tau),
                  window.alpha, window.tau,
                  if (tau.window.auto) " [auto from alpha]" else ""))
  }

  return(c(bd[["b"]], log(bd[["delta"]])))
}

MargModeFinder <- function(obs, mut.cn.dat, SSNV_model, SCNA_model, b.res=0.1, d.res=0.01, verbose=FALSE)
{
   fit_func = function(b.grid, d.grid, i, j, obs, SCNA_model)
   {
      cur.par <- c(b.grid[i], d.grid[j])
      val = RunOpt(cur.par, obs, SCNA_model, verbose=verbose) 
      res =  c(val[[1]], val[[2]], val[[3]])
      return(res)
   }

  # b.grid <- seq( SCNA_model[["kDom1"]][1], SCNA_model[["kDom1"]][2], b.res)  # b.res=0.125
  # d.grid <- seq( SCNA_model[["kDom2"]][1], SCNA_model[["kDom2"]][2], d.res)  # d.res = 0.125

  b.grid <- log( seq(from=exp(SCNA_model[["kDom1"]][1]), to=exp(SCNA_model[["kDom1"]][2]), by=b.res) )
  d.grid <- log( seq(from=exp(SCNA_model[["kDom2"]][1]), to=exp(SCNA_model[["kDom2"]][2]), by=d.res) )

  n.b <- length(b.grid)
  n.d <- length(d.grid)

  if (verbose) {print(paste("Provisional mode sweep across", n.b * n.d, "b-delta grid points:"))}
  
  mode.tab <- array(NA, dim = c(n.b * n.d, 3))
  #  for (i in seq_len(n.b)) 
  opt_res = foreach( i = 1:n.b ) %dopar%
  {
    res= matrix( NA, nrow=n.d, ncol=3)
    for (j in seq_len(n.d)) 
    {
      #       val = RunOpt(cur.par, obs, SCNA_model, verbose=verbose) 
      #       res[j,] =  c(val[[1]], val[[2]], val[[3]])
      res[j,] = fit_func( b.grid, d.grid, i, j, obs, SCNA_model )
    }
    return(res)
  }
  if (verbose) { cat("\n") }
  
  for( i in 1:n.b )
  {
    res = opt_res[[i]] 
    #       cur.par <- c(b.grid[i], d.grid[j])
    #       res <- RunOpt(cur.par, obs, SCNA_model, verbose=verbose)
    if (!identical(res, NA))
    {
      #       mode.tab[(i - 1) * n.d + j, ] <- c(res[[1]], res[[2]], res[[3]])
      for( j in seq_len(n.d)) {
        mode.tab[(i - 1) * n.d + j, ] = res[j,]
      }
    }
  }

  ## Try 1d opt for pure tumors
  delta_dom = log(c(1 / SCNA_model[["kTauDom"]][2] - 0.05, 1))
  res_1d = run_1d_opt(obs, SCNA_model, delta_dom, d_res=d.res, verbose=verbose)
  if (verbose) {
    print("1d mode opt: ")
    print(res_1d)
  }
  mode.tab = rbind(mode.tab, res_1d)
    
  # try opt on SNVs only for diploid tumors
  if ( FALSE &   !identical(mut.cn.dat, NA)) {
    alpha_dom = c(0.1, 1) 
    res_snv_only = run_diploid_snv_purity_opt(obs, mut.cn.dat, SSNV_model, alpha_dom, verbose=verbose)
    
    if (!is.na(res_snv_only)) {
      if (verbose) {
        print("SNV only opt: ")
        print(res_snv_only)
      }
      mode.tab = rbind(mode.tab, res_snv_only)
    }
  }
    
  ## find unique modes in table
  ix <- !is.na(mode.tab[, 1])
  mode.list <- mode.tab[ix, c(1, 2), drop = FALSE]
  umodes = unique(round(mode.list, 2))
  
  if (verbose) {
    print(paste(nrow(umodes), " unique modes found", sep = ""))
  }

  keep = umodes[, 1] >= SCNA_model[["kDom1"]][1] &  umodes[, 1] <= SCNA_model[["kDom1"]][2] &
         umodes[, 2] >= SCNA_model[["kDom2"]][1] &  umodes[, 2] <= SCNA_model[["kDom2"]][2] 
  umodes <- umodes[keep, , drop=FALSE ]

  if (verbose) {
    print(paste(nrow(umodes), " modes in b / delta range.", sep = ""))
  }

 ## now make sure the grid of b / delta intersections is filled in
  new.mode.tab = fill_mode_neighborhoods( umodes, SCNA_model )
  if (verbose) { print(paste(nrow(new.mode.tab), " modes after neighborhood closure", sep = "")) }

## Fill in provisional LL of each mode
  LL = rep(NA, nrow(new.mode.tab))
  for( i in 1:nrow(new.mode.tab))
  {
    b = new.mode.tab[i,1]
    delta = exp(new.mode.tab[i,2])
    LL[i] = provisional_SCNA_LL( obs, b, delta, SCNA_model )
  }
  new.mode.tab = cbind( new.mode.tab, LL )


  ht = hclust( dist(new.mode.tab[,c(1,2)]))
  assign = cutree( ht, h=0.2 )
  k = length(unique(assign))
  clust.mode.tab = matrix( NA, nrow=k, ncol=2 )

  clust_ids = unique(assign)
  for( i in 1:k )
  {
    c.ix = which(assign==clust_ids[i])

## pick representative with best provis LL-score for each cluster
    c.best = which.max(new.mode.tab[c.ix,3])
    if(length(c.ix)==1) { c.best=1 }

    clust.mode.tab[i,] = new.mode.tab[ c.ix[c.best],  c(1,2) ]
  }
  if (verbose) { print(paste(nrow(clust.mode.tab), " modes after clustering", sep = "")) }
  colnames(clust.mode.tab) = c("b", "delta")

  return( clust.mode.tab )
}


fill_mode_neighborhoods = function( mode.tab, SCNA_model )
{
## tabs are b, log(delta)  pairs

  add_neighbors = function( mode.tab, SCNA_model )
  {
    new.modes = matrix( NA, ncol=2, nrow=0 )

    for( i in 1:nrow( mode.tab ) )
    {
      b = mode.tab[i, 1] 
      delta = exp( mode.tab[i, 2] )

    ## b/ delta multiples
      delta.set <- c(delta / 2, 2 * delta)
      b.set <- c(b - 2 * delta, b + 2 * delta)
 
      new.tab = rbind( 
                       cbind( b, log(delta.set) ), 
                       cbind( b.set, log(delta) ),
                       cbind( b.set, log(delta.set)),
                       cbind( rev(b.set), log(delta.set)),
                       cbind( b.set, rev(log(delta.set)))
                     )

      new.modes = rbind( new.modes, new.tab )
    }

    umodes = unique(round(new.modes, 2))
    keep = umodes[, 1] >= SCNA_model[["kDom1"]][1] &  umodes[, 1] <= SCNA_model[["kDom1"]][2] &
           umodes[, 2] >= SCNA_model[["kDom2"]][1] &  umodes[, 2] <= SCNA_model[["kDom2"]][2] 

    umodes = umodes[ keep, , drop = FALSE]
    new.tab = rbind(mode.tab, umodes)

    new.tab = unique( new.tab )
    return(new.tab)
  } 


  n_modes = nrow(mode.tab)
  while( 1 )
  {   
#     cat("+")
     mode.tab = add_neighbors( mode.tab, SCNA_model )

     new_n_modes = nrow( mode.tab ) 

     break
#     if( new_n_modes == n_modes ) { break }
     n_modes = new_n_modes
  }
  cat("\n")
 
  return( mode.tab )
}


RunOpt <- function(cur.par, obs, SCNA_model, eval.hessian=FALSE, verbose=FALSE) {
  LL <- CombLL(cur.par, obs = obs, SCNA_model = SCNA_model )
    
  if (is.finite(LL)) {
    res <- optim(par = cur.par, fn = CombLL, gr = NULL, method = "Nelder-Mead", 
                 obs = obs, SCNA_model=SCNA_model, control = list(maxit = 1000),
                 hessian = eval.hessian)
    res <- list(res[["par"]][1], res[["par"]][2], -res[["value"]],
                res[["hessian"]])
    if (verbose) {
      cat(".")
    }
  } else {
    if (verbose) {
      cat("!")
      stop()
    }
    res <- NA
  }
  
  return(res)
}

run_1d_opt = function(obs, SCNA_model, dom, d_res, verbose=FALSE) 
{
  comb_1d_ll = function(par, obs, SCNA_model, dom) {

    delta = par

    if( !is.finite(delta) ) { stop("Non-finite delta!") }

#    if(delta<dom[1] | delta>dom[2]) { return( Inf ) }
    if( delta<dom[1] ) { delta = dom[1] }
    if( delta>dom[2] ) { delta = dom[2] }

    delta = exp(delta)
    b = 0
#    LL = SCNA_model_loglik(obs, b, delta, SCNA_model, Provisional=TRUE )
    LL = provisional_SCNA_LL( obs, b, delta, SCNA_model )
    
    if( !is.finite(LL)) { stop() }

    return(-LL)
  }
  
  d_grid = seq(dom[1], dom[2], d_res)
  mode_tab = array(NA, dim=c( length(d_grid), 3))
  
  for (i in seq_along(d_grid)) {
    opt = try( nlm(f=comb_1d_ll, p=d_grid[i], obs=obs, SCNA_model=SCNA_model, dom=dom) )
    if( class(opt)=="try-error" ) { next }
    
    mode_tab[i, ] = c(0, opt$estimate, -opt$minimum)
    if (verbose) {
      cat("-")
    }
  }

  n.ix = is.na(mode_tab[,1])
  mode_tab = mode_tab[ !n.ix,, drop=FALSE ]

  if (verbose) {
    cat("\n")
  }
  
  return(mode_tab)
}

CombLL <- function(par, obs, SCNA_model )
{
  if (any(is.na(par))) {
      cat("$")
      return(Inf)
  }
    
  b <- par[1]
  delta <- exp(par[2])

#  alpha <- par[1]
#  tau <- par[2]

#   res = get_b_and_delta( alpha, tau )
#   b = res[["b"]]
#   delta = res[["delta"]]

#  print(str(obs))
#  print(str(b))
#  print(str(delta))
#  print(str(SCNA_model))
  
  LL = provisional_SCNA_LL( obs, b, delta, SCNA_model )
    
  if( !is.finite(LL)) { stop("Not a finite LL") }

  return(-LL)
}

provisional_SCNA_LL = function( obs, b, delta, SCNA_model )
{
   clonal_seg_LL_mat = get_clonal_seg_LL_mat( obs, b, delta, SCNA_model )
## don't add outlier prob to segs below comb(0)
if( FALSE )
{
   neg.ix = obs$d.tx < b
   out_LL = rep(log(1e-3), length(obs$d.tx))
   out_LL[neg.ix] = log(1e-5)

   seg_LL = LogAdd( cbind(clonal_seg_LL_mat, out_LL) )   ## add an outlier state.. 
   LL = sum(seg_LL)
} else {
      LL = sum( LogAdd(clonal_seg_LL_mat) )
   }

   return(LL)
}

get_b_and_delta = function( alpha, tau )
{
  delta <- alpha / (2 * (1 - alpha) + alpha * tau)
  b <- 2 * (1 - alpha) / (2 * (1 - alpha) + alpha * tau)

  return(list("b"=b, "delta"=delta))
}


GetAlphaAndTau <- function(b, delta) {
   alpha = (2 * delta) / (2 * delta + b)
   tau = -(b - 1) / delta

   return(list(alpha=alpha, tau=tau)) 
}


