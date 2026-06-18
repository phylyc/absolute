
fit_modes_SCNA_models = function( seg.obj, mode.tab, SCNA_model, mut.cn.dat, chr.arms.dat, verbose=FALSE )
{
  Q = SCNA_model[["kQ"]]
  n.modes = nrow(mode.tab)
  obs <- seg.obj[["obs.scna"]]

  if (verbose) {
    print(paste0("Optimizing SCNA_model | comb) for ", n.modes, " modes: "))
  }
  
## TODO: absorb into SCNA_model
  log_ccf_dens = array(NA, dim=c(n.modes, length(obs[["W"]]), length(SCNA_model[["ccf_grid"]])) )
  dimnames(log_ccf_dens)[[3]] = SCNA_model[["ccf_grid"]]
  subclonal_scna_tab = array(NA, dim=c(n.modes, length(obs[["W"]]), 7) )
  dimnames(subclonal_scna_tab)[[3]] = c("CCF_hat", "CI95_low", "CI95_high", "subclonal_ix", "Pr_subclonal", "qs", "qc" )

  N_tot_seg = length(seg.obj$obs.total.scna$d)
  total_subclonal_scna_tab = array(NA, dim=c(n.modes, N_tot_seg, 7) )
  dimnames(subclonal_scna_tab)[[3]] = c("CCF_hat", "CI95_low", "CI95_high", "subclonal_ix", "Pr_subclonal", "qs", "qc" )
  total_log_ccf_dens = array(NA, dim=c(n.modes, N_tot_seg, length(SCNA_model[["ccf_grid"]])) )
  dimnames(log_ccf_dens)[[3]] = SCNA_model[["ccf_grid"]]



   

  new_cols = c("genome mass", "sigma.h.hat", "theta.z.hat", "sigma.A.hat", "theta.Q.hat", "lambda.hat", "theta.0", "frac.het", "SCNA_LL", "SCNA_min_chrarm_events", "entropy", "Kar_LL", "WGD", "combined_LL", "SSNV_LL", "SCNA_Theta_integral" )

  old_cols = colnames(mode.tab)
  mode.tab <- cbind(mode.tab, matrix(NA, nrow=nrow(mode.tab), ncol=length(new_cols))) 
  colnames(mode.tab) = c( old_cols, new_cols ) 
##
  mode_SCNA_models = list()

## get genome ev params:
  WGD0_Prs = get_WGD0_CN_prob_vector(SCNA_model)
  WGD1_Prs = get_WGD1_CN_prob_vector(SCNA_model)

  init = SCNA_model_init(SCNA_model)

  if( FALSE )  ## DEBUG
  {
     modes_DP_res = list()
     for( i in 1:n.modes )
     {
       if(verbose) {
         cat("\n")
          print( paste("Optimizing PP mode #", i, sep=""))
       }
        modes_DP_res[[i]] = SCNA_model_calc_CCF_DP_loglik( obs, mode.tab[i, "b"], mode.tab[i, "delta"], init, verbose=verbose )
     }
  }
  else
  {
     modes_DP_res = foreach (i = 1:n.modes) %dopar%
     {
       if(verbose) {
         cat("\n")
          print( paste("Optimizing PP mode #", i, sep=""))
       }

       res = SCNA_model_calc_CCF_DP_loglik( obs, mode.tab[i, "b"], mode.tab[i, "delta"], init, verbose=verbose )
       return(res)
     }
  }

  ## fill in other fields 
  for (i in 1:n.modes) 
  {
    delta <- mode.tab[i, "delta"]
    b <- mode.tab[i, "b"]

#    mode_SCNA_models[[i]] = SCNA_model_calc_CCF_DP_loglik( obs, b, delta, init, verbose=verbose )
    mode_SCNA_models[[i]] = modes_DP_res[[i]]
    mode_SCNA_models[[i]] = calc_mode_seg_tabs( seg.obj, mode_SCNA_models[[i]], b, delta, chr.arms.dat )

## Annotate SCNA clonality summary for SSNV models
## Note - these functions only use the mut.cn.dat to disallow clonal homozygous calls if there is 1 > alt SSNV read in the seg
    if (seg.obj[["copy_num_type"]] == "allelic") {
      if( !identical(mut.cn.dat, NA))
      {
         missing.AS.seg.ix = apply( is.na(mut.cn.dat[,c("A1.ix", "A2.ix")]), 1, any )
         if( any(!missing.AS.seg.ix))
         {
            allelic_res = allelic_get_subclonal_SCNA_info( obs, b, delta, mode_SCNA_models[[i]], mut.cn.dat[!missing.AS.seg.ix,, drop=FALSE] )
         }
      }
      else
      {
         allelic_res = allelic_get_subclonal_SCNA_info( obs, b, delta, mode_SCNA_models[[i]], mut.cn.dat )
      }

      subclonal_scna_tab[i,,] = allelic_res[["subclonal_scna_tab"]]
      log_ccf_dens[i,,] = allelic_res[["log_ccf_dens"]]

    ## compute ev score/data
      mode_SCNA_models[[i]][["DP_CN_chrarm_states"]] = get_post_DP_chrarm_states( seg.obj, mode_SCNA_models[[i]], subclonal_scna_tab[i,,], chr.arms.dat )
      mode_SCNA_models[[i]][["SCNA_minev_chrarm_result"]] = compute_chrarm_ev_score( mode_SCNA_models[[i]][["DP_CN_chrarm_states"]], mode_SCNA_models[[i]],  WGD0_Prs, WGD1_Prs )
    }
    else
    {
      ## total CR chr-arm minimum-event (parsimony) score from arm-level modal total CN.
      ## This is the total-CN analog of the allelic event score and is what lets
      ## WeighSampleModes break the WGD/ploidy degeneracy for total CR.
      arm_CN = total_get_chrarm_modal_CN( mode_SCNA_models[[i]][["chr.arm.tab"]] )
      frac_het = sum( obs[["W"]] * mode_SCNA_models[[i]][["seg.z.tab"]] )
      mode_SCNA_models[[i]][["SCNA_minev_chrarm_result"]] = total_compute_chrarm_ev_score( arm_CN, frac_het, mode_SCNA_models[[i]] )
    }

    # allelic_get_subclonal_SCNA_info sets "tot.xxx" states which are needed going forward.
    # In allelic mode total_get_subclonal_SCNA_info reads those tot.* fields; in total CR
    # mode it reads the plain (total) fields (see total_get_subclonal_SCNA_info).
    total_res = total_get_subclonal_SCNA_info( seg.obj[["obs.total.scna"]], b, delta, mode_SCNA_models[[i]], mut.cn.dat )

#    res = get_subclonal_SCNA_info( obs, b, delta, mode_SCNA_models[[i]], mut.cn.dat )

    total_subclonal_scna_tab[i,,] = total_res[["subclonal_scna_tab"]]
    total_log_ccf_dens[i,,] = total_res[["log_ccf_dens"]]

    ## In total CR mode the canonical (primary) subclonal slots hold the total-CN results, so
    ## downstream consumers (plots, extraction, SSNV fit) work through a single uniform slot.
    if (seg.obj[["copy_num_type"]] == "total") {
      subclonal_scna_tab[i,,] = total_res[["subclonal_scna_tab"]]
      log_ccf_dens[i,,] = total_res[["log_ccf_dens"]]
    }

    SCNA_model[["WGD"]] =  SCNA_model[["SCNA_minev_chrarm_result"]][["WGD"]]  ## override provisional estimate
    mode.tab = fill_mode.tab_row( seg.obj, mode_SCNA_models[[i]], obs, b, delta, mode.tab, i )
  }

  subclonal_SCNA_res = list(
    subclonal_SCNA_tab = subclonal_scna_tab, log_CCF_dens = log_ccf_dens,
    total_subclonal_scna_tab = total_subclonal_scna_tab, total_log_ccf_dens = total_log_ccf_dens
  )

  return( list(mode.tab=mode.tab, mode_SCNA_models=mode_SCNA_models, subclonal_SCNA_res=subclonal_SCNA_res, mode.flag=NA ) )
}




WeighSampleModes <- function(mode.res) 
{
## combined various scores
  mode.tab = mode.res[["mode.tab"]]

  if (!all(is.na(mode.tab[,"SCNA_min_chrarm_events"]))) {
    ## Only use SCNA ev score
    mode.res[["mode.tab"]][, "combined_LL"] = mode.tab[,"SCNA_min_chrarm_events"]
  } else {
    ## Only use SCNA LL score
    mode.res[["mode.tab"]][, "combined_LL"] = mode.tab[,"SCNA_LL"] # + mode.tab[,"Kar_LL"] + mode.tab[,"SSNV_LL"]
  }

  LL = mode.res[["mode.tab"]][, "combined_LL"]
  if( !all( is.finite(LL)) ) { stop("Non-finite mode combined_LL!") }

  dens = exp(LL - LogAdd(LL))
  mode.res[["mode.tab"]] <- cbind(mode.res[["mode.tab"]], dens)
  ix = order(LL, decreasing=TRUE )
  
  mode.res <- ReorderModeRes(mode.res, ix)
  
  return(mode.res)
}

fill_mode.tab_row = function( seg.obj, SCNA_model, obs, b, delta, mode.tab, i )
{
   Q = SCNA_model[["kQ"]]
## This is leftover from ancient code - copies fields from fit SCNA_model into mode.tab.   Would it be easier to just save a list of SCNA_model objects??
    Theta_hat = SCNA_model_Theta_tx(SCNA_model[["Theta"]])  ## recover natural units
    mode.tab[i, "theta.0"] <- SCNA_model[["theta.0"]]
    mode.tab[i, "theta.Q.hat"] <- Theta_hat["theta.Q"]
    mode.tab[i, "sigma.A.hat"] <- Theta_hat["sigma.A"]
    mode.tab[i, "sigma.h.hat"] <- SCNA_model[["sigma.h"]]
    mode.tab[i, "lambda.hat"] <- Theta_hat["lambda"]
    mode.tab[i, "SCNA_LL"] <- SCNA_model[["LL"]]
    mode.tab[i, "SCNA_min_chrarm_events"] = SCNA_model[["SCNA_minev_chrarm_result"]][["score"]]

#    mode.tab[i, "theta.z.hat"] <-  sum(SCNA_model[["mix.w"]][c("unif","exp")])
    mode.tab[i, "theta.z.hat"] = sum( SCNA_model[["seg.qz.tab"]][,(Q + 1)] * obs$W)
    ## compute % non-clonal genome
    frac.het = sum(obs[["W"]] * SCNA_model[["seg.z.tab"]])  
    if( !is.finite(frac.het)) { stop() }

    mode.tab[i, "frac.het"] <- frac.het
    ## calculate allelic-balance
    if( seg.obj[["copy_num_type"]] == "allelic" )  {
      mode.tab[i,"genome mass"] <- 2 * sum(c((1:Q)-1) * colSums( SCNA_model[["seg.q.tab"]] * obs[["W"]]))
    }
    if( seg.obj[["copy_num_type"]] == "total" ) {
      mode.tab[i,"genome mass"] <- 1 * sum(c((1:Q)-1) * colSums( SCNA_model[["seg.q.tab"]] * obs[["W"]]))
    }

    mode.tab[i, "WGD"] = SCNA_model[["WGD"]]

    ## weighted entropy average over segs
    mode.tab[i, "entropy"] <- CalcFitEntropy(obs, SCNA_model[["seg.qz.tab"]] )

    return(mode.tab)
}

ReorderModeRes <- function(mode.res, ix, DROP=FALSE) 
{
   mode.res[["mode_SCNA_models"]] = mode.res[["mode_SCNA_models"]][ix] 
   mode.res[["mode.tab"]] <- mode.res[["mode.tab"]][ix,, drop=DROP]
   mode.res[["mode.clust.p"]] <- mode.res[["mode.clust.p"]][ix , , drop=DROP]  ## From Kar model

# ? does this crash with no MAF??
   mode.res[["subclonal_SCNA_res"]][["subclonal_SCNA_tab"]] = mode.res[["subclonal_SCNA_res"]][["subclonal_SCNA_tab"]][ix, , , drop=DROP]
   mode.res[["subclonal_SCNA_res"]][["log_CCF_dens"]] = mode.res[["subclonal_SCNA_res"]][["log_CCF_dens"]][ix, , , drop=DROP]

   ## keep the parallel total-CN tables aligned with their modes too
   mode.res[["subclonal_SCNA_res"]][["total_subclonal_scna_tab"]] = mode.res[["subclonal_SCNA_res"]][["total_subclonal_scna_tab"]][ix, , , drop=DROP]
   mode.res[["subclonal_SCNA_res"]][["total_log_ccf_dens"]] = mode.res[["subclonal_SCNA_res"]][["total_log_ccf_dens"]][ix, , , drop=DROP]

   ## only exists if MAF supplied.
   if (!is.null(mode.res[["modeled.muts"]])) 
   {
      mode.res[["SSNV.ccf.dens"]] = mode.res[["SSNV.ccf.dens"]][ix,,,drop=DROP]
      mode.res[["modeled.muts"]] <- mode.res[["modeled.muts"]][ix]
      mode.res[["mode_SSNV_models"]] = mode.res[["mode_SSNV_models"]][ix]
   }
   mode.res[["mode.posts"]] <- mode.res[["mode.posts"]][ix]

   return(mode.res)
}


## dont call non-aneuploid if mut_dat is present
GetCallStatus <- function(mode.res, seg.w) {
  status <- "called"
  q.tab <- mode.res[["mode_SCNA_models"]][[1]][["seg.qz.tab"]]
  q.tab <- q.tab[, c(1:  (ncol(q.tab)) )]
  
  Q <- ncol(q.tab) - 1
  max.q <- (apply(q.tab, 1, which.max))

  peak_masses <- rep(0, Q)
  for (i in 1:Q) {
    ix <- which(max.q==i) 
    peak_masses[i] <- sum(seg.w[ix])
  }
  
  b <- mode.res[["mode.tab"]][1,"b"]
  
  ## don't count 0-CN state if b is too small - could be due to IBD, not LOH
  if (peak_masses[1] < 0.01 & b < 0.15) {
    peak_masses[1] <- 0
  }

  six <- order(peak_masses, decreasing=TRUE)

  if (peak_masses[six[3]] < 0.0001) {
    if (peak_masses[six[2]] < 0.0001) { 
      if (mode.res[["mode.tab"]][1,"sigma.h.hat"] < 0.02 &
          is.null(mode.res[["modeled.muts"]])) { 
        status <- "non-aneuploid" 
      }
      if (mode.res[["mode.tab"]][1,"sigma.h.hat"] >= 0.02 ) {
        status <- "low purity"
      }
    }
  }

   if (mode.res[["mode.tab"]][1, "entropy"] > 0.2) {
     status <- "high entropy"
   }

   if (mode.res[["mode.tab"]][1, "frac.het"] > 0.2) {
     status <- "high non-clonal"
   }

   return(status)
}
