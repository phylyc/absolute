ApplySSNVModel <- function(mode.res, mut.cn.dat, SSNV_model, verbose=FALSE)
{
  if (verbose) {
    print(paste("Evaluating ", nrow(mut.cn.dat),
                " mutations over ", nrow(mode.res[["mode.tab"]]),
                " purity/ploidy modes: ", 
                sep = ""))
  }

  n_modes = nrow(mode.res[["mode.tab"]])
  ccf_grid = SSNV_model[["ccf_grid"]]
  mode.res[["SSNV.ccf.dens"]] = array( NA, dim=c(n_modes, nrow(mut.cn.dat), length(ccf_grid)) )
  dimnames(mode.res[["SSNV.ccf.dens"]])[[3]] = ccf_grid
  
  mode.res[["modeled.muts"]] = list()
  mode.res[["mode_SSNV_models"]] = list()

  for (j in 1:n_modes) 
  {
    alpha <- mode.res[["mode.tab"]][j, "alpha"]

    seg.q.tab = mode.res[["mode_SCNA_models"]][[j]] [["seg.q.tab"]]
    subclonal_scna_tab = mode.res[["subclonal_SCNA_res"]][["subclonal_SCNA_tab"]][j, , ]
    SCNA_log_ccf_dens = mode.res[["subclonal_SCNA_res"]][["log_CCF_dens"]][j, , ]

    # tot.seg.q.tab =  mode.res[["mode_SCNA_models"]][[j]] [["tot.seg.q.tab"]]
    # total_subclonal_scna_tab = mode.res[["subclonal_SCNA_res"]][["total_subclonal_scna_tab"]][j, , ]
    # total_scna_log_ccf_dens = mode.res[["subclonal_SCNA_res"]][["total_log_ccf_dens"]][j, , ]

    
    res = FitPPModeSomaticMuts(
      SSNV_model, mode.res[["mode_SCNA_models"]][[j]], mut.cn.dat, mode.res$mode.tab[j, ],
      subclonal_scna_tab, SCNA_log_ccf_dens, seg.q.tab,
      # total_subclonal_scna_tab, total_scna_log_ccf_dens, tot.seg.q.tab, verbose=verbose
    )
        
    res[["modeled.muts"]][["purity"]] = alpha
    res[["modeled.muts"]][["SSNV_skew"]] = SSNV_skew
    modeled.muts = res[["modeled.muts"]]

    mode.res[["modeled.muts"]][[j]] = modeled.muts

    som.theta.q.map <- res[["som.theta.q.map"]]
    mode.res[["SSNV.ccf.dens"]][j,,] <- res[["ccf.dens"]] 
#    rownames(modeled.muts)=NULL
#    mode.res[["modeled.muts"]][[j]] = data.frame(modeled.muts, "purity" = alpha, "SSNV_skew"=SSNV_model[["SSNV_skew"]], stringsAsFactors=FALSE, check.names=FALSE )

    mode.res[["mode_SSNV_models"]][[j]] = res[["mode_SSNV_models"]]
    
    mode.res[["mode.tab"]][j, "SSNV_LL"] <- sum(modeled.muts[, "LL"], na.rm = TRUE) + dDirichlet(som.theta.q.map, SSNV_model[["kPiSomThetaQ"]], log.p = TRUE)
 
    if (verbose) { cat(".") }
  }
  if (verbose) {
    cat("\n")
  }
  
#  new.ll <- mode.res[["mode.tab"]][, "combined_LL"] +
#            mode.res[["mode.tab"]][, "SSNV_LL"]
#  mode.res[["mode.tab"]][, "combined_LL"] <- new.ll

  return(mode.res)
}




FitPPModeSomaticMuts <- function(
  SSNV_model, SCNA_model, mut.cn.dat, mode_info,
  subclonal_scna_tab, scna_log_ccf_dens, seg.q.tab,
  total_subclonal_scna_tab, total_SCNA_log_ccf_dens, tot.seg.q.tab, verbose=FALSE)
{
  clonal_scna_mut_ix = !get_subclonal_scna_mut_ix(mut.cn.dat, subclonal_scna_tab)

  res = get_subclonal_scna_tab(mut.cn.dat[!clonal_scna_mut_ix,], subclonal_scna_tab, seg.q.tab )
  subclonal_scna_tab = res[["subclonal_scna_tab"]]
  subclonal.mut.tab = as.data.frame(matrix( NA, nrow=nrow(mut.cn.dat), ncol=ncol(res[["subclonal.mut.tab"]]) ))
  colnames(subclonal.mut.tab) = colnames(res[["subclonal.mut.tab"]])
  subclonal.mut.tab[!clonal_scna_mut_ix,] = res[["subclonal.mut.tab"]]



  N_SSNV = nrow(mut.cn.dat)

  total.clonal.mut.tab = total_get_muts_nearest_clonal_scna(mut.cn.dat, SCNA_model[["seg.q.tab"]], SCNA_model[["kQ"]])
  colnames(total.clonal.mut.tab) = c("total_q_hat")
  allelic.clonal.mut.tab = matrix(NA, nrow=N_SSNV, ncol=3)
  colnames(allelic.clonal.mut.tab) = c("q_hat", "HS_q_hat_1", "HS_q_hat_2")

  if (("A1.ix" %in% colnames(mut.cn.dat)) && ("A2.ix" %in% colnames(mut.cn.dat))) {
    missing.AS.seg.ix = apply( is.na(mut.cn.dat[,c("A1.ix", "A2.ix")]), 1, any )
    if( any(!missing.AS.seg.ix))
    {
       allelic.clonal.mut.tab[!missing.AS.seg.ix,] = allelic_get_muts_nearest_clonal_scna(mut.cn.dat[!missing.AS.seg.ix,, drop=FALSE], SCNA_model[["seg.q.tab"]], SCNA_model[["kQ"]])
    }
  }

#  clonal.mut.tab = get_muts_nearest_clonal_scna(mut.cn.dat, seg.q.tab, SSNV_model[["kQ"]])

  mut.modeled.cn = cbind(allelic.clonal.mut.tab, total.clonal.mut.tab, subclonal.mut.tab, "clonal_scna_mut_ix"=clonal_scna_mut_ix)

  ## Total-CR muts have no allele-specific modal CN, so allelic.clonal.mut.tab leaves q_hat NA.
  ## Use the total-CN modal call (total_q_hat) as q_hat so the modeled.muts carried into the
  ## ABS MAF and the SSNV plot panels (VAF / multiplicity / CCF) treat q_hat uniformly across
  ## modes. Gated on the absence of allele-specific seg indices => allelic runs are unchanged.
  if (!all(c("A1.ix", "A2.ix") %in% colnames(mut.cn.dat))) {
    na.q = is.na(mut.modeled.cn[, "q_hat"])
    mut.modeled.cn[na.q, "q_hat"] = mut.modeled.cn[na.q, "total_q_hat"]
  }

  fit_res = fit_SSNV_model( cbind(mut.cn.dat, mut.modeled.cn), mode_info, SSNV_model, subclonal_scna_tab, scna_log_ccf_dens )
  som_theta_q_map = fit_res$som_theta_Q_MAP
  post_prs = fit_res$post_Prs
  ssnv.ccf.dens = fit_res[["ssnv.ccf.dens"]]
     
  ## Subclonal SCNA?
  var_classes = ClassifySomaticVariants(post_prs, 0.5)
  mut_mult_res = calc_CCF_95CI( cbind(mut.cn.dat, mut.modeled.cn), ssnv.ccf.dens, SSNV_model )

  detection_power = mode_SSNV_pow_calc( SSNV_model, mut.cn.dat, mut.modeled.cn, mode_info["alpha"] )
  detection_power_for_single_read = mode_SSNV_pow_calc( SSNV_model, mut.cn.dat, mut.modeled.cn, mode_info["alpha"], single_read=TRUE )

  modeled.muts <- cbind( mut.modeled.cn, post_prs, var_classes, mut_mult_res, fit_res[["H.ev"]], detection_power=detection_power, detection_power_for_single_read=detection_power_for_single_read )
#  modeled.muts <- cbind( mut.modeled.cn, post_prs, fit_res[["H.ev"]] )


  return(list( "ccf.dens" = ssnv.ccf.dens, "modeled.muts" = modeled.muts, "som.theta.q.map" = som_theta_q_map, "mode_SSNV_models"=fit_res[["SSNV_model"]]))
}



 
