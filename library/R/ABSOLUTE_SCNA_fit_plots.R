GetModeColors <- function() {
  return(c("springgreen3", "dodgerblue", "maroon1", "goldenrod1",
           "orangered2", 4, "salmon", "yellowgreen", "cyan",
           "sienna", "rosybrown", rep(1,100) ))
}

PlotModes_layout = function()
{
  NPLOT=6
  layout( mat=matrix( 1:(5*NPLOT), nrow=5, ncol=NPLOT, byrow=TRUE), widths = c(1, 2.5, rep(1, NPLOT-2)), heights=1 )

  par(las=1)
  par( tcl = -0.15 )
  par( mgp=c(1.5, 0.3, 0))
#  par( mar=c(2.5, 4.6, 2, 1))
#  par( cex=0.6 )
}

PlotModes <- function(segobj, chr.arms.dat, n.print = NA, called.mode.ix=NA, verbose=FALSE, plot.total.CN=FALSE)
{
#  Q = dim(segobj[["mode.res"]][["seg.q.tab"]])[3]
  Q = segobj[["mode.res"]][["mode_SCNA_models"]][[1]][["kQ"]] 

  alpha.dom <- c(0, 1)
  tau.dom <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["kTauDom"]]
  mode.colors <- GetModeColors()
  max_CR <- 2.25
  binW <- 0.025
  N_PLOTS=6  ##  Per row

  mode.tab <- segobj[["mode.res"]][["mode.tab"]]

  obs <- ExtractSampleObs(segobj)
  allele.segs = get_hom_pairs_segtab( segobj )

  SN <- segobj[["sample.name"]]
  model.id <- segobj[["group"]]

  if(is.na(n.print)) { n.print = nrow(mode.tab) }
  n.plot <- min(n.print, NROW(mode.tab))
  mode.tab <- mode.tab[c(1:n.plot), , drop = FALSE]
    
  for (i in seq_len(n.plot))
  {
      SCNA_model = segobj[["mode.res"]][["mode_SCNA_models"]][[i]] 
      SSNV_model = segobj[["mode.res"]][["mode_SSNV_models"]][[i]] 

      res = get_b_and_delta(  mode.tab[i, "alpha"], mode.tab[i, "tau"] )
      delta = res$delta
      b = res$b
##
      tree_clust = resort_tree_clusters( SCNA_model, SCNA_model[["seg_CCF_DP"]][["tree_clust"]] )
      SCNA_model[["seg_CCF_DP"]][["tree_clust"]] = tree_clust
      SCNA_model[["seg_CCF_DP"]][["seg_clust_tab"]] = get_seg_clust_tab( SCNA_model )
##


## 1- alpha vs tau
      modes_purity_ploidy_plot(mode.tab, mode.colors, alpha.dom, tau.dom, SN, segobj, called.mode.ix,
                               call.status=segobj[["mode.res"]][["call.status"]], model.id=model.id, mode.focus.ix=i)


      res = get_b_and_delta(  mode.tab[i, "alpha"], mode.tab[i, "tau"] )
      delta = res$delta
      b = res$b
      mode.info <- mode.tab[i, ]
      comb <- GetCopyRatioComb(Q, delta, b, obs[["error.model"]])
      
      if (!is.null(SCNA_model[["ab.tab"]])) {
        comb.ab <- SCNA_model[["ab.tab"]]
      } else {
        comb.ab <- NA
      }
      
      unif=SCNA_model[["seg.Wu.tab"]]
      exp = rep(0, length(unif))
      clonal=1-(exp+unif)
      clonal[clonal<0]=0   # round-off error

#      clonal_seg_colors = rgb( "red"=exp, "green"=clonal, "blue"=unif )
      colpal = colorRampPalette( c("magenta", "deepskyblue"))(1000)
#      colpal = colorRampPalette( c("darkslateblue", "coral"))(1000)
      pal.idx = floor(clonal * 999) + 1
      clonal_seg_colors = colpal[pal.idx] 

      Wq0 = get_comb_Wq0( obs[["e.cr"]], comb, log(mode.info["theta.Q.hat"]), mode.info["theta.0"] )

      if( segobj[["copy_num_type"]] == "allelic" ) { copy_ratio_label="Allelic copy number" }
      if( segobj[["copy_num_type"]] == "total" ) { copy_ratio_label="Total copy ratio" }


    ## Original plot style with seg hist and comb fit
#      old.mar <- par("mar")
#      if (plot.model.fit) {
#         par(mar = old.mar + c(0, 0, 1, 0))  ## add margin line to top
#      }
#      plot_ABS_seg_hist(obs[["d.tx"]], obs[["W"]], copy_ratio_label, clonal_seg_colors, max_CR, sideways=sideways )
#      plot_ABS_comb_fit( comb, sideways, col=mode.color, max_CR=max_CR, plot.model.fit=TRUE, mode.info=mode.info, Wq0=Wq0, comb.ab=comb.ab )
#      par(mar = old.mar)

    # New plot version with genome-plot and sideways hist summary
# 2 and 3 - genome and seghist
      PlotHscrAndSeghist( allele.segs, total.seg.dat=allele.segs, clonal_seg_colors, chr.arms.dat, max_CR, plot.hist=TRUE, plot.abs.fit=TRUE, plot.total.CN=plot.total.CN, comb=comb, mode.info=mode.info, Wq0=Wq0, comb.ab=comb.ab, fit.color=mode.colors[i], plot.seg.sem=TRUE )
      title("Seg clonality", line = 0.5, cex.main = 0.9)

  # 4 - SCNAs CCF
#      SCNA_CCF_plot( segobj, mode.ix=i )

  # 5 - SCNAs DP predictive density
      barplot( SCNA_model[["seg_CCF_DP"]][["CCF_DP_dens_est"]][-1], space=0, border=FALSE, col=1 )
     ## add KDE for comparison
      gr = c(1:length(SCNA_model[["collapsed_CCF_KDE"]]))[-1]
      lines( gr, SCNA_model[["collapsed_CCF_KDE"]][-1], col="red", lwd=1, lty=2 )
      pred_95_CI = SCNA_model[["seg_CCF_DP"]][["predictive_density_summary"]][c(1,3),-1]
      pred_median = SCNA_model[["seg_CCF_DP"]][["predictive_density_summary"]][2,-1]

      lines( gr, pred_95_CI[1,], col="blue", lwd=1 )
      lines( gr, pred_95_CI[2,], col="blue", lwd=1 )
      lines( gr, pred_median, col="blue", lwd=1, lty=2 )

      title("DP predictive density", line = 0.5, cex.main = 0.9)

      nc = ncol(SCNA_model[["collapsed_CCF_dens"]])
      GRID = cumsum( c(0, rep(1/(nc-1),(nc-1))))

      DP.out.ix = apply( SCNA_model[["seg.ix.tab"]][ , c("amp.ix", "neg.ix", "clonal.ix", "high.sem.ix") ], 1, any )
      clonal.ix = SCNA_model[["seg.ix.tab"]][ , "clonal.ix"] 
      nix = apply( SCNA_model[["seg.ix.tab"]][ , c("amp.ix", "neg.ix", "high.sem.ix", "clonal.ix") ], 1, any )

      before_dens =  SCNA_model[["collapsed_CCF_dens"]][ !nix, ]
      after_dens =  SCNA_model[["seg_CCF_DP"]][["collapsed_DP_CCF_dens"]][ !nix,]
      SID = ""

      mut_cols = SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["assign"]]
     
      if(all(is.na(mut_cols))) 
      { 
         mut_cols = rep(1, length(mut_cols)) 
         clust_cols = 1
      }
      else 
      {
         clust_cols = c(1:max(mut_cols, na.rm=TRUE))
      }

    ## expected CCF value of each cluster
# does not work as expected >> ?
#      E_clust_CCF = apply( GRID * SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["CCF_dens"]], 1, sum )
#      pal.idx = floor( (1-E_clust_CCF) * 999) + 1
#      seg_clust_CCF_colors = colpal[pal.idx] 

#      seg_assign = SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["assign"]]
#      n.ix = is.na(seg_assign)
#      mut_cols = rep(NA, length(seg_assign) )
#      mut_cols[!n.ix] = seg_clust_CCF_colors[ seg_assign[!n.ix] ]

      if (length(c(before_dens)) > 1) {
        YMAX = max( c(before_dens) )
      } else {
        YMAX = 0
      }
      if( nrow(before_dens) > 0 ) 
      {
         sample_trans_ccf_plot( mut_cols[!nix], GRID, before_dens, SID, YMAX )
        title("Before clustering", line = 0.5, cex.main = 0.9)
      } else{ frame() }
      if( nrow(after_dens) > 0 )
      {
         sample_trans_ccf_plot( mut_cols[!nix], GRID, after_dens, SID, YMAX )
        title("After clustering", line = 0.5, cex.main = 0.9)
      } else{ frame() }

## new row
      PpModeScorerBarplot(mode.tab, mode.colors, obs, n.plot)
      # frame()
## color genome by seg clust
#      if( length(mut_cols) != nrow(allele.segs)) { stop("wrong # of cols/segs") }

      PlotHscrAndSeghist( allele.segs, total.seg.dat=allele.segs, mut_cols, chr.arms.dat, max_CR, plot.hist=TRUE, plot.abs.fit=TRUE, plot.total.CN=plot.total.CN, comb=comb, mode.info=NA, Wq0=Wq0, comb.ab=comb.ab, fit.color=mode.colors[i], plot.seg.sem=TRUE )
      title("Seg clusters", line = 0.5, cex.main = 0.9)

      seg_LL = SCNA_model[["seg_LL"]] 
      # seg_CN_LL = SCNA_model[["seg_CN_LL"]]
      use.pal = (colorRampPalette(c("red", "yellow", "green")))  (1000)
# color by seg LL
      seg_LL_cols = get_seg_colors(seg_LL, use.pal=use.pal)
      PlotHscrAndSeghist(allele.segs, total.seg.dat=allele.segs, seg_LL_cols, chr.arms.dat, max_CR=5.0, plot.genome=FALSE, plot.hist=TRUE, plot.abs.fit=TRUE, plot.total.CN=plot.total.CN, comb=comb, mode.info=NA, Wq0=Wq0, comb.ab=comb.ab, fit.color=mode.colors[i], plot.seg.sem=TRUE, )
      title("Seg LL", line = 0.5, cex.main = 0.9)

#  CCF of hard-clusters
      clust_dens = SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["CCF_dens"]] 
      YMAX = max( apply( clust_dens, 1, max, na.rm=TRUE), na.rm=TRUE )
      sample_trans_ccf_plot( clust_cols, GRID, clust_dens, "", YMAX )

   if( !is.null(SCNA_model[["DP_res"]]) )
   {
## Plot of DP K prior and sampling-distr
#      N_iter = length(SCNA_model[["DP_res"]][["DP_post"]][["K"]]) - 1
#      N_burn= floor(N_iter/2)
      N_burn= SCNA_model[["N_DP_burn"]]   
      DP_post_K_plot( SCNA_model[["DP_res"]], N_burn )
   }
   else { frame() }


# new row:
   frame()
 ## plot genome colored according to seg.ix.tab
    seg_ix_colors = rep( NA, nrow(SCNA_model[["seg.ix.tab"]]) )
    no.ix = !apply( SCNA_model[["seg.ix.tab"]], 1, any )
    seg_ix_colors[ no.ix ] = 1     ## No flags
    seg_ix_colors[ !no.ix ] = apply( SCNA_model[["seg.ix.tab"]][ !no.ix , c("amp.ix", "neg.ix", "high.sem.ix"), drop=FALSE ], 1, which.max ) + 1
    seg_ix_colors[ SCNA_model[["seg.ix.tab"]][, 5] ] = max(seg_ix_colors, na.rm=TRUE)+1
    PlotHscrAndSeghist( allele.segs, total.seg.dat=allele.segs, seg_ix_colors, chr.arms.dat, max_CR, plot.hist=TRUE, plot.abs.fit=TRUE, plot.total.CN=plot.total.CN, comb=comb, plot.seg.sem=TRUE )

   ## SCNA ev summary
    if (plot.total.CN) {
      ## total CR: show the chr-arm minimum-event score per WGD hypothesis (lower = preferred).
      arm_ev_result = SCNA_model[["SCNA_minev_chrarm_result"]]
      e_wgd = arm_ev_result[["e_wgd"]]
      if (!is.null(e_wgd) && all(is.finite(e_wgd))) {
        bcols = rep("grey70", length(e_wgd)); bcols[which.min(e_wgd)] = "firebrick"
        barplot( e_wgd, names=c("WGD0","WGD1","WGD2")[seq_along(e_wgd)],
                 ylab="chr-arm events (lower=preferred)", col=bcols )
        title(paste0("selected WGD", arm_ev_result[["WGD"]]), line=0.5, cex.main=0.9)
      } else { frame() }
      frame()
    } else {
      arm_ev_result = SCNA_model[["SCNA_minev_chrarm_result"]]
      barplot( arm_ev_result[["WGD_Pr"]], names=c("WGD0", "WGD1", "WGD2"), ylab="Probability", ylim=c(0,1) )

  #    barplot( c(arm_ev_result[["ancestral_LL"]], arm_ev_result[["subclonal_LL"]],  arm_ev_result[["num_NA"]]), names=c("ancestral CN", "subclonal CN", "num_NA"), ylab="Log likelihood"  )
      ev_plot_dat = arm_ev_result[ c("average_WGD_score", "log_multinom_coef", "SC_ev_steps_LL", "neg_arm_LL", "subclone_CCF_volume") ]
      barplot( unlist(ev_plot_dat), names=names(ev_plot_dat), las=2, ylab="Log likelihood" )
    }

#"WGD_Pr"=WGD_Pr, "ancestral_LL"=WGD_ev_result, "subclonal_LL"=log_multinom_coef, "num_NA"=num_NA_chrarms
    frame()

# new row
    frame()
 ## plot genome with total segs 
#    tcols = rep( 1, nrow(segobj[["total.seg.dat"]]) )
#    PlotHscrAndSeghist( allele.segs, seg_LL_cols, chr.arms.dat, max_CR=5, plot.hist=TRUE, plot.abs.fit=TRUE, comb=comb, plot.seg.sem=TRUE,
# plot.total.CN=TRUE, tot.seg.colors=tcols ) 

## plot d0_segs
    if (!plot.total.CN) {
      segs_d0 = deconstruct_SCNAs( SCNA_model, obs, allele.segs, b, delta )
      AS.seg.ix = allele.segs[, c("seg.ix.1", "seg.ix.2")]
      d0.allele.segs = allele.segs

      d0.allele.segs[,"A1.Seg.CN"] = segs_d0[ AS.seg.ix[,1] ]
      d0.allele.segs[,"A2.Seg.CN"] = segs_d0[ AS.seg.ix[,2] ]

      PlotHscrAndSeghist( d0.allele.segs, total.seg.dat=allele.segs, mut_cols, chr.arms.dat, max_CR=4.25, plot.hist=TRUE, plot.abs.fit=FALSE, plot.total.CN=plot.total.CN, comb=comb, plot.seg.sem=FALSE, y.lab=copy_ratio_label )

      frame(); frame(); frame()

  # new row
      frame()
    ## Genome plot
      PlotHscrAndSeghist( d0.allele.segs, total.seg.dat=allele.segs, mut_cols, chr.arms.dat, max_CR=2.5, plot.abs.fit=FALSE, plot.total.CN=plot.total.CN, comb=comb, plot.seg.sem=FALSE, y.lab=copy_ratio_label )

      if(!is.null(segobj[["mode.res"]][["modeled.muts"]][[i]]) && !all(is.na(segobj[["mode.res"]][["modeled.muts"]][[i]][,"ccf_hat"])) )  ## protect against edge case of all muts on homozygously del SCNAs
      {
         mut.cn.dat <- segobj[["mut.cn.dat"]]
         modeled <- segobj[["mode.res"]][["modeled.muts"]][[i]]
         ## carry the per-mode ploidy (tau) onto the mutation table so the density plots can
         ## reference it (e.g. the alpha/tau VAF line in Mut_AF_plot); mode.tab is the only place
         ## ploidy lives -- modeled.muts itself only carries purity + SSNV_skew.
         modeled.mut.dat <- cbind(mut.cn.dat, modeled, tau = mode.tab[i, "tau"])

         max_SSNVs_plot = 500

       ## plot SSNVs on genome
         # SSNV_cols = c("dodgerblue", "darkgrey", "seagreen3")   ## SC, clonal, mult>1
         SSNV_cols = c("magenta", "olivedrab", "dodgerblue")   ## SC, clonal, mult>1
         plot_SSNVs_on_genome( SSNV_model, SSNV_cols, modeled.mut.dat, segobj, i, mode.colors[i], chr.arms.dat, max_SSNVs_plot=max_SSNVs_plot, verbose=verbose)

    ## Now plot SSNVs densities ... 2 plots
          PlotSomaticMutDensities(modeled.mut.dat, segobj, i,
                                  mode.colors[i], min.cov=3, max_SSNVs_plot=max_SSNVs_plot, verbose=verbose)
  #  2 plots
          plot_SSNV_on_subclonal_SCNA_CCF_summaries(modeled.mut.dat, segobj, i, min.cov=3, max_SSNVs_plot=max_SSNVs_plot, verbose=verbose)
      }
      else { for(i in 1:4){ frame() } }

    } else {
      ## total CR: no allele-specific SCNA deconstruction, but draw the SSNV density panels
      ## (point-mutation VAF, multiplicity, and clonal-SSNV CCF) from the total-CN SSNV fit.
      ## The "SSNVs on subclonal SCNAs" panel stays blank (no allele-resolved subclonal model).
      ## Cell budget for this branch is 11 (the 5x6 mode grid); keep it balanced with frames.
      if( !is.null(segobj[["mode.res"]][["modeled.muts"]][[i]]) && !all(is.na(segobj[["mode.res"]][["modeled.muts"]][[i]][,"ccf_hat"])) )
      {
         mut.cn.dat <- segobj[["mut.cn.dat"]]
         modeled <- segobj[["mode.res"]][["modeled.muts"]][[i]]
         ## carry the per-mode ploidy (tau) onto the mutation table so the density plots can
         ## reference it (e.g. the alpha/tau VAF line in Mut_AF_plot); mode.tab is the only place
         ## ploidy lives -- modeled.muts itself only carries purity + SSNV_skew.
         modeled.mut.dat <- cbind(mut.cn.dat, modeled, tau = mode.tab[i, "tau"])
         max_SSNVs_plot = 500

         for(k in 1:7){ frame() }
       ## 2 plots: point-mutation VAF + multiplicity
         PlotSomaticMutDensities(modeled.mut.dat, segobj, i, mode.colors[i], min.cov=3, max_SSNVs_plot=max_SSNVs_plot, verbose=verbose)
       ## 2 plots: clonal-SSNV CCF + (blank) SSNVs-on-subclonal-SCNA
         plot_SSNV_on_subclonal_SCNA_CCF_summaries(modeled.mut.dat, segobj, i, min.cov=3, max_SSNVs_plot=max_SSNVs_plot, verbose=verbose)
      }
      else { for(k in 1:11){ frame() } }
    }
  }  ## modes

}


modes_purity_ploidy_plot <- function(mode.tab, mode.colors, alpha.dom, tau.dom, sample.name, 
                      seg.dat, called.mode.ix, debug.info=TRUE, call.status = NA, model.id = NA, mode.focus.ix=1)
{
  mode.colors <- mode.colors[c(1:NROW(mode.tab))]
    
  a1 <- mode.tab[mode.focus.ix, "alpha"]
  t1 <- mode.tab[mode.focus.ix, "tau"]
  den <- 2 * (1 - a1) + a1 * t1
  delta <- a1 / den
  b <- 2 * (1 - a1) / den
    
  delta.set <- c(delta / 2, delta, 2 * delta)
  b.set <- c(b, b - 2 * delta, b + 2 * delta)
  b.set <- b.set[b.set >= 0 & b.set < 1]
    
  if( !is.na(mode.focus.ix) ) {
     pp.grid.col = mode.colors[mode.focus.ix]
  } else { pp.grid.col = "grey30" }

  old.mar <- par("mar")
  if (debug.info) {
    ## add margin line to top
    par(mar = old.mar + c(0, 0, 1, 0))  
  }
  
  ylab <- substitute( paste( "Fraction cancer nuclei", (hat(alpha)), sep="") )
  xlab <- substitute( paste( "Ploidy", (hat(tau)), sep="") )

  plot(0, type = "n", ylab = ylab, xlab = xlab,
       xlim = tau.dom, ylim = alpha.dom, main = "", bty = "n", las = 1,
       cex = 1.5, xaxt = "n")
  
  at <- seq(2, max(tau.dom), by = 2)
  axis(side = 1, at = at, labels = paste(at, "N", sep = ""))
  
  for (i in seq_along(b.set)) {
    curve(AlphaBFunc(x, b.set[i]), from = tau.dom[1], to = tau.dom[2],
          col = pp.grid.col, add = TRUE, lty = 3)
  }
  
  for (i in seq_along(delta.set)) {
    curve(AlphaFunc(x, delta.set[i]), from = tau.dom[1], to = tau.dom[2],
          col = pp.grid.col, add = TRUE, lty = 3, n = 500)
  }
  
  points(mode.tab[, "tau"], mode.tab[, "alpha"], col = "black", bg = mode.colors, 
         pch = 21, cex = par("cex") * 2.0, lwd = 0.5)
  
  if (!is.na(called.mode.ix)) {
    ix <- called.mode.ix
    points(mode.tab[ix, "tau"], mode.tab[ix, "alpha"], col = "red", pch = "*", 
           cex = par("cex") * 4)
  }
  
  if (debug.info) {
    title(paste("Solution", mode.focus.ix), line = 3, col.main = pp.grid.col, cex.main = 2)
    # mtext(call.status, line = 0, adj = 0)
    mtext(sample.name, side=3, line=1, adj=0, cex=par("cex") )
    mtext(paste("Disease:", model.id), side=3, line=0, adj=0, cex=par("cex"))
  }
  
  if ((seg.dat[["platform"]] == "SNP_6.0") && debug.info) {
    hscn.params <- seg.dat[["error.model"]]
    mtext(substitute(paste(sigma[eta] == x, ",", sigma[epsilon] == y, sep = ""), 
                     list(x = round(hscn.params[["sigma.eta"]], 3),
                          y = round(hscn.params[["sigma.nu"]], 3))),
          line = 1.25, adj = 0)
  }
  
  par(mar = old.mar)
}


PpModeScorerBarplot <- function(mode.tab, mode.colors, obs, n.print) {
  n.print <- min(n.print, nrow(mode.tab))
  mode.colors <- mode.colors[c(1:NROW(mode.tab))]
  
  use_cols = c()
  use_names = c()

  if (("SCNA_LL" %in% colnames(mode.tab)) &&
      (!any(is.na(mode.tab[1:n.print, "SCNA_LL"])))) {
    use_cols = c(use_cols, "SCNA_LL")
    use_names = c(use_names, "SCNAs")
  }

  if(( "SCNA_min_chrarm_events" %in% colnames(mode.tab)) &&
      (!any(is.na(mode.tab[1:n.print, "SCNA_LL"])))) {
    use_cols = c(use_cols, "SCNA_min_chrarm_events")
    use_names = c(use_names, "CN Ev")
  }


#  if (("Kar_LL" %in% colnames(mode.tab)) &&
#      (!any(is.na(mode.tab[1:n.print, "Kar_LL"])))) {
#    use_cols = c(use_cols, "Kar_LL")
#    use_names = c(use_names, "karyotype")
#  }

  if (("SSNV_LL" %in% colnames(mode.tab)) &&
      (!any(is.na(mode.tab[1:n.print, "SSNV_LL"])))) {
    use_cols = c(use_cols, "SSNV_LL")
    use_names = c(use_names, "SSNVs")
  }

   mat <- mode.tab[1:n.print, use_cols, drop = FALSE]
  colnames(mat)=use_names

  
  ix <- apply(!is.finite(mat), 1, sum) > 0
  mat <- mat[!ix, , drop=FALSE]
  
#  mat <- mat - min(mat)
  if( nrow(mat) > 1 )
  {
    mat <- t( t(mat) - apply(mat, 2, min) )
    mat <- mat + max(mat) * 0.1
  }
  mat <- cbind(mat, rowMeans(mat))
  colnames(mat)[ncol(mat)] <- "combined"
  
#  cols = mode.colors[(1:n.print)[!ix]]
  cols = mode.colors[ !ix ]
  barplot(mat, beside = TRUE, col = cols, axes = FALSE, las = 2,
          ylab = "", space = c(0, 2), cex.names = par("cex.axis"))
  
  mtext("Log-likelihood", side = 2, line = 1, las = 3, cex = par("cex") * par("cex.axis"))
#  mtext("Model-based evaluation", side = 3, line = 0, adj = 0, cex = par("cex.axis"))
  
  axis(side=2, labels=FALSE)
}



plot_ABS_seg_hist = function(d, W, copy_ratio_label, seg_colors, min_CR=0, max_CR=Inf, bin.w=0.025, sideways=FALSE, label.genomic.fraction.axis=TRUE, genomic_fraction_label="Genomic fraction" )
{
  ix <- (d >= min_CR) & (d < max_CR)

  PlotSeglenHist(d[ix], W[ix], seg_colors[ix], 
                 bin.w = bin.w, x.max = max_CR, x.min=min_CR,
                 data.whiskers = FALSE, xlab = copy_ratio_label, ylab = genomic_fraction_label, xlim=c(min_CR, max_CR), x.axis=FALSE,
                 sideways = sideways, y.ax.labs=label.genomic.fraction.axis )
}


## Comb-level tick label. In total CR mode the axis is the total copy ratio, so a comb line at
## absolute integer copy number `cn` sits at copy ratio cn/2 (germline diploid CN2 = CR 1; this
## also makes the single genome-wide comb read correctly over male X/Y, where germline CN1 = CR 1).
## Integers are shown plainly; half-integers as a compact nicefrac-style super/subscript fraction
## (base-R plotmath has no \nicefrac, so emulate it as {}^num "/" {}[den]). Allelic mode is
## unchanged: it keeps the raw integer copy number.
comb_level_label = function( cn, cr.labels )
{
  if (!cr.labels) { return(cn) }
  if (cn %% 2 == 0) { return(cn %/% 2) }
  return( bquote( "  " * {}^.(cn) * "/" * {}[2] ) )
}

plot_ABS_comb_fit = function( sideways, col, max_CR, comb, mode.info, Wq0, comb.ab, plot.model.fit=FALSE, cr.labels=FALSE )
{
  if (plot.model.fit)
  {
     msg = get_SCNA_model_fit_plot_strings( mode.info )

     for( i in 1:length(msg)) { 
        mtext(msg[[i]], line = 3-i+1, cex = 0.75, adj = 0)
     }
  }

#  extra <- (comb[2] - comb[1]) / 2
  for (q in seq_along(comb)) 
  {
    if (comb[q] > max_CR) {
      break
    }
  }
  max.q <- q 
  
  for (q in seq_along(comb))
  {
    if (q > max.q) {
      next
    }
    
    ## ABS CN at each level
    if (!sideways) 
    {
      abline(v = comb[q], lwd = 0.5, lty = 3, col = col)
      lines( x=c(comb[q], comb[q]), y=c(0,Wq0[q]), lwd=0.5, lty=1, col="black")
    }
    else
    {
      abline(h = comb[q], lwd = 0.5, lty = 3, col = col)
      lines( y=c(comb[q], comb[q]), x=c(0,Wq0[q]), lwd=0.5, lty=1, col="black")
    }
    
    side <- ifelse(!sideways, 3, 4)
    
    if ((q - 1 < 10 | (q - 1)%%2 == 1) & comb[q] < max_CR ) {
      mtext(text = comb_level_label(q - 1, cr.labels), side = side, at = comb[q], col = col, line = 0.2, cex = par("cex") * par("cex.axis"))
    }
    
    ## % AB in each level
    if (!all(is.na(comb.ab))) {
      q.ab <- round(comb.ab[q], 2)
      if (is.finite(q.ab) & q.ab > 0) 
      {
         mtext(text = paste("(", q.ab, ")", sep = ""), side=side, at = comb[q], col = 1, line = 1, cex = 0.5)
      }
    }
  }
}


  
AlphaBFunc <- function(tau, b) {
  alpha <- 2 * (1 - b) / (tau * b - 2 * b + 2)
  return(alpha)
}


AlphaFunc <- function(tau, delta) {
  alpha <- 2 * delta / (1 - delta * tau + 2 * delta)
  alpha[alpha > 1] <- NA
  return(alpha)
}


SCNA_CCF_plot = function( seg.dat, mode.ix )
{
 ## Plot SCNA rescaling to cancer cell-fraction
#  SCNA_cols = c("cyan4", "chocolate4") 
  SCNA_cols = c("blue", "green") 

  sc_tab = seg.dat$mode.res$subclonal_SCNA_res$subclonal_SCNA_tab[mode.ix, , ]
  SCNA_pr_subclonal = sc_tab[, "Pr_subclonal"]
  SCNA_pr_clonal = 1 - SCNA_pr_subclonal
  SCNA_pr_clonal[SCNA_pr_clonal < 0] = 0 ## round-off error
  SCNA_ccf_dens = exp( seg.dat$mode.res$subclonal_SCNA_res$log_CCF_dens[mode.ix, , ] )
  ccf_grid = as.numeric(colnames(SCNA_ccf_dens))

#  SCNA_ccf_dens = SCNA_ccf_dens[ SCNA_pr_clonal < 0.5, , drop=FALSE]
#  SCNA_pr_clonal = SCNA_pr_clonal[SCNA_pr_clonal<0.5]

#  SCNA_ymax = max(SCNA_ccf_dens, na.rm=TRUE)
  SCNA_ymax = max(  apply( SCNA_ccf_dens, 1, max, na.rm=TRUE), na.rm=TRUE )
  y_lim = SCNA_ymax
  scale_total = y_lim/2

  SCNA_mut_cols = get_mut_cols( SCNA_cols, SCNA_pr_clonal )

# SCNA  CCFs
   draw_grid_mut_densities( SCNA_ccf_dens, ccf_grid, col="red", mut_colors=SCNA_mut_cols, scale_total, draw_total=TRUE, draw_indv=TRUE, add=FALSE, y_lim=y_lim )


   legend( x='topleft', legend=c( "subclonal CN", "clonal CN"), col=SCNA_cols, lty=1, lwd=1.5, bty="n")
}


get_seg_colors = function(color.by, use.pal=NA, color.range=NA)
{
  if(any(is.na(use.pal))) { use.pal <- heat.colors(1000) }

  col.scale <- length(use.pal)
  
  if (is.na(color.range)) {
    color.range <- range(color.by, na.rm=TRUE)
  }

  if( color.range[1] == color.range[2] )
  {
    # Edge-case where all colors are the same
     mid.idx = floor( length(use.pal)/2 )
     cols = rep( use.pal[mid.idx], length(color.by))
  }
  else
  {
     pal.idx <- floor((color.by - color.range[1]) /
                      (color.range[2]-color.range[1]) * (col.scale-1))
     cols = use.pal[pal.idx]
  }
 
  return(cols)
}
