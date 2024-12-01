PlotModes_review_layout = function()
{
  NCOL=5
  NROW=6
  layout( mat=matrix( 1:(NROW*NCOL), nrow=NROW, ncol=NCOL, byrow=TRUE), widths = c(1, 2.5, rep(1, NCOL-2)), heights=1 )

  par(las=1)
  par( tcl = -0.15 )
  par( mgp=c(1.5, 0.3, 0))
  par( mar=c(2.5, 4.6, 2, 1))
  par( cex=0.6 )
}

PlotMode_review_summary = function(segobj, chr.arms.dat, n.print = NA, called.mode.ix=NA, max_SSNVs_plot=500, verbose=FALSE)
{
  Q = segobj[["mode.res"]][["mode_SCNA_models"]][[1]][["kQ"]] 

  alpha.dom <- c(0, 1)
  tau.dom <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["kTauDom"]]
  mode.colors <- GetModeColors()
  max_CR <- 2.25
  binW <- 0.025

  mode.tab <- segobj[["mode.res"]][["mode.tab"]]

  obs <- ExtractSampleObs(segobj)
  allele.segs = get_hom_pairs_segtab( segobj )

  SN <- segobj[["sample.name"]]
  model.id <- segobj[["group"]]

  if(is.na(n.print)) { n.print = nrow(mode.tab) }
  n.plot <- min(n.print, NROW(mode.tab))
  mode.tab <- mode.tab[c(1:n.plot), , drop = FALSE]


## special 1st row with just mode-ranking barplot
  frame()
  PpModeScorerBarplot(mode.tab, mode.colors, obs, n.plot)
  mtext( SN, side=3, adj=0, cex=par("cex") )

  frame(); frame(); frame()
    
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
      segs_d0 = deconstruct_SCNAs( SCNA_model, obs, allele.segs, b, delta )
##

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

      colpal = colorRampPalette( c("magenta", "deepskyblue"))(1000)
#      colpal = colorRampPalette( c("darkslateblue", "coral"))(1000)
      pal.idx = floor(clonal * 999) + 1
      clonal_seg_colors = colpal[pal.idx] 

      Wq0 = get_comb_Wq0( obs[["e.cr"]], comb, log(mode.info["theta.Q.hat"]), mode.info["theta.0"] )

      if( segobj[["copy_num_type"]] == "allelic" ) { copy_ratio_label="Allelic copy ratio" }
      if( segobj[["copy_num_type"]] == "total" ) { copy_ratio_label="Total copy ratio" }


## 1- alpha vs tau
      modes_purity_ploidy_plot(mode.tab, mode.colors, alpha.dom, tau.dom, SN, segobj, called.mode.ix,
                               call.status=segobj[["mode.res"]][["call.status"]], model.id=model.id, mode.focus.ix=i)

# New plot version with genome-plot and sideways hist summary
# 2 and 3 - genome and seghist
## color genome by seg clust
 ## expected CCF value of each cluster
      nc = ncol(SCNA_model[["collapsed_CCF_dens"]])
      GRID = cumsum( c(0, rep(1/(nc-1),(nc-1))))

      E_clust_CCF = apply( GRID * SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["CCF_dens"]], 1, sum )
      pal.idx = floor( (1-E_clust_CCF) * 999) + 1
      seg_clust_CCF_colors = colpal[pal.idx] 

      seg_assign = SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["assign"]]
      n.ix = is.na(seg_assign)
      mut_cols = rep(NA, length(seg_assign) )
      mut_cols[!n.ix] = seg_clust_CCF_colors[ seg_assign[!n.ix] ]
      if( length(mut_cols) != 2*nrow(allele.segs)) { stop("wrong # of cols/segs") }

      PlotHscrAndSeghist( allele.segs, mut_cols, chr.arms.dat, max_CR, plot.hist=TRUE, plot.abs.fit=TRUE, comb=comb, mode.info=mode.info, Wq0=Wq0, comb.ab=comb.ab, fit.color=mode.colors[i], plot.seg.sem=TRUE )


## 2 plots for before / after seg CCF DP
      nix = apply( SCNA_model[["seg.ix.tab"]][ , c("amp.ix", "neg.ix", "high.sem.ix", "clonal.ix") ], 1, any )
      before_dens =  SCNA_model[["collapsed_CCF_dens"]][ !nix, ]
      after_dens =  SCNA_model[["seg_CCF_DP"]][["collapsed_DP_CCF_dens"]][ !nix,]
      SID = ""

      YMAX = max( c(before_dens) )
      if( nrow(before_dens) > 0 ) 
      {
         sample_trans_ccf_plot( mut_cols[!nix], GRID, before_dens, SID, YMAX )
      } else{ frame() }
      if( nrow(after_dens) > 0 )
      {
         sample_trans_ccf_plot( mut_cols[!nix], GRID, after_dens, SID, YMAX )
      } else{ frame() }
##


## new row
   frame()

  ## Genome plot at decon CN w. SSNV 95CIs overlayed
    AS.seg.ix = allele.segs[, c("seg.ix.1", "seg.ix.2")]
    d0.allele.segs = allele.segs
    d0.allele.segs[,"A1.Seg.CN"] = segs_d0[ AS.seg.ix[,1] ]
    d0.allele.segs[,"A2.Seg.CN"] = segs_d0[ AS.seg.ix[,2] ]

## can't overlay SSNVs if you plot.hist=TRUE :(
    PlotHscrAndSeghist( d0.allele.segs, mut_cols, chr.arms.dat, max_CR=2.5, plot.hist=FALSE, plot.abs.fit=FALSE, comb=comb, plot.seg.sem=FALSE, y.lab="Copy number" )
    frame()

    if(!is.null(segobj[["mut.cn.dat"]]) & !all(is.na(segobj[["mode.res"]][["modeled.muts"]][[i]][,"ccf_hat"])) )  ## protect against edge case of all muts on homozygously del SCNAs
    {
       mut.cn.dat <- segobj[["mut.cn.dat"]]
       modeled <- segobj[["mode.res"]][["modeled.muts"]][[i]]
       modeled.mut.dat <- cbind(mut.cn.dat, modeled)

     ## plot SSNVs on genome
       SSNV_cols = c("dodgerblue", "darkgrey", "seagreen3")   ## SC, clonal, mult>1
       plot_SSNVs_on_genome( SSNV_model, SSNV_cols, modeled.mut.dat, segobj, i, mode.colors[i], max_SSNVs_plot=max_SSNVs_plot, verbose=verbose)

## plot
  ## Now plot SSNVs densities ... 2 plots 
        PlotSomaticMutDensities(modeled.mut.dat, segobj, i,
                                mode.colors[i], min.cov=3, max_SSNVs_plot=max_SSNVs_plot, verbose=verbose)
    }
    else { for(i in 1:2){ frame() } }

  }  ## modes
}












# 1 solution per row
dens_PlotMode_review_summary = function(segobj, chr.arms.dat, n.print = NA, called.mode.ix=NA, max_SSNVs_plot=500, verbose=FALSE)
{
  Q = segobj[["mode.res"]][["mode_SCNA_models"]][[1]][["kQ"]] 
  alpha.dom <- c(0, 1)
  tau.dom <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["kTauDom"]]
  mode.colors <- GetModeColors()
  max_CR <- 2.25
  binW <- 0.025

  mode.tab <- segobj[["mode.res"]][["mode.tab"]]

  obs <- ExtractSampleObs(segobj)
  allele.segs = get_hom_pairs_segtab( segobj )

  SN <- segobj[["sample.name"]]
  model.id <- segobj[["group"]]

  if(is.na(n.print)) { n.print = nrow(mode.tab) }
  n.plot <- min(n.print, NROW(mode.tab))
  mode.tab <- mode.tab[c(1:n.plot), , drop = FALSE]


## special 1st row with just mode-ranking barplot
  frame()
  PpModeScorerBarplot(mode.tab, mode.colors, obs, n.plot)
  mtext( SN, side=3, adj=0, cex=par("cex") )

  frame(); frame(); frame()
    
  for (i in seq_len(n.plot))
  {
      SCNA_model = segobj[["mode.res"]][["mode_SCNA_models"]][[i]] 
      SSNV_model = segobj[["mode.res"]][["mode_SSNV_models"]][[i]] 

      res = get_b_and_delta(  mode.tab[i, "alpha"], mode.tab[i, "tau"] )
      delta = res$delta
      b = res$b

      tree_clust = resort_tree_clusters( SCNA_model, SCNA_model[["seg_CCF_DP"]][["tree_clust"]] )
      SCNA_model[["seg_CCF_DP"]][["tree_clust"]] = tree_clust
      SCNA_model[["seg_CCF_DP"]][["seg_clust_tab"]] = get_seg_clust_tab( SCNA_model )
      segs_d0 = deconstruct_SCNAs( SCNA_model, obs, allele.segs, b, delta )

      mode.info <- mode.tab[i, ]
      comb <- GetCopyRatioComb(Q, delta, b, obs[["error.model"]])
      
      if (!is.null(SCNA_model[["ab.tab"]])) {
        comb.ab <- SCNA_model[["ab.tab"]]
      } else {
        comb.ab <- NA
      }
      
      Wq0 = get_comb_Wq0( obs[["e.cr"]], comb, log(mode.info["theta.Q.hat"]), mode.info["theta.0"] )

      if( segobj[["copy_num_type"]] == "allelic" ) { copy_ratio_label="Allelic copy ratio" }
      if( segobj[["copy_num_type"]] == "total" ) { copy_ratio_label="Total copy ratio" }


## 1- alpha vs tau
      modes_purity_ploidy_plot(mode.tab, mode.colors, alpha.dom, tau.dom, SN, segobj, called.mode.ix,
                               call.status=segobj[["mode.res"]][["call.status"]], model.id=model.id, mode.focus.ix=i)

# New plot version with genome-plot and sideways hist summary
# 2 and 3 - genome and seghist
## color genome by seg clust
 ## expected CCF value of each cluster
      nc = ncol(SCNA_model[["collapsed_CCF_dens"]])
      GRID = cumsum( c(0, rep(1/(nc-1),(nc-1))))

      colpal = colorRampPalette( c("magenta", "deepskyblue"))(1000)
      E_clust_CCF = apply( GRID * SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["CCF_dens"]], 1, sum )
      pal.idx = floor( (1-E_clust_CCF) * 999) + 1
      seg_clust_CCF_colors = colpal[pal.idx] 

      seg_assign = SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["assign"]]
      n.ix = is.na(seg_assign)
      mut_cols = rep(NA, length(seg_assign) )
      mut_cols[!n.ix] = seg_clust_CCF_colors[ seg_assign[!n.ix] ]

      if( length(mut_cols) != 2*nrow(allele.segs)) { stop("wrong # of cols/segs") }

      PlotHscrAndSeghist( allele.segs, mut_cols, chr.arms.dat, max_CR, plot.hist=TRUE, plot.abs.fit=TRUE, comb=comb, mode.info=mode.info, Wq0=Wq0, comb.ab=comb.ab, fit.color=mode.colors[i], plot.seg.sem=TRUE )


    if(!is.null(segobj[["mut.cn.dat"]]) & !all(is.na(segobj[["mode.res"]][["modeled.muts"]][[i]][,"ccf_hat"])) )  ## protect against edge case of all muts on homozygously del SCNAs
    {
       mut.cn.dat <- segobj[["mut.cn.dat"]]
       modeled <- segobj[["mode.res"]][["modeled.muts"]][[i]]
       modeled.mut.dat <- cbind(mut.cn.dat, modeled)

     ## plot SSNVs on genome
       SSNV_cols = c("dodgerblue", "darkgrey", "seagreen3")   ## SC, clonal, mult>1
#       plot_SSNVs_on_genome( SSNV_model, SSNV_cols, modeled.mut.dat, segobj, i, mode.colors[i], max_SSNVs_plot=max_SSNVs_plot, verbose=verbose)

## plot
  ## Now plot SSNVs densities ... 2 plots 
        PlotSomaticMutDensities(modeled.mut.dat, segobj, i,
                                mode.colors[i], min.cov=3, max_SSNVs_plot=max_SSNVs_plot, verbose=verbose)
    }
    else { for(i in 1:2){ frame() } }

  }  ## modes
}

